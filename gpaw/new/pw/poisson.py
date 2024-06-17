from functools import cached_property
from math import pi

import numpy as np
from ase.units import Bohr, Ha
from gpaw.core import PWDesc, UGDesc
from gpaw.core.plane_waves import PWArray
from gpaw.new.poisson import PoissonSolver
from scipy.special import erf


def make_poisson_solver(pw: PWDesc,
                        grid: UGDesc,
                        charge: float,
                        strength: float = 1.0,
                        dipolelayer: bool = False,
                        **kwargs) -> PoissonSolver:
    if charge != 0.0 and not grid.pbc_c.any():
        return ChargedPWPoissonSolver(pw, grid, charge, strength, **kwargs)

    ps = PWPoissonSolver(pw, charge, strength)

    if dipolelayer:
        return DipoleLayerPWPoissonSolver(ps, grid, **kwargs)
    assert not kwargs
    return ps


class PWPoissonSolver(PoissonSolver):
    def __init__(self,
                 pw: PWDesc,
                 charge: float = 0.0,
                 strength: float = 1.0):
        self.pw = pw
        self.charge = charge
        self.strength = strength

        self.ekin_g = pw.ekin_G.copy()
        if pw.comm.rank == 0:
            # Avoid division by zero:
            self.ekin_g[0] = 1.0

    def __str__(self) -> str:
        txt = ('poisson solver:\n'
               f'  ecut: {self.pw.ecut * Ha}  # eV\n')
        if self.strength != 1.0:
            txt += f'  strength: {self.strength}\n'
        if self.charge != 0.0:
            txt += f'  uniform background charge: {self.charge}  # electrons\n'
        return txt

    def solve(self,
              vHt_g: PWArray,
              rhot_g: PWArray) -> float:
        """Solve Poisson equeation.

        Places result in vHt_g ndarray.
        """
        epot = self._solve(vHt_g, rhot_g)
        return epot

    def _solve(self,
               vHt_g,
               rhot_g) -> float:
        vHt_g.data[:] = 2 * pi * self.strength * rhot_g.data
        if self.pw.comm.rank == 0:
            # Use uniform backgroud charge in case we have a charged system:
            vHt_g.data[0] = 0.0
        if not isinstance(self.ekin_g, vHt_g.xp.ndarray):
            self.ekin_g = vHt_g.xp.array(self.ekin_g)
        vHt_g.data /= self.ekin_g
        epot = 0.5 * vHt_g.integrate(rhot_g)
        return epot


class ChargedPWPoissonSolver(PWPoissonSolver):
    def __init__(self,
                 pw: PWDesc,
                 grid: UGDesc,
                 charge: float,
                 strength: float = 1.0,
                 alpha: float = None,
                 eps: float = 1e-5):
        """Reciprocal-space Poisson solver for charged molecules.

        * Add a compensating Guassian-shaped charge to the density
          in order to make the total charge neutral (placed in the
          middle of the unit cell

        * Solve Poisson equation.

        * Correct potential so that it has the correct 1/r
          asymptotic behavior

        * Correct energy to remove the artificial interaction with
          the compensation charge

        Parameters
        ----------
        pw: PWDesc
        grid: UGDesc
        charge: float
        strength: float
        alpha: float
        eps: float

        Attributes
        ----------
        alpha : float
        charge_g : np.ndarray
            Guassian-shaped charge in reciprocal space
        potential_g : PWArray
             Potential in reciprocal space created by charge_g
        """
        super().__init__(pw, charge, strength)

        if alpha is None:
            # Shortest distance from center to edge of cell:
            rcut = 0.5 / (pw.icell**2).sum(axis=1).max()**0.5

            # Make sure e^(-alpha*rcut^2)=eps:
            alpha = -rcut**-2 * np.log(eps)

        self.alpha = alpha

        center_v = pw.cell_cv.sum(axis=0) / 2
        G2_g = 2 * pw.ekin_G
        G_gv = pw.G_plus_k_Gv
        self.charge_g = np.exp(-1 / (4 * alpha) * G2_g +
                               1j * (G_gv @ center_v))
        self.charge_g *= charge / pw.dv

        R_Rv = grid.xyz()
        d_R = ((R_Rv - center_v)**2).sum(axis=3)**0.5
        potential_R = grid.empty()

        # avoid division by 0
        zero_indx = d_R == 0
        d_R[zero_indx] = 1
        potential_R.data[:] = charge * erf(alpha**0.5 * d_R) / d_R
        # at zero we should have:
        # erf(alpha**0.5 * d_R) / d_R = alpha**0.5 * 2 / sqrt(pi)
        potential_R.data[zero_indx] = charge * alpha**0.5 * 2 / np.sqrt(pi)
        self.potential_g = potential_R.fft(pw=pw)

    def __str__(self) -> str:
        txt, x, _ = super().__str__().rsplit('\n', 2)
        assert x.startswith('  uniform background charge:')
        txt += (
            '\n  # using Gaussian-shaped compensation charge: e^(-alpha r^2)\n'
            f'  alpha: {self.alpha}   # bohr^-2')
        return txt

    def _solve(self,
               vHt_g,
               rhot_g) -> float:
        neutral_g = rhot_g.copy()
        neutral_g.data += self.charge_g

        if neutral_g.desc.comm.rank == 0:
            error = neutral_g.data[0]  # * self.pd.gd.dv
            assert error.imag == 0.0, error
            assert abs(error.real) < 0.00001, error
            neutral_g.data[0] = 0.0

        vHt_g.data[:] = 2 * pi * neutral_g.data
        vHt_g.data /= self.ekin_g
        epot = 0.5 * vHt_g.integrate(neutral_g)
        epot -= self.potential_g.integrate(rhot_g)
        epot -= self.charge**2 * (self.alpha / 2 / pi)**0.5
        vHt_g.data -= self.potential_g.data
        return epot


class DipoleLayerPWPoissonSolver(PoissonSolver):
    def __init__(self,
                 ps: PWPoissonSolver,
                 grid: UGDesc,
                 width: float = 1.0):  # Ångström
        self.ps = ps
        self.grid = grid
        self.width = width / Bohr
        (self.axis,) = np.where(~grid.pbc_c)[0]
        self.correction = np.nan

    def solve(self,
              vHt_g: PWArray,
              rhot_g: PWArray) -> float:
        epot = self.ps.solve(vHt_g, rhot_g)

        dip_v = -rhot_g.moment()
        c = self.axis
        L = self.grid.cell_cv[c, c]
        self.correction = 2 * np.pi * dip_v[c] * L / self.grid.volume
        vHt_g.data -= 2 * self.correction * self.sawtooth_g.data
        return epot + 2 * np.pi * dip_v[c]**2 / self.grid.volume

    def dipole_layer_correction(self) -> float:
        return self.correction

    @cached_property
    def sawtooth_g(self) -> PWArray:
        grid = self.grid
        if grid.comm.rank == 0:
            c = self.axis
            L = grid.cell_cv[c, c]
            w = self.width / 2
            assert w < L / 2, (w, L, c)
            gc = int(w / L * grid.size_c[c])
            x = np.linspace(0, L, grid.size_c[c], endpoint=False)
            sawtooth = x / L - 0.5
            a = 1 / L - 0.75 / w
            b = 0.25 / w**3
            sawtooth[:gc] = x[:gc] * (a + b * x[:gc]**2)
            sawtooth[-gc:] = -sawtooth[gc:0:-1]
            sawtooth_r = grid.new(comm=None).empty()
            shape = [1, 1, 1]
            shape[c] = -1
            sawtooth_r.data[:] = sawtooth.reshape(shape)
            sawtooth_g = sawtooth_r.fft(pw=self.ps.pw.new(comm=None)).data
        else:
            sawtooth_g = None

        result_g = self.ps.pw.empty()
        result_g.scatter_from(sawtooth_g)
        return result_g
