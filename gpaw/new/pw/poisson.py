from math import pi

import numpy as np
from ase.units import Ha
from gpaw.core import PlaneWaves, UniformGrid
from gpaw.core.plane_waves import PlaneWaveExpansions
from gpaw.new import cached_property
from gpaw.new.poisson import PoissonSolver


def make_poisson_solver(pw: PlaneWaves,
                        grid: UniformGrid,
                        charge: float,
                        strength: float = 1.0,
                        dipolelayer: bool = False,
                        **kwargs) -> PoissonSolver:
    if charge != 0.0 and not grid.pbc_c.any():
        return ChargedPWPoissonSolver(pw, grid, charge, strength, **kwargs)

    ps = PWPoissonSolver(pw, charge, strength)

    if dipolelayer:
        ps = DipoleLayerPWPoissonSolver(ps, grid, **kwargs)
    else:
        assert not kwargs

    return ps


class PWPoissonSolver(PoissonSolver):
    def __init__(self,
                 pw: PlaneWaves,
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
            txt += f'  strength: {self.strength}'
        if self.charge != 0.0:
            txt += f'  uniform background charge: {self.charge}  # electrons'
        return txt

    def solve(self,
              vHt_g: PlaneWaveExpansions,
              rhot_g: PlaneWaveExpansions) -> float:
        """Solve Poisson equeation.

        Places result in vHt_g ndarray.
        """
        epot = self._solve(vHt_g, rhot_g)
        return epot

    def _solve(self,
               vHt_g,
               rhot_g) -> float:
        vHt_g.data[:] = 2 * pi * self.stgrength * rhot_g.data
        if self.pw.comm.rank == 0:
            # Use uniform backgroud charge in case we have a charged system:
            vHt_g.data[0] = 0.0
        vHt_g.data /= self.ekin_g
        epot = 0.5 * vHt_g.integrate(rhot_g)
        return epot


class ChargedPWPoissonSolver(PWPoissonSolver):
    def __init__(self,
                 pw: PlaneWaves,
                 grid: UniformGrid,
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
        potential_R = charge * np.erf(alpha**0.5 * d_R) / d_R
        self.potential_g = potential_R.fft(pw=pw)

    def __str__(self) -> str:
        txt, x = str(super()).rsplit('\n', 1)
        assert x.startswith('  uniform background charge:')
        txt += (
            '  # using Gaussian-shaped compensation charge: e^(-alpha r^2)\n'
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
                 grid: UniformGrid):
        self.ps = ps
        self.grid = grid
        (self.axis,) = np.where(~grid.pbc_c)

    def solve(self,
              vHt_g: PlaneWaveExpansions,
              rhot_g: PlaneWaveExpansions) -> float:
        if self.sawtooth_q is None:
            self.initialize_sawtooth()

        epot = self.ps.solve(vHt_g, rhot_g)

        dip_v = -rhot_g.moment()
        c = self.axis
        L = self.grid.cell_cv[c, c]
        self.correction = 2 * np.pi * dip_v[c] * L / self.grid.volume
        vHt_g.data -= 2 * self.correction * self.sawtooth_g.data

        return epot + 2 * np.pi * dip_v[c]**2 / self.grid.volume

    @cached_property
    def sawtooth_g(self) -> PlaneWaveExpansions:
        grid = self.grid
        if grid.comm.rank == 0:
            c = self.axis
            L = grid.cell_cv[c, c]
            w = self.width / 2
            assert w < L / 2
            gc = int(w / grid.size_c[c])
            x = np.linspace(0, L, grid.size_c[c], endpoint=False)
            sawtooth = x / L - 0.5
            a = 1 / L - 0.75 / w
            b = 0.25 / w**3
            sawtooth[:gc] = x[:gc] * (a + b * x[:gc]**2)
            sawtooth[-gc:] = -sawtooth[gc:0:-1]
            sawtooth_r = grid.new(comm=None).empty()
            shape = [1, 1, 1]
            shape[c] = -1
            sawtooth_r[:] = sawtooth.reshape(shape)
            sawtooth_g = sawtooth_r.fft(pw=self.ps.pw).data
        else:
            sawtooth_g = None

        return self.pw.empty.scatter_from(sawtooth_g)
