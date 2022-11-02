import os
import sys
import pickle
from math import pi

import numpy as np
from ase.units import Hartree, Bohr

import gpaw.mpi as mpi

from gpaw.response.coulomb_kernels import get_coulomb_kernel
from gpaw.response.wstc import WignerSeitzTruncatedCoulomb
from gpaw.response.density_kernels import get_density_xc_kernel
from gpaw.response.chi0 import Chi0Calculator, new_frequency_descriptor
from gpaw.response.pair import get_gs_and_context, NoCalculatorPairDensity


class DielectricFunctionCalculator:
    def __init__(self, chi0calc, truncation):
        from gpaw.response.pw_parallelization import Blocks1D
        self.chi0calc = chi0calc

        self.truncation = truncation
        self.context = chi0calc.context
        self.wd = chi0calc.wd
        self.blocks1d = Blocks1D(self.context.world, len(self.wd))

    @property
    def gs(self):
        return self.chi0calc.gs

    def calculate_chi0(self, q_c, spin='all'):
        """Calculates the response function.

        Calculate the response function for a specific momentum.

        q_c: [float, float, float]
            The momentum wavevector.
        spin : str or int
            If 'all' then include all spins.
            If 0 or 1, only include this specific spin.
            (not used in transverse reponse functions)
        """

        chi0 = self.chi0calc.calculate(q_c, spin)
        chi0_wGG = chi0.distribute_frequencies()
        self.context.write_timer()
        return chi0.pd, chi0_wGG, chi0.chi0_wxvG, chi0.chi0_wvv

    def read(self, name):
        self.context.print('Reading from', name)
        with open(name, 'rb') as fd:
            omega_w, pd, chi0_wGG, chi0_wxvG, chi0_wvv = pickle.load(fd)
            for omega in self.wd.omega_w:
                assert np.any(np.abs(omega - omega_w) < 1e-8)

            wmin = np.argmin(np.abs(np.min(self.wd.omega_w) - omega_w))
            world = self.context.world

            nw = len(omega_w)
            nG = pd.ngmax

            blocks1d = self.blocks1d

            mynw = blocks1d.blocksize

            if chi0_wGG is not None:
                # Old file format:
                chi0_wGG = chi0_wGG[wmin + blocks1d.a:blocks1d.b].copy()
            else:
                if world.rank == 0:
                    chi0_wGG = np.empty((mynw, nG, nG), complex)
                    for _ in range(wmin):
                        pickle.load(fd)
                    for chi0_GG in chi0_wGG:
                        chi0_GG[:] = pickle.load(fd)
                    tmp_wGG = np.empty((mynw, nG, nG), complex)
                    w1 = mynw
                    for rank in range(1, world.size):
                        w2 = min(w1 + mynw, nw)
                        for w in range(w2 - w1):
                            tmp_wGG[w] = pickle.load(fd)
                        world.send(tmp_wGG[:w2 - w1], rank)
                        w1 = w2
                else:
                    chi0_wGG = np.empty((self.blocks1d.nlocal, nG, nG),
                                        complex)
                    world.receive(chi0_wGG, 0)

            if chi0_wvv is not None:
                chi0_wxvG = chi0_wxvG[wmin:wmin + nw]
                chi0_wvv = chi0_wvv[wmin:wmin + nw]

        return pd, chi0_wGG, chi0_wxvG, chi0_wvv

    def collect(self, a_w):
        return self.blocks1d.collect(a_w)

    def get_frequencies(self):
        """ Return frequencies that Chi is evaluated on"""
        return self.wd.omega_w * Hartree

    def get_chi(self, xc='RPA', q_c=[0, 0, 0], spin='all',
                direction='x', return_VchiV=True, q_v=None,
                rshelmax=-1, rshewmin=None):
        """ Returns v^1/2 chi v^1/2 for the density response and chi for the
        spin response. The truncated Coulomb interaction is included as
        v^-1/2 v_t v^-1/2. This is in order to conform with
        the head and wings of chi0, which is treated specially for q=0.

        spin : str or int
            If 'all' then include all spins.
            If 0 or 1, only include this specific spin.
            (not used in transverse reponse functions)
        rshelmax : int or None
            Expand kernel in real spherical harmonics inside augmentation
            spheres. If None, the kernel will be calculated without
            augmentation. The value of rshelmax indicates the maximum index l
            to perform the expansion in (l < 6).
        rshewmin : float or None
            If None, the kernel correction will be fully expanded up to the
            chosen lmax. Given as a float, (0 < rshewmin < 1) indicates what
            coefficients to use in the expansion. If any coefficient
            contributes with less than a fraction of rshewmin on average,
            it will not be included.
        """
        pd, chi0_wGG, chi0_wxvG, chi0_wvv = self.calculate_chi0(q_c, spin)

        N_c = self.gs.kd.N_c

        Kbare_G = get_coulomb_kernel(pd,
                                     N_c,
                                     truncation=None,
                                     q_v=q_v)
        vsqr_G = Kbare_G**0.5
        nG = len(vsqr_G)

        if self.truncation is not None:
            if self.truncation == 'wigner-seitz':
                self.wstc = WignerSeitzTruncatedCoulomb(pd.gd.cell_cv, N_c)
            else:
                self.wstc = None
            Ktrunc_G = get_coulomb_kernel(pd,
                                          N_c,
                                          truncation=self.truncation,
                                          wstc=self.wstc,
                                          q_v=q_v)
            K_GG = np.diag(Ktrunc_G / Kbare_G)
        else:
            K_GG = np.eye(nG, dtype=complex)

        if pd.kd.gamma:
            if isinstance(direction, str):
                d_v = {'x': [1, 0, 0],
                       'y': [0, 1, 0],
                       'z': [0, 0, 1]}[direction]
            else:
                d_v = direction
            d_v = np.asarray(d_v) / np.linalg.norm(d_v)
            W = self.blocks1d.myslice
            chi0_wGG[:, 0] = np.dot(d_v, chi0_wxvG[W, 0])
            chi0_wGG[:, :, 0] = np.dot(d_v, chi0_wxvG[W, 1])
            chi0_wGG[:, 0, 0] = np.dot(d_v, np.dot(chi0_wvv[W], d_v).T)

        if xc != 'RPA':
            Kxc_GG = get_density_xc_kernel(pd,
                                           self.gs, self.context,
                                           functional=xc,
                                           chi0_wGG=chi0_wGG)
            K_GG += Kxc_GG / vsqr_G / vsqr_G[:, np.newaxis]

        # Invert Dyson eq.
        chi_wGG = []
        for chi0_GG in chi0_wGG:
            """v^1/2 chi0 V^1/2"""
            chi0_GG[:] = chi0_GG * vsqr_G * vsqr_G[:, np.newaxis]
            chi_GG = np.dot(np.linalg.inv(np.eye(nG) -
                                          np.dot(chi0_GG, K_GG)),
                            chi0_GG)
            if not return_VchiV:
                chi0_GG /= vsqr_G * vsqr_G[:, np.newaxis]
                chi_GG /= vsqr_G * vsqr_G[:, np.newaxis]
            chi_wGG.append(chi_GG)

        if len(chi_wGG):
            chi_wGG = np.array(chi_wGG)
        else:
            chi_wGG = np.zeros((0, nG, nG), complex)

        return pd, chi0_wGG, np.array(chi_wGG)

    def get_dynamic_susceptibility(self, xc='ALDA', q_c=[0, 0, 0],
                                   q_v=None,
                                   rshelmax=-1, rshewmin=None,
                                   filename='chiM_w.csv'):
        """Calculate the dynamic susceptibility.

        Returns macroscopic(could be generalized?) dynamic susceptibility:
        chiM0_w, chiM_w = DielectricFunction.get_dynamic_susceptibility()
        """

        pd, chi0_wGG, chi_wGG = self.get_chi(xc=xc, q_c=q_c,
                                             rshelmax=rshelmax,
                                             rshewmin=rshewmin,
                                             return_VchiV=False)

        rf0_w = np.zeros(len(chi_wGG), dtype=complex)
        rf_w = np.zeros(len(chi_wGG), dtype=complex)

        for w, (chi0_GG, chi_GG) in enumerate(zip(chi0_wGG, chi_wGG)):
            rf0_w[w] = chi0_GG[0, 0]
            rf_w[w] = chi_GG[0, 0]

        rf0_w = self.collect(rf0_w)
        rf_w = self.collect(rf_w)

        if filename is not None and mpi.rank == 0:
            write_response_function(filename, self.wd.omega_w * Hartree,
                                    rf0_w, rf_w)

        return rf0_w, rf_w

    def get_dielectric_matrix(self, xc='RPA', q_c=[0, 0, 0],
                              direction='x', symmetric=True,
                              calculate_chi=False, q_v=None,
                              add_intraband=False):
        r"""Returns the symmetrized dielectric matrix.

        ::

            \tilde\epsilon_GG' = v^{-1/2}_G \epsilon_GG' v^{1/2}_G',

        where::

            epsilon_GG' = 1 - v_G * P_GG' and P_GG'

        is the polarization.

        ::

            In RPA:   P = chi^0
            In TDDFT: P = (1 - chi^0 * f_xc)^{-1} chi^0

        in addition to RPA one can use the kernels, ALDA, Bootstrap and
        LRalpha (long-range kerne), where alpha is a user specified parameter
        (for example xc='LR0.25')

        The head of the inverse symmetrized dielectric matrix is equal
        to the head of the inverse dielectric matrix (inverse dielectric
        function)"""

        pd, chi0_wGG, chi0_wxvG, chi0_wvv = self.calculate_chi0(q_c)

        if add_intraband:
            print('add_intraband=True is not supported at this time')
            raise NotImplementedError

        N_c = self.gs.kd.N_c
        if self.truncation == 'wigner-seitz':
            self.wstc = WignerSeitzTruncatedCoulomb(pd.gd.cell_cv, N_c)
        else:
            self.wstc = None
        K_G = get_coulomb_kernel(pd,
                                 N_c,
                                 truncation=self.truncation,
                                 wstc=self.wstc,
                                 q_v=q_v)**0.5
        nG = len(K_G)

        if pd.kd.gamma:
            if isinstance(direction, str):
                d_v = {'x': [1, 0, 0],
                       'y': [0, 1, 0],
                       'z': [0, 0, 1]}[direction]
            else:
                d_v = direction

            d_v = np.asarray(d_v) / np.linalg.norm(d_v)
            W = self.blocks1d.myslice
            chi0_wGG[:, 0] = np.dot(d_v, chi0_wxvG[W, 0])
            chi0_wGG[:, :, 0] = np.dot(d_v, chi0_wxvG[W, 1])
            chi0_wGG[:, 0, 0] = np.dot(d_v, np.dot(chi0_wvv[W], d_v).T)
            if q_v is not None:
                print('Restoring q dependence of head and wings of chi0')
                chi0_wGG[:, 1:, 0] *= np.dot(q_v, d_v)
                chi0_wGG[:, 0, 1:] *= np.dot(q_v, d_v)
                chi0_wGG[:, 0, 0] *= np.dot(q_v, d_v)**2

        if xc != 'RPA':
            Kxc_GG = get_density_xc_kernel(pd,
                                           self.gs, self.context,
                                           functional=xc,
                                           chi0_wGG=chi0_wGG)

        if calculate_chi:
            chi_wGG = []

        for chi0_GG in chi0_wGG:
            if xc == 'RPA':
                P_GG = chi0_GG
            else:
                P_GG = np.dot(np.linalg.inv(np.eye(nG) -
                                            np.dot(chi0_GG, Kxc_GG)),
                              chi0_GG)
            if symmetric:
                e_GG = np.eye(nG) - P_GG * K_G * K_G[:, np.newaxis]
            else:
                K_GG = (K_G**2 * np.ones([nG, nG])).T
                e_GG = np.eye(nG) - P_GG * K_GG
            if calculate_chi:
                K_GG = np.diag(K_G**2)
                if xc != 'RPA':
                    K_GG += Kxc_GG
                chi_wGG.append(np.dot(np.linalg.inv(np.eye(nG) -
                                                    np.dot(chi0_GG, K_GG)),
                                      chi0_GG))
            chi0_GG[:] = e_GG

        # chi0_wGG is now the dielectric matrix
        if calculate_chi:
            if len(chi_wGG):
                chi_wGG = np.array(chi_wGG)
            else:
                chi_wGG = np.zeros((0, nG, nG), complex)

        if not calculate_chi:
            return chi0_wGG
        else:
            # chi_wGG is the full density response function..
            return pd, chi0_wGG, chi_wGG

    def get_dielectric_function(self, xc='RPA', q_c=[0, 0, 0], q_v=None,
                                direction='x', filename='df.csv'):
        """Calculate the dielectric function.

        Returns dielectric function without and with local field correction:
        df_NLFC_w, df_LFC_w = DielectricFunction.get_dielectric_function()
        """
        e_wGG = self.get_dielectric_matrix(xc, q_c, direction, q_v=q_v)
        df_NLFC_w = np.zeros(len(e_wGG), dtype=complex)
        df_LFC_w = np.zeros(len(e_wGG), dtype=complex)

        for w, e_GG in enumerate(e_wGG):
            df_NLFC_w[w] = e_GG[0, 0]
            df_LFC_w[w] = 1 / np.linalg.inv(e_GG)[0, 0]

        df_NLFC_w = self.collect(df_NLFC_w)
        df_LFC_w = self.collect(df_LFC_w)

        if filename is not None and mpi.rank == 0:
            write_response_function(filename, self.wd.omega_w * Hartree,
                                    df_NLFC_w, df_LFC_w)

        return df_NLFC_w, df_LFC_w

    def get_macroscopic_dielectric_constant(self, xc='RPA',
                                            direction='x', q_v=None):
        """Calculate macroscopic dielectric constant.

        Returns eM_NLFC and eM_LFC.

        Macroscopic dielectric constant is defined as the real part
        of dielectric function at w=0.

        Parameters:

        eM_LFC: float
            Dielectric constant without local field correction. (RPA, ALDA)
        eM2_NLFC: float
            Dielectric constant with local field correction.
        """

        self.context.print('', flush=False)
        self.context.print('%s Macroscopic Dielectric Constant:' % xc)

        df_NLFC_w, df_LFC_w = self.get_dielectric_function(
            xc=xc,
            filename=None,
            direction=direction,
            q_v=q_v)
        eps0 = np.real(df_NLFC_w[0])
        eps = np.real(df_LFC_w[0])
        self.context.print('  %s direction' % direction, flush=False)
        self.context.print('    Without local field: %f' % eps0, flush=False)
        self.context.print('    Include local field: %f' % eps)

        return eps0, eps

    def get_eels_spectrum(self, xc='RPA', q_c=[0, 0, 0],
                          direction='x', filename='eels.csv'):
        r"""Calculate EELS spectrum. By default, generate a file 'eels.csv'.

        EELS spectrum is obtained from the imaginary part of the
        density response function as, EELS(\omega) = - 4 * \pi / q^2 Im \chi.
        Returns EELS spectrum without and with local field corrections:

        df_NLFC_w, df_LFC_w = DielectricFunction.get_eels_spectrum()
        """

        # Calculate V^1/2 \chi V^1/2
        pd, Vchi0_wGG, Vchi_wGG = self.get_chi(xc=xc, q_c=q_c,
                                               direction=direction)

        # Calculate eels = -Im 4 \pi / q^2  \chi
        eels_NLFC_w = -(1. / (1. - Vchi0_wGG[:, 0, 0])).imag
        eels_LFC_w = - (Vchi_wGG[:, 0, 0]).imag

        # Collect frequencies
        eels_NLFC_w = self.collect(eels_NLFC_w)
        eels_LFC_w = self.collect(eels_LFC_w)

        # Write to file
        if filename is not None and self.context.world.rank == 0:
            omega_w = self.wd.omega_w
            write_response_function(filename, omega_w * Hartree,
                                    eels_NLFC_w, eels_LFC_w)

        return eels_NLFC_w, eels_LFC_w

    def get_polarizability(self, xc='RPA', direction='x', q_c=[0, 0, 0],
                           filename='polarizability.csv'):
        r"""Calculate the polarizability alpha.
        In 3D the imaginary part of the polarizability is related to the
        dielectric function by Im(eps_M) = 4 pi * Im(alpha). In systems
        with reduced dimensionality the converged value of alpha is
        independent of the cell volume. This is not the case for eps_M,
        which is ill-defined. A truncated Coulomb kernel will always give
        eps_M = 1.0, whereas the polarizability maintains its structure.

        By default, generate a file 'polarizability.csv'. The five columns are:
        frequency (eV), Real(alpha0), Imag(alpha0), Real(alpha), Imag(alpha)
        alpha0 is the result without local field effects and the
        dimension of alpha is \AA to the power of non-periodic directions
        """

        cell_cv = self.gs.gd.cell_cv
        pbc_c = self.gs.pbc

        if pbc_c.all():
            V = 1.0
        else:
            V = np.abs(np.linalg.det(cell_cv[~pbc_c][:, ~pbc_c]))

        if not self.truncation:
            """Standard expression for the polarizability"""
            df0_w, df_w = self.get_dielectric_function(xc=xc,
                                                       q_c=q_c,
                                                       filename=None,
                                                       direction=direction)
            alpha_w = V * (df_w - 1.0) / (4 * pi)
            alpha0_w = V * (df0_w - 1.0) / (4 * pi)
        else:
            # Since eps_M = 1.0 for a truncated Coulomb interaction, it does
            # not make sense to apply it here. Instead one should define the
            # polarizability by
            #
            #     alpha * eps_M^{-1} = -L / (4 * pi) * <v_ind>
            #
            # where <v_ind> = 4 * pi * \chi / q^2 is the averaged induced
            # potential (relative to the strength of the  external potential).
            # With the bare Coulomb potential, this expression is equivalent to
            # the standard one. In a 2D system \chi should be calculated with a
            # truncated Coulomb potential and eps_M = 1.0

            self.context.print('Using truncated Coulomb interaction')

            pd, chi0_wGG, chi_wGG = self.get_chi(xc=xc,
                                                 q_c=q_c,
                                                 direction=direction)
            alpha_w = -V * (chi_wGG[:, 0, 0]) / (4 * pi)
            alpha0_w = -V * (chi0_wGG[:, 0, 0]) / (4 * pi)

            alpha_w = self.collect(alpha_w)
            alpha0_w = self.collect(alpha0_w)

        # Convert to external units
        hypervol = Bohr**(sum(~pbc_c))
        alpha0_w *= hypervol
        alpha_w *= hypervol

        # Write results file
        if filename is not None and self.context.world.rank == 0:
            omega_w = self.wd.omega_w
            write_response_function(filename, omega_w * Hartree,
                                    alpha0_w, alpha_w)

        return alpha0_w, alpha_w

    def check_sum_rule(self, spectrum=None):
        """Check f-sum rule.

        It takes the y of a spectrum as an entry and it checks its integral.

        spectrum: np.ndarray
            Input spectrum

        Note: not tested for spin response
        """

        assert (self.wd.omega_w[1:] - self.wd.omega_w[:-1]).ptp() < 1e-10

        if spectrum is None:
            raise ValueError('No spectrum input ')
        dw = self.wd.omega_w[1] - self.wd.omega_w[0]
        N1 = 0
        for iw in range(len(spectrum)):
            w = iw * dw
            N1 += spectrum[iw] * w
        N1 *= dw * self.gs.volume / (2 * pi**2)

        self.context.print('', flush=False)
        self.context.print('Sum rule:', flush=False)
        nv = self.gs.nvalence
        self.context.print('N1 = %f, %f  %% error' %
                           (N1, (N1 - nv) / nv * 100))


class DielectricFunction(DielectricFunctionCalculator):
    """This class defines dielectric function related physical quantities."""

    def __init__(self, calc, *,
                 frequencies=None,
                 domega0=None,  # deprecated
                 omega2=None,  # deprecated
                 omegamax=None,  # deprecated
                 ecut=50,
                 hilbert=True,
                 nbands=None, eta=0.2, ftol=1e-6, threshold=1,
                 intraband=True, nblocks=1, world=mpi.world, txt=sys.stdout,
                 truncation=None, disable_point_group=False,
                 disable_time_reversal=False,
                 integrationmode=None, rate=0.0,
                 eshift=0.0):
        """Creates a DielectricFunction object.

        calc: str
            The groundstate calculation file that the linear response
            calculation is based on.
        frequencies:
            Input parameters for frequency_grid.
            Can be array of frequencies to evaluate the response function at
            or dictionary of paramaters for build-in nonlinear grid
            (see :ref:`frequency grid`).
        ecut: float
            Plane-wave cut-off.
        hilbert: bool
            Use hilbert transform.
        nbands: int
            Number of bands from calc.
        eta: float
            Broadening parameter.
        ftol: float
            Threshold for including close to equally occupied orbitals,
            f_ik - f_jk > ftol.
        threshold: float
            Threshold for matrix elements in optical response perturbation
            theory.
        intraband: bool
            Include intraband transitions.
        world: comm
            mpi communicator.
        nblocks: int
            Split matrices in nblocks blocks and distribute them G-vectors or
            frequencies over processes.
        txt: str
            Output file.
        truncation: str
            'wigner-seitz' for Wigner Seitz truncated Coulomb.
            '2D, 1D or 0d for standard analytical truncation schemes.
            Non-periodic directions are determined from k-point grid
        eshift: float
            Shift unoccupied bands
        """

        gs, context = get_gs_and_context(calc, txt, world, timer=None)
        nbands = nbands or gs.bd.nbands

        wd = new_frequency_descriptor(gs, context, nbands, frequencies,
                                      domega0=domega0,
                                      omega2=omega2, omegamax=omegamax)

        pair = NoCalculatorPairDensity(
            gs=gs, context=context, threshold=threshold, nblocks=nblocks)

        chi0calc = Chi0Calculator(
            wd=wd,
            pair=pair,
            ecut=ecut, nbands=nbands, eta=eta,
            hilbert=hilbert,
            ftol=ftol,
            intraband=intraband,
            disable_point_group=disable_point_group,
            disable_time_reversal=disable_time_reversal,
            integrationmode=integrationmode,
            rate=rate, eshift=eshift
        )

        super().__init__(chi0calc=chi0calc, truncation=truncation)


def write_response_function(filename, omega_w, rf0_w, rf_w):
    with open(filename, 'w') as fd:
        for omega, rf0, rf in zip(omega_w, rf0_w, rf_w):
            if rf0_w.dtype == complex:
                print('%.6f, %.6f, %.6f, %.6f, %.6f' %
                      (omega, rf0.real, rf0.imag, rf.real, rf.imag),
                      file=fd)
            else:
                print('%.6f, %.6f, %.6f' % (omega, rf0, rf), file=fd)


def read_response_function(filename):
    """Read a stored response function file"""
    d = np.loadtxt(filename, delimiter=',')
    omega_w = np.array(d[:, 0], float)

    if d.shape[1] == 3:
        # Real response function
        rf0_w = np.array(d[:, 1], float)
        rf_w = np.array(d[:, 2], float)
    elif d.shape[1] == 5:
        rf0_w = np.array(d[:, 1], complex)
        rf0_w.imag = d[:, 2]
        rf_w = np.array(d[:, 3], complex)
        rf_w.imag = d[:, 4]
    else:
        raise ValueError(f'Unexpected array dimension {d.shape}')

    return omega_w, rf0_w, rf_w
