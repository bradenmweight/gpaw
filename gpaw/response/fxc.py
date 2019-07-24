# -*- coding: utf-8
"""This module calculates XC kernels for response function calculations.
"""
from __future__ import print_function

import numpy as np

from ase.utils.timing import Timer
from ase.units import Ha
'''
from scipy.special import spherical_jn

import gpaw.mpi as mpi
from gpaw.xc import XC
from gpaw.sphere.lebedev import weight_n, R_nv
'''
from gpaw.io.tar import Reader
'''
from gpaw.spherical_harmonics import Yarr

from ase.utils import convert_string_to_fd
from ase.utils.timing import Timer, timer
from ase.units import Bohr, Ha
'''

from gpaw.response.kxc import AdiabaticSusceptibilityFXC
from gpaw.response.tms import find_goldstone_scaling


'''
def get_fxc(fxc, calc, response='susceptibility', mode='pw',
            world=mpi.world, txt='-', timer=None, **kwargs):
    """Factory function getting an initiated version of the fxc class."""
    functional = fxc
    fxc = create_fxc(functional, response, mode)
    return fxc(functional, calc, world=world, txt=txt, timer=timer, **kwargs)


def create_fxc(functional, response, mode):
    """Creator component for the FXC classes."""
    # Only one kind of response and mode is supported for now
    if functional in ['ALDA_x', 'ALDA_X', 'ALDA']:
        if response == 'susceptibility' and mode == 'pw':
            return AdiabaticSusceptibilityFXC
    raise ValueError(functional, response, mode)


class FXC:
    """General class to calculate exchange-correlation kernels."""

    def __init__(self, functional, calc, world=mpi.world, txt='-', timer=None):
        """
        Parameters
        ----------
        functional : str
            xc-functional
        kslrf : KohnShamLinearResponseFunction instance
        calc : obj
            GPAW calculator object of ground state calculation
        world : mpi.world
        txt : str or filehandle
            defines output file through ase.utils.convert_string_to_fd
        """
        self.functional = functional
        self.calc = calc

        self.world = world
        self.fd = convert_string_to_fd(txt, world)

        # Timer
        if timer is None:
            self.timer = Timer()
            self.write_timer = True
        else:
            self.timer = timer
            self.write_timer = False

    def calculate(self, *args, **kwargs):
        raise NotImplementedError

    def is_calculated(self, *args, **kwargs):
        # Read/write has not been implemented
        return False

    def read(self, *args, **kwargs):
        raise NotImplementedError

    def write(self, *args, **kwargs):
        # Not implemented
        pass


class AdiabaticSusceptibilityFXC(FXC):
    """Class for calculating adiabatic exchange correlation kernels for
    susceptibility calculations."""

    def __init__(self, functional, calc, world=mpi.world, txt='-', timer=None,
                 **kwargs):
        FXC.__init__(self, functional, calc, world=world, txt=txt, timer=timer)

        # The two calculators might be merged into one?
        # It would be much easier, if all spincomponents could be evaluated
        # by a single calculator class also supporting Goldstone scaling.
        # For now, we will use the existing code and make a hack.
        self.dcalculator = AdiabaticDensityKernelCalculator(world=self.world,
                                                            txt=self.fd,
                                                            timer=self.timer,
                                                            **kwargs)

        self.tcalculator = ATKCHack(world=self.world, txt=self.fd,
                                    timer=self.timer, **kwargs)

    def calculate(self, spincomponent, pd, kslrf=None, chiks_wGG=None):
        """Goldstone scaling is not yet implemented."""
        calculator = self.create_calculator(spincomponent)
        Kxc_GG = calculator(self.calc, pd, self.functional,
                            kslrf=kslrf, chiks_wGG=chiks_wGG)  # Goldstone stuff should be rewritten XXX

        if self.write_timer:
            self.timer.write(self.fd)

        return Kxc_GG

    def create_calculator(self, spincomponent):
        """Creator component."""
        if spincomponent in ['00', 'uu', 'dd']:
            return self.dcalculator
        elif spincomponent in ['+-', '-+']:
            return self.tcalculator
        else:
            raise ValueError(spincomponent)
'''


def get_xc_kernel(pd, chi0, functional='ALDA', kernel='density',
                  rshe=0.99,
                  chi0_wGG=None,
                  fxc_scaling=None,
                  density_cut=None,
                  spinpol_cut=None):
    """
    Factory function that calls the relevant functions below
    """

    if kernel == 'density':
        return get_density_xc_kernel(pd, chi0, functional=functional,
                                     rshe=rshe,
                                     chi0_wGG=chi0_wGG,
                                     density_cut=density_cut)
    elif kernel in ['+-', '-+']:
        # Currently only collinear adiabatic xc kernels are implemented
        # for which the +- and -+ kernels are the same
        return get_transverse_xc_kernel(pd, chi0, functional=functional,
                                        rshe=rshe,
                                        chi0_wGG=chi0_wGG,
                                        fxc_scaling=fxc_scaling,
                                        density_cut=density_cut,
                                        spinpol_cut=spinpol_cut)
    else:
        raise ValueError('%s kernels not implemented' % kernel)


def get_density_xc_kernel(pd, chi0, functional='ALDA',
                          rshe=0.99,
                          chi0_wGG=None,
                          density_cut=None):
    """
    Density-density xc kernels
    Factory function that calls the relevant functions below
    """

    calc = chi0.calc
    fd = chi0.fd
    nspins = len(calc.density.nt_sG)
    assert nspins == 1

    if functional[0] == 'A':
        # Standard adiabatic kernel
        print('Calculating %s kernel' % functional, file=fd)
        Kcalc = AdiabaticSusceptibilityFXC(calc, functional,
                                           world=chi0.world, txt=fd,
                                           timer=chi0.timer,
                                           rshe=rshe,
                                           density_cut=density_cut)
        Kxc_GG = Kcalc('00', pd)
        if pd.kd.gamma:
            Kxc_GG[0, :] = 0.0
            Kxc_GG[:, 0] = 0.0
        Kxc_sGG = np.array([Kxc_GG])
    elif functional[0] == 'r':
        # Renormalized kernel
        print('Calculating %s kernel' % functional, file=fd)
        Kxc_sGG = calculate_renormalized_kernel(pd, calc, functional, fd)
    elif functional[:2] == 'LR':
        print('Calculating LR kernel with alpha = %s' % functional[2:],
              file=fd)
        Kxc_sGG = calculate_lr_kernel(pd, calc, alpha=float(functional[2:]))
    elif functional == 'DM':
        print('Calculating DM kernel', file=fd)
        Kxc_sGG = calculate_dm_kernel(pd, calc)
    elif functional == 'Bootstrap':
        print('Calculating Bootstrap kernel', file=fd)
        Kxc_sGG = get_bootstrap_kernel(pd, chi0, chi0_wGG, fd)
    else:
        raise ValueError('density-density %s kernel not'
                         + ' implemented' % functional)

    return Kxc_sGG[0]


def get_transverse_xc_kernel(pd, chi0, functional='ALDA_x',
                             rshe=0.99,
                             chi0_wGG=None,
                             fxc_scaling=None,
                             density_cut=None,
                             spinpol_cut=None):
    """ +-/-+ xc kernels
    Currently only collinear ALDA kernels are implemented
    Factory function that calls the relevant functions below
    """

    calc = chi0.calc
    fd = chi0.fd
    nspins = len(calc.density.nt_sG)
    assert nspins == 2

    if functional in ['ALDA_x', 'ALDA_X', 'ALDA']:
        # Adiabatic kernel
        print("Calculating transverse %s kernel" % functional, file=fd)
        Kcalc = AdiabaticSusceptibilityFXC(calc, functional,
                                           world=chi0.world, txt=fd,
                                           timer=chi0.timer,
                                           rshe=rshe,
                                           density_cut=density_cut,
                                           spinpol_cut=spinpol_cut)
    else:
        raise ValueError("%s spin kernel not implemented" % functional)

    Kxc_GG = Kcalc('+-', pd)

    if fxc_scaling is not None:
        assert isinstance(fxc_scaling[0], bool)
        if fxc_scaling[0]:
            if fxc_scaling[1] is None:
                assert pd.kd.gamma
                print('Finding rescaling of kernel to fulfill the '
                      'Goldstone theorem', file=chi0.fd)
                fxc_scaling[1] = find_goldstone_scaling(chi0.omega_w,
                                                        -chi0_wGG, Kxc_GG,
                                                        world=chi0.world)

            assert isinstance(fxc_scaling[1], float)
            Kxc_GG *= fxc_scaling[1]

    return Kxc_GG


'''
def find_Goldstone_scaling(pd, chi0, chi0_wGG, Kxc_GG):
    """ Find a scaling of the kernel to move the magnon peak to omeaga=0. """
    # q should be gamma - scale to hit Goldstone
    assert pd.kd.gamma

    fd = chi0.fd
    omega_w = chi0.omega_w
    wgs = np.abs(omega_w).argmin()

    if not np.allclose(omega_w[wgs], 0., rtol=1.e-8):
        raise ValueError("Frequency grid needs to include"
                         + " omega=0. to allow Goldstone scaling")

    fxcs = 1.
    print("Finding rescaling of kernel to fulfill the Goldstone theorem",
          file=fd)

    world = chi0.world
    # Only one rank, rgs, has omega=0 and finds rescaling
    nw = len(omega_w)
    mynw = (nw + world.size - 1) // world.size
    rgs, mywgs = wgs // mynw, wgs % mynw
    fxcsbuf = np.empty(1, dtype=float)
    if world.rank == rgs:
        chi0_GG = chi0_wGG[mywgs]
        chi_GG = np.dot(np.linalg.inv(np.eye(len(chi0_GG)) -
                                      np.dot(chi0_GG, Kxc_GG * fxcs)),
                        chi0_GG)
        # Scale so that kappaM=0 in the static limit (omega=0)
        kappaM = (chi0_GG[0, 0] / chi_GG[0, 0]).real
        # If kappaM > 0, increase scaling (recall: kappaM ~ 1 - Kxc Re{chi_0})
        scaling_incr = 0.1 * np.sign(kappaM)
        while abs(kappaM) > 1.e-7 and abs(scaling_incr) > 1.e-7:
            fxcs += scaling_incr
            if fxcs <= 0.0 or fxcs >= 10.:
                raise Exception('Found an invalid fxc_scaling of %.4f' % fxcs)

            chi_GG = np.dot(np.linalg.inv(np.eye(len(chi0_GG)) -
                                          np.dot(chi0_GG, Kxc_GG * fxcs)),
                            chi0_GG)
            kappaM = (chi0_GG[0, 0] / chi_GG[0, 0]).real

            # If kappaM changes sign, change sign and refine increment
            if np.sign(kappaM) != np.sign(scaling_incr):
                scaling_incr *= -0.2
        fxcsbuf[:] = fxcs

    # Broadcast found rescaling
    world.broadcast(fxcsbuf, rgs)
    fxcs = fxcsbuf[0]

    return fxcs
'''


'''
class AdiabaticKernelCalculator:
    """ Adiabatic kernels with PAW """

    def __init__(self, world=mpi.world, txt='-', timer=None,
                 rshe=0.99, ecut=None, **unused):
        """
        rshe : float or None
            Expand kernel in real spherical harmonics inside augmentation
            spheres. If None, the kernel will be calculated without
            augmentation. The value of rshe (0<rshe<1) sets a convergence
            criteria for the expansion in real spherical harmonics.
        """
        self.world = world
        self.fd = convert_string_to_fd(txt, world)

        # Timer
        if timer is None:
            self.timer = Timer()
            self.write_timer = True
        else:
            self.timer = timer
            self.write_timer = False

        self.ecut = ecut

        self.permitted_functionals = []
        self.functional = None

        self.rshe = rshe is not None

        if self.rshe:
            self.rshecc = rshe
            self.dfSns_g = None
            self.dfSns_gL = None
            self.dfmask_g = None
            self.rsheconvmin = None

            self.rsheL_M = None

    @timer('Calculate XC kernel')
    def __call__(self, calc, pd, functional, **unused):
        print('', file=self.fd)
        print('Excecuting a response.fxc.AdiabaticKernelCalculator.__call__()',
              file=self.fd)

        assert functional in self.permitted_functionals
        self.functional = functional
        
        # Get the spin density we need and allocate fxc
        n_sG = self.get_density_on_grid(calc, pd)
        fxc_G = np.zeros(np.shape(n_sG[0]))

        print('    Calculating fxc on real space grid', file=self.fd)
        self.add_fxc(pd.gd, n_sG, fxc_G)

        # Fourier transform to reciprocal space
        Kxc_GG = self.ft_from_grid(fxc_G, pd)

        if self.rshe:  # Do PAW correction to Fourier transformed kernel
            KxcPAW_GG = self.calculate_kernel_paw_correction(calc, pd)
            Kxc_GG += KxcPAW_GG

        if self.write_timer:
            self.timer.write(self.fd)

        return Kxc_GG / pd.gd.volume
            
    def get_density_on_grid(self, calc, pd):
        """Get the spin density on coarse real-space grid.
        
        Returns
        -------
        nt_sG or n_sG : nd.array
            Spin density on coarse real-space grid. If not self.rshe, use
            the PAW corrected all-electron spin density.
        """
        if self.rshe:
            return calc.density.nt_sG  # smooth density
        
        print('    Calculating all-electron density', file=self.fd)
        with self.timer('Calculating all-electron density'):
            return calc.density.get_all_electron_density(atoms=calc.atoms,
                                                         gridrefinement=1)[0]

    @timer('Fourier transform of kernel from real-space grid')
    def ft_from_grid(self, fxc_G, pd):
        print('    Fourier transforming kernel from real-space grid',
              file=self.fd)
        nG = pd.gd.N_c
        nG0 = nG[0] * nG[1] * nG[2]

        tmp_g = np.fft.fftn(fxc_G) * pd.gd.volume / nG0

        # The unfolding procedure could use vectorization and parallelization.
        # This remains a slow step for now. XXX
        Kxc_GG = np.zeros((pd.ngmax, pd.ngmax), dtype=complex)
        for iG, iQ in enumerate(pd.Q_qG[0]):
            iQ_c = (np.unravel_index(iQ, nG) + nG // 2) % nG - nG // 2
            for jG, jQ in enumerate(pd.Q_qG[0]):
                jQ_c = (np.unravel_index(jQ, nG) + nG // 2) % nG - nG // 2
                ijQ_c = (iQ_c - jQ_c)
                if (abs(ijQ_c) < nG // 2).all():
                    Kxc_GG[iG, jG] = tmp_g[tuple(ijQ_c)]

        return Kxc_GG

    @timer('Calculate PAW corrections to kernel')
    def calculate_kernel_paw_correction(self, calc, pd):
        print("    Calculating PAW corrections to the kernel",
              file=self.fd)

        # Allocate array and distribute plane waves
        npw = pd.ngmax
        KxcPAW_GG = np.zeros((npw, npw), dtype=complex)
        G_myG = self._distribute_correction(npw)

        # Calculate (G-G') reciprocal space vectors, their length and direction
        dG_myGGv, dG_myGG, dGn_myGGv = self._calculate_dG(pd, G_myG)

        # Calculate PAW correction to each augmentation sphere (to each atom)
        R_av = calc.atoms.positions / Bohr
        for a, R_v in enumerate(R_av):
            # Calculate dfxc on Lebedev quadrature and radial grid
            # Please note: Using the radial grid descriptor with add_fxc
            # might give problems beyond ALDA
            df_ng, Y_nL, rgd = self._calculate_dfxc(calc, a)

            # Calculate the surface norm square of df
            self.dfSns_g = self._ang_int(df_ng ** 2)
            # Reduce radial grid by excluding points where dfSns_g = 0
            df_ng, r_g, dv_g = self._reduce_radial_grid(df_ng, rgd)

            # Expand correction in real spherical harmonics
            df_gL = self._perform_rshe(a, df_ng, Y_nL)
            # Reduce expansion by removing coefficients that are zero
            df_gM, L_M, l_M = self._reduce_rsh_expansion(df_gL)

            # Expand plane wave differences (G-G')
            (ii_MmyGG,
             j_gMmyGG,
             Y_MmyGG) = self._expand_plane_waves(dG_myGG, dGn_myGGv,
                                                 r_g, L_M, l_M)

            # Perform integration
            with self.timer('Integrate PAW correction'):
                coefatomR_GG = np.exp(-1j * np.inner(dG_myGGv, R_v))
                coefatomang_MGG = ii_MmyGG * Y_MmyGG
                coefatomrad_MGG = np.tensordot(j_gMmyGG * df_gL[:, L_M,
                                                                np.newaxis,
                                                                np.newaxis],
                                               dv_g, axes=([0, 0]))
                coefatom_GG = np.sum(coefatomang_MGG * coefatomrad_MGG, axis=0)
                KxcPAW_GG[G_myG] += coefatom_GG * coefatomR_GG

        self.world.sum(KxcPAW_GG)

        return KxcPAW_GG

    def _distribute_correction(self, npw):
        """Distribute correction"""
        nGpr = (npw + self.world.size - 1) // self.world.size
        Ga = min(self.world.rank * nGpr, npw)
        Gb = min(Ga + nGpr, npw)

        return range(Ga, Gb)

    def _calculate_dG(self, pd, G_myG):
        """Calculate (G-G') reciprocal space vectors,
        their length and direction"""
        npw = pd.ngmax
        G_Gv = pd.get_reciprocal_vectors()

        # Calculate bare dG
        dG_myGGv = np.zeros((len(G_myG), npw, 3))
        for v in range(3):
            dG_myGGv[:, :, v] = np.subtract.outer(G_Gv[G_myG, v], G_Gv[:, v])

        # Find length of dG and the normalized dG
        dG_myGG = np.linalg.norm(dG_myGGv, axis=2)
        dGn_myGGv = np.zeros_like(dG_myGGv)
        mask0 = np.where(dG_myGG != 0.)
        dGn_myGGv[mask0] = dG_myGGv[mask0] / dG_myGG[mask0][:, np.newaxis]

        return dG_myGGv, dG_myGG, dGn_myGGv

    def _get_densities_in_augmentation_sphere(self, calc, a):
        """Get the all-electron and smooth spin densities inside the
        augmentation spheres.

        Returns
        -------
        n_sLg : nd.array
            all-electron density
        nt_sLg : nd.array
            smooth density
        (s=spin, L=(l,m) spherical harmonic index, g=radial grid index)
        """
        setup = calc.wfs.setups[a]
        n_qg = setup.xc_correction.n_qg
        nt_qg = setup.xc_correction.nt_qg
        nc_g = setup.xc_correction.nc_g
        nct_g = setup.xc_correction.nct_g

        D_sp = calc.density.D_asp[a]
        B_pqL = setup.xc_correction.B_pqL
        D_sLq = np.inner(D_sp, B_pqL.T)
        nspins = len(D_sp)

        n_sLg = np.dot(D_sLq, n_qg)
        nt_sLg = np.dot(D_sLq, nt_qg)

        # Add core density
        n_sLg[:, 0] += np.sqrt(4. * np.pi) / nspins * nc_g
        nt_sLg[:, 0] += np.sqrt(4. * np.pi) / nspins * nct_g

        return n_sLg, nt_sLg

    @timer('Calculate PAW correction inside augmentation spheres')
    def _calculate_dfxc(self, calc, a):
        """Calculate the difference between fxc of the all-electron spin
        density and fxc of the smooth spin density.

        Returns
        -------
        df_ng : nd.array
            (f_ng - ft_ng) where (n=Lebedev index, g=radial grid index)
        Y_nL : nd.array
            real spherical harmonics on Lebedev quadrature
        rgd : GridDescriptor
            non-linear radial grid descriptor
        """
        # Extract spin densities from ground state calculation
        n_sLg, nt_sLg = self._get_densities_in_augmentation_sphere(calc, a)

        setup = calc.wfs.setups[a]
        Y_nL = setup.xc_correction.Y_nL
        rgd = setup.xc_correction.rgd
        f_g = rgd.zeros()
        ft_g = rgd.zeros()
        df_ng = np.array([rgd.zeros() for n in range(len(R_nv))])
        for n, Y_L in enumerate(Y_nL):
            f_g[:] = 0.
            n_sg = np.dot(Y_L, n_sLg)
            self.add_fxc(rgd, n_sg, f_g)

            ft_g[:] = 0.
            nt_sg = np.dot(Y_L, nt_sLg)
            self.add_fxc(rgd, nt_sg, ft_g)

            df_ng[n, :] = f_g - ft_g

        return df_ng, Y_nL, rgd

    def _ang_int(self, f_nA):
        """ Make surface integral on a sphere using Lebedev quadrature """
        return 4. * np.pi * np.tensordot(weight_n, f_nA, axes=([0], [0]))

    def _reduce_radial_grid(self, df_ng, rgd):
        """Reduce the radial grid, by excluding points where dfSns_g = 0,
        in order to avoid excess computation. Only points after the outermost
        point where dfSns_g is non-zero will be excluded.

        Returns
        -------
        df_ng : nd.array
            df_ng on reduced radial grid
        r_g : nd.array
            radius of each point on the reduced radial grid
        dv_g : nd.array
            volume element of each point on the reduced radial grid
        """
        # Find PAW correction range
        self.dfmask_g = np.where(self.dfSns_g > 0.)
        ng = np.max(self.dfmask_g) + 1

        # Integrate only r-values inside augmentation sphere
        df_ng = df_ng[:, :ng]

        r_g = rgd.r_g[:ng]
        dv_g = rgd.dv_g[:ng]

        return df_ng, r_g, dv_g

    @timer('Expand PAW correction in real spherical harmonics')
    def _perform_rshe(self, a, df_ng, Y_nL):
        """Perform expansion of dfxc in real spherical harmonics. Note that the
        Lebedev quadrature, which is used to calculate the expansion
        coefficients, is exact to order l=11. This implies that functions
        containing angular components l<=5 can be expanded exactly.

        Returns
        -------
        df_gL : nd.array
            dfxc in g=radial grid index, L=(l,m) spherical harmonic index
        """
        L_L = []
        l_L = []
        nL = min(Y_nL.shape[1], 36)
        # The convergence of the expansion is tracked through the
        # surface norm square of df
        self.dfSns_gL = np.repeat(self.dfSns_g,
                                  nL).reshape(self.dfSns_g.shape[0], nL)
        # Initialize convergence criteria
        self.rsheconvmin = 0.

        # Add real spherical harmonics to fulfill convergence criteria.
        df_gL = np.zeros((df_ng.shape[1], nL))
        l = 0
        while self.rsheconvmin < self.rshecc:
            if l > int(np.sqrt(nL) - 1):
                raise Exception('Could not expand %.f of' % self.rshecc
                                + ' PAW correction to atom %d in ' % a
                                + 'real spherical harmonics up to '
                                + 'order l=%d' % int(np.sqrt(nL) - 1))

            L_L += range(l**2, l**2 + 2 * l + 1)
            l_L += [l] * (2 * l + 1)

            self._add_rshe_coefficients(df_ng, df_gL, Y_nL, l)
            self._evaluate_rshe_convergence(df_gL)

            print('    At least a fraction of '
                  + '%.8f' % self.rsheconvmin
                  + ' of the PAW correction to atom %d could be ' % a
                  + 'expanded in spherical harmonics up to l=%d' % l,
                  file=self.fd)
            l += 1

        return df_gL

    def _add_rshe_coefficients(self, df_ng, df_gL, Y_nL, l):
        """
        Adds the l-components in the real spherical harmonic expansion
        of df_ng to df_gL.
        Assumes df_ng to be a real function.
        """
        nm = 2 * l + 1
        L_L = np.arange(l**2, l**2 + nm)
        df_ngm = np.repeat(df_ng, nm, axis=1).reshape((*df_ng.shape, nm))
        Y_ngm = np.repeat(Y_nL[:, L_L],
                          df_ng.shape[1], axis=0).reshape((*df_ng.shape, nm))
        df_gL[:, L_L] = self._ang_int(Y_ngm * df_ngm)

    def _evaluate_rshe_coefficients(self, f_gL):
        """
        Checks weither some (l,m)-coefficients are very small for all g,
        in that case they can be excluded from the expansion.
        """
        fc_L = np.sum(f_gL[self.dfmask_g]**2 / self.dfSns_gL[self.dfmask_g],
                      axis=0)
        self.rsheL_M = fc_L > (1. - self.rshecc) * 1.e-3

    def _evaluate_rshe_convergence(self, f_gL):
        """ The convergence of the real spherical harmonics expansion is
        tracked by comparing the surface norm square calculated using the
        expansion and the full result.
        
        Find also the minimal fraction of f_ng captured by the expansion
        in real spherical harmonics f_gL.
        """
        self._evaluate_rshe_coefficients(f_gL)

        rsheconv_g = np.ones(f_gL.shape[0])
        dfSns_g = np.sum(f_gL[:, self.rsheL_M]**2, axis=1)
        rsheconv_g[self.dfmask_g] = dfSns_g[self.dfmask_g]
        rsheconv_g[self.dfmask_g] /= self.dfSns_g[self.dfmask_g]

        self.rsheconvmin = np.min(rsheconv_g)

    def _reduce_rsh_expansion(self, df_gL):
        """Reduce the composite index L=(l,m) to M, which indexes non-zero
        coefficients in the expansion only.

        Returns
        -------
        df_gM : nd.array
            PAW correction in reduced rsh index
        L_M : nd.array
            L=(l,m) spherical harmonics indices in reduced rsh index
        l_M : list
            l spherical harmonics indices in reduced rsh index
        """
        # Recreate l_L array
        nL = df_gL.shape[1]
        l_L = []
        for l in range(int(np.sqrt(nL))):
            l_L += [l] * (2 * l + 1)
            
        # Filter away unused (l,m)-coefficients
        L_M = np.where(self.rsheL_M)[0]
        l_M = [l_L[L] for L in L_M]
        df_gM = df_gL[:, L_M]

        return df_gM, L_M, l_M

    @timer('Expand plane waves')
    def _expand_plane_waves(self, dG_myGG, dGn_myGGv, r_g, L_M, l_M):
        """Expand plane waves in spherical Bessel functions and real spherical
        harmonics:
                         l
                     __  __
         -iK.r       \   \      l             ^     ^
        e      = 4pi /   /  (-i)  j (|K|r) Y (K) Y (r)
                     ‾‾  ‾‾        l        lm    lm
                     l  m=-l

        Returns
        -------
        ii_MmyGG : nd.array
            (-i)^l for used (l,m) coefficients M
        j_gMmyGG : nd.array
            j_l(|dG|r) for used (l,m) coefficients M
        Y_MmyGG : nd.array
                 ^
            Y_lm(K) for used (l,m) coefficients M
        """
        # Setup arrays to fully vectorize computations
        nM = len(L_M)
        (r_gMmyGG, l_gMmyGG,
         dG_gMmyGG) = [a.reshape(len(r_g), nM, *dG_myGG.shape)
                       for a in np.meshgrid(r_g, l_M, dG_myGG.flatten(),
                                            indexing='ij')]

        with self.timer('Compute spherical bessel functions'):
            # Slow step. If it ever gets too slow, one can use the same
            # philosophy as _ft_from_grid, where dG=(G-G') results are
            # "unfolded" from a fourier transform to all unique K=dG
            # reciprocal lattice vectors. It should be possible to vectorize
            # the unfolding procedure to make it fast.
            j_gMmyGG = spherical_jn(l_gMmyGG, dG_gMmyGG * r_gMmyGG)

        Y_MmyGG = Yarr(L_M, dGn_myGGv)
        ii_MK = (-1j) ** np.repeat(l_M,
                                   np.prod(dG_myGG.shape))
        ii_MmyGG = ii_MK.reshape((nM, *dG_myGG.shape))

        return ii_MmyGG, j_gMmyGG, Y_MmyGG

    def add_fxc(self, gd, n_sg, fxc_g):
        raise NotImplementedError


class AdiabaticDensityKernelCalculator(AdiabaticKernelCalculator):

    def __init__(self, world=mpi.world, txt='-', timer=None,
                 rshe=0.99, ecut=None,
                 density_cut=None, **unused):
        """
        density_cut : float
            cutoff density below which f_xc is set to zero
        """
        AdiabaticKernelCalculator.__init__(self, world=world, txt=txt,
                                           timer=timer,
                                           rshe=rshe, ecut=ecut)

        self.density_cut = density_cut

        self.permitted_functionals += ['ALDA_x', 'ALDA_X', 'ALDA']

    def __call__(self, calc, pd, functional):

        Kxc_GG = AdiabaticKernelCalculator.__call__(self, calc, pd,
                                                    functional)

        if pd.kd.gamma:
            Kxc_GG[0, :] = 0.0
            Kxc_GG[:, 0] = 0.0

        return Kxc_GG

    def add_fxc(self, gd, n_sG, fxc_G):
        """
        Calculate fxc, using the cutoffs from input above

        ALDA_x is an explicit algebraic version
        ALDA_X uses the libxc package
        """

        _calculate_fxc = self._calculate_fxc
        density_cut = self.density_cut

        # Mask small n
        n_G = np.sum(n_sG, axis=0)
        if density_cut:
            npos_G = np.abs(n_G) > density_cut
        else:
            npos_G = np.full(np.shape(n_G), True, np.array(True).dtype)

        # Calculate fxc
        fxc_G[npos_G] += _calculate_fxc(gd, n_sG)[npos_G]

    def _calculate_fxc(self, gd, n_sG):
        if self.functional == 'ALDA_x':
            n_G = np.sum(n_sG, axis=0)
            fx_G = -1. / 3. * (3. / np.pi)**(1. / 3.) * n_G**(-2. / 3.)
            return fx_G
        else:
            fxc_sG = np.zeros_like(n_sG)
            xc = XC(self.functional[1:])
            xc.calculate_fxc(gd, n_sG, fxc_sG)

            return fxc_sG[0]


class AdiabaticTransverseKernelCalculator(AdiabaticKernelCalculator):

    def __init__(self, world=mpi.world, txt='-', timer=None,
                 rshe=0.99, ecut=None,
                 density_cut=None, spinpol_cut=None, **unused):
        """
        density_cut : float
            cutoff density below which f_xc is set to zero
        spinpol_cut : float
            cutoff spin polarization. Below, f_xc is evaluated in zeta=0 limit
        """
        AdiabaticKernelCalculator.__init__(self, world=world, txt=txt,
                                           timer=timer,
                                           rshe=rshe, ecut=ecut)

        self.density_cut = density_cut
        self.spinpol_cut = spinpol_cut

        self.permitted_functionals += ['ALDA_x', 'ALDA_X', 'ALDA']

    def add_fxc(self, gd, n_sG, fxc_G):
        """
        Calculate fxc, using the cutoffs from input above

        ALDA_x is an explicit algebraic version
        ALDA_X uses the libxc package
        """

        _calculate_pol_fxc = self._calculate_pol_fxc
        _calculate_unpol_fxc = self._calculate_unpol_fxc
        spinpol_cut = self.spinpol_cut
        density_cut = self.density_cut

        # Mask small zeta
        n_G, m_G = None, None
        if spinpol_cut is not None:
            m_G = n_sG[0] - n_sG[1]
            n_G = n_sG[0] + n_sG[1]
            zetasmall_G = np.abs(m_G / n_G) < spinpol_cut
        else:
            zetasmall_G = np.full(np.shape(n_sG[0]), False,
                                  np.array(False).dtype)

        # Mask small n
        if density_cut:
            if n_G is None:
                n_G = n_sG[0] + n_sG[1]
            npos_G = np.abs(n_G) > density_cut
        else:
            npos_G = np.full(np.shape(n_sG[0]), True, np.array(True).dtype)

        # Don't use small zeta limit if n is small
        zetasmall_G = np.logical_and(zetasmall_G, npos_G)

        # In small zeta limit, use unpolarized fxc
        if zetasmall_G.any():
            if n_G is None:
                n_G = n_sG[0] + n_sG[1]
            fxc_G[zetasmall_G] += _calculate_unpol_fxc(gd, n_G)[zetasmall_G]

        # Set fxc to zero if n is small
        allfine_G = np.logical_and(np.invert(zetasmall_G), npos_G)

        # Above both spinpol_cut and density_cut calculate polarized fxc
        if m_G is None:
            m_G = n_sG[0] - n_sG[1]
        fxc_G[allfine_G] += _calculate_pol_fxc(gd, n_sG, m_G)[allfine_G]

    def _calculate_pol_fxc(self, gd, n_sG, m_G):
        """ Calculate polarized fxc """

        assert np.shape(m_G) == np.shape(n_sG[0])

        if self.functional == 'ALDA_x':
            fx_G = - (6. / np.pi)**(1. / 3.) \
                * (n_sG[0]**(1. / 3.) - n_sG[1]**(1. / 3.)) / m_G
            return fx_G
        else:
            v_sG = np.zeros(np.shape(n_sG))
            xc = XC(self.functional[1:])
            xc.calculate(gd, n_sG, v_sg=v_sG)

            return (v_sG[0] - v_sG[1]) / m_G

    def _calculate_unpol_fxc(self, gd, n_G):
        """ Calculate unpolarized fxc """
        fx_G = - (3. / np.pi)**(1. / 3.) * 2. / 3. * n_G**(-2. / 3.)
        if self.functional in ('ALDA_x', 'ALDA_X'):
            return fx_G
        else:
            # From Perdew & Wang 1992
            A = 0.016887
            a1 = 0.11125
            b1 = 10.357
            b2 = 3.6231
            b3 = 0.88026
            b4 = 0.49671

            rs_G = 3. / (4. * np.pi) * n_G**(-1. / 3.)
            X_G = 2. * A * (b1 * rs_G**(1. / 2.)
                            + b2 * rs_G + b3 * rs_G**(3. / 2.) + b4 * rs_G**2.)
            ac_G = 2. * A * (1 + a1 * rs_G) * np.log(1. + 1. / X_G)

            fc_G = 2. * ac_G / n_G

            return fx_G + fc_G


class ATKCHack(AdiabaticTransverseKernelCalculator):
    """Hack to enable Goldstone scaling. Should go into
    AdiabaticTransverseKernelCalculator, when old format is removed."""
    def __init__(self, *args, fxc_scaling=None, **kwargs):
        AdiabaticTransverseKernelCalculator.__init__(self, *args, **kwargs)
        self.fxc_scaling = fxc_scaling

    def __call__(self, calc, pd, functional, kslrf=None, chiks_wGG=None):
        Kxc_GG = AdiabaticTransverseKernelCalculator.__call__(self, calc, pd,
                                                              functional)
        fxc_scaling = self.fxc_scaling

        if fxc_scaling is not None:
            assert isinstance(fxc_scaling[0], bool)
            if fxc_scaling[0]:
                if fxc_scaling[1] is None:
                    # Correct with new sign convention XXX
                    fxc_scaling[1] = find_Goldstone_scaling(pd, kslrf,
                                                            -chiks_wGG, Kxc_GG)

                assert isinstance(fxc_scaling[1], float)
                Kxc_GG *= fxc_scaling[1]

        return Kxc_GG
'''


def calculate_renormalized_kernel(pd, calc, functional, fd):
    """Renormalized kernel"""

    from gpaw.xc.fxc import KernelDens
    kernel = KernelDens(calc,
                        functional,
                        [pd.kd.bzk_kc[0]],
                        fd,
                        calc.wfs.kd.N_c,
                        None,
                        ecut=pd.ecut * Ha,
                        tag='',
                        timer=Timer())

    kernel.calculate_fhxc()
    r = Reader('fhxc_%s_%s_%s_%s.gpw' %
               ('', functional, pd.ecut * Ha, 0))
    Kxc_sGG = np.array([r.get('fhxc_sGsG')])

    v_G = 4 * np.pi / pd.G2_qG[0]
    Kxc_sGG[0] -= np.diagflat(v_G)

    if pd.kd.gamma:
        Kxc_sGG[:, 0, :] = 0.0
        Kxc_sGG[:, :, 0] = 0.0

    return Kxc_sGG


def calculate_lr_kernel(pd, calc, alpha=0.2):
    """Long range kernel: fxc = \alpha / |q+G|^2"""

    assert pd.kd.gamma

    f_G = np.zeros(len(pd.G2_qG[0]))
    f_G[0] = -alpha
    f_G[1:] = -alpha / pd.G2_qG[0][1:]

    return np.array([np.diag(f_G)])


def calculate_dm_kernel(pd, calc):
    """Density matrix kernel"""

    assert pd.kd.gamma

    nv = calc.wfs.setups.nvalence
    psit_nG = np.array([calc.wfs.kpt_u[0].psit_nG[n]
                        for n in range(4 * nv)])
    vol = np.linalg.det(calc.wfs.gd.cell_cv)
    Ng = np.prod(calc.wfs.gd.N_c)
    rho_GG = np.dot(psit_nG.conj().T, psit_nG) * vol / Ng**2

    maxG2 = np.max(pd.G2_qG[0])
    cut_G = np.arange(calc.wfs.pd.ngmax)[calc.wfs.pd.G2_qG[0] <= maxG2]

    G_G = pd.G2_qG[0]**0.5
    G_G[0] = 1.0

    Kxc_GG = np.diagflat(4 * np.pi / G_G**2)
    Kxc_GG = np.dot(Kxc_GG, rho_GG.take(cut_G, 0).take(cut_G, 1))
    Kxc_GG -= 4 * np.pi * np.diagflat(1.0 / G_G**2)

    return np.array([Kxc_GG])


def get_bootstrap_kernel(pd, chi0, chi0_wGG, fd):
    """ Bootstrap kernel (see below) """

    if chi0.world.rank == 0:
        chi0_GG = chi0_wGG[0]
        if chi0.world.size > 1:
            # If size == 1, chi0_GG is not contiguous, and broadcast()
            # will fail in debug mode.  So we skip it until someone
            # takes a closer look.
            chi0.world.broadcast(chi0_GG, 0)
    else:
        nG = pd.ngmax
        chi0_GG = np.zeros((nG, nG), complex)
        chi0.world.broadcast(chi0_GG, 0)

    return calculate_bootstrap_kernel(pd, chi0_GG, fd)


def calculate_bootstrap_kernel(pd, chi0_GG, fd):
    """Bootstrap kernel PRL 107, 186401"""

    if pd.kd.gamma:
        v_G = np.zeros(len(pd.G2_qG[0]))
        v_G[0] = 4 * np.pi
        v_G[1:] = 4 * np.pi / pd.G2_qG[0][1:]
    else:
        v_G = 4 * np.pi / pd.G2_qG[0]

    nG = len(v_G)
    K_GG = np.diag(v_G)

    fxc_GG = np.zeros((nG, nG), dtype=complex)
    dminv_GG = np.zeros((nG, nG), dtype=complex)

    for iscf in range(120):
        dminvold_GG = dminv_GG.copy()
        Kxc_GG = K_GG + fxc_GG

        chi_GG = np.dot(np.linalg.inv(np.eye(nG, nG)
                                      - np.dot(chi0_GG, Kxc_GG)), chi0_GG)
        dminv_GG = np.eye(nG, nG) + np.dot(K_GG, chi_GG)

        alpha = dminv_GG[0, 0] / (K_GG[0, 0] * chi0_GG[0, 0])
        fxc_GG = alpha * K_GG
        print(iscf, 'alpha =', alpha, file=fd)
        error = np.abs(dminvold_GG - dminv_GG).sum()
        if np.sum(error) < 0.1:
            print('Self consistent fxc finished in %d iterations !' % iscf,
                  file=fd)
            break
        if iscf > 100:
            print('Too many fxc scf steps !', file=fd)

    return np.array([fxc_GG])
