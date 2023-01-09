"""Contains methods for calculating local LR-TDDFT kernels."""

from pathlib import Path
from functools import partial

import numpy as np
from scipy.special import spherical_jn

from ase.units import Bohr

from gpaw.xc import XC
from gpaw.spherical_harmonics import Yarr
from gpaw.sphere.lebedev import weight_n, R_nv

from gpaw.response import ResponseGroundStateAdapter, ResponseContext, timer
from gpaw.response.chiks import ChiKS
from gpaw.response.goldstone import get_goldstone_scaling
from gpaw.response.localft import (LocalFTCalculator,
                                   add_LDA_dens_fxc, add_LSDA_trans_fxc)


class FXCScaling:
    """Helper for scaling fxc kernels."""

    def __init__(self, mode, lambd=None):
        self.mode = mode
        self.lambd = lambd

    def get_scaling(self, *args):
        if self.lambd is None:
            self.lambd = self.calculate_scaling(*args)

        return self.lambd

    def calculate_scaling(self, chiks, Kxc_GG):
        if chiks.spincomponent in ['+-', '-+']:
            return get_goldstone_scaling(self.mode, chiks, Kxc_GG)
        else:
            raise ValueError('No scaling method implemented for '
                             f'spincomponent={chiks.spincomponent}')


class FXCFactory:
    """Exchange-correlation kernel factory."""

    def __init__(self,
                 gs: ResponseGroundStateAdapter,
                 context: ResponseContext):
        self.gs = gs
        self.context = context

    def __call__(self, fxc, chiks: ChiKS,
                 calculator={'method': 'old',
                             'rshelmax': -1,
                             'rshewmin': None},
                 filename=None,
                 fxc_scaling=None):
        """Get the fxc kernel Kxc_GG.

        Parameters
        ----------
        fxc : str
            Approximation to the (local) xc kernel.
            Choices: ALDA, ALDA_X, ALDA_x
        calculator : dict
            Parameters to set up the FXCCalculator. The 'method' key
            determines what calculator is initilized and remaining parameters
            are passed to the calculator as key-word arguments.
        filename : str
            Store a calculated kernel as a .npy file buffer with the given file
            name. For subsequent calls to the kernel factory with the same file
            name, the factory will read the kernel from the file, instead of
            recalculating it.
            NB: If you want to calculate a NEW kernel (according to some change
            in input parameters), you will have to use a NEW filename.
        fxc_scaling : None or FXCScaling
        """
        if self.file_buffer_exists(filename):
            Kxc_GG = self.read(filename)
        else:
            Kxc_GG = self.calculate(fxc, chiks, calculator=calculator)
            self.write(Kxc_GG, filename)

        if fxc_scaling is not None:
            lambd = fxc_scaling.get_scaling(chiks, Kxc_GG)
            self.context.print(r'Rescaling the xc-kernel by a factor of λ='
                               f'{lambd}')
            Kxc_GG *= lambd

        return Kxc_GG

    @staticmethod
    def file_buffer_exists(filename):
        if filename is None:
            return False
        return Path(filename).is_file()

    @staticmethod
    def write(Kxc_GG, filename):
        if filename is not None:
            np.save(filename, Kxc_GG)

    @staticmethod
    def read(filename):
        return np.load(filename)

    def calculate(self, fxc, chiks, *, calculator):
        """In-place calculation of the bare (unscaled) fxc kernel Kxc_GG."""
        assert isinstance(calculator, dict) and 'method' in calculator

        # Generate the desired calculator
        calc_kwargs = calculator.copy()
        method = calc_kwargs.pop('method')
        fxc_calculator = self.get_fxc_calculator(method=method, **calc_kwargs)

        Kxc_GG = fxc_calculator(fxc, chiks.spincomponent, chiks.pd)

        return Kxc_GG

    def get_fxc_calculator(self, *, method, **calc_kwargs):
        """Factory function for initializing fxc calculators."""
        fxc_calculator_cls = self.get_fxc_calculator_cls(method)

        return fxc_calculator_cls(self.gs, self.context, **calc_kwargs)

    @staticmethod
    def get_fxc_calculator_cls(method):
        """Factory function for selecting fxc calculators."""
        if method == 'old':
            return OldAdiabaticFXCCalculator
        elif method == 'new':
            return NewAdiabaticFXCCalculator

        raise ValueError(f'Invalid fxc calculator method {method}')


class NewAdiabaticFXCCalculator:
    """Calculator for adiabatic local exchange-correlation kernels."""

    def __init__(self, gs, context, localft_calc: LocalFTCalculator):
        """Contruct the fxc calculator based on a local FT calculator."""
        self.localft_calc = localft_calc

        self.context = localft_calc.context

    @timer('Calculate XC kernel')
    def __call__(self, fxc, spincomponent, pd):
        """Calculate the fxc kernel."""
        # Calculate fxc(G)
        add_fxc = self.create_add_fxc(fxc, spincomponent)
        fxc_G = self.localft_calc(pd, add_fxc)

        # Unfold the kernel to Kxc_GG' = fxc(G-G')
        Kxc_GG = self.unfold_kernel(pd, fxc_G)

        return Kxc_GG

    @staticmethod
    def create_add_fxc(fxc, spincomponent):
        """Creator component to set up the right calculation."""
        assert fxc in ['ALDA_x', 'ALDA_X', 'ALDA']

        if spincomponent in ['00', 'uu', 'dd']:
            add_fxc = partial(add_LDA_dens_fxc, fxc=fxc)
        elif spincomponent in ['+-', '-+']:
            add_fxc = partial(add_LSDA_trans_fxc, fxc=fxc)
        else:
            raise ValueError(spincomponent)

        return add_fxc

    @timer('Unfold kernel')
    def unfold_kernel(self, pd, fxc_G):
        """
        Some documentation here!                                               XXX
        """
        # Calculate (G-G') reciprocal space vectors
        dG_GGv = self.calculate_dG(pd)

        # Reshape to composite K = (G, G') index
        dG_Kv = dG_GGv.reshape(-1, dG_GGv.shape[-1])

        # Find unique dG-vectors
        dG_dGv, dG_K = np.unique(dG_Kv, return_inverse=True, axis=0)
        ndG = len(dG_dGv)

        # Create mapping from (G-G') index dG to fxc(G) index G
        dG_G = self.map_to_dG(pd, dG_dGv)

        # Unfold fxc(G) to fxc(G-G')
        fxc_dG = np.zeros(ndG, dtype=complex)
        fxc_dG[dG_G] = fxc_G

        # Unfold fxc(G-G') to Kxc_GG'
        Kxc_GG = fxc_dG[dG_K].reshape(dG_GGv.shape[:2])

        return Kxc_GG

    @staticmethod
    def calculate_dG(pd):
        """
        Some documentation here!                                               XXX
        """
        nG = pd.ngmax
        G_Gv = pd.get_reciprocal_vectors(add_q=False)

        dG_GGv = np.zeros((nG, nG, 3))
        for v in range(3):
            dG_GGv[:, :, v] = np.subtract.outer(G_Gv[:, v], G_Gv[:, v])

        return dG_GGv
        
    @staticmethod
    def map_to_dG(pd, dG_dGv):
        """
        Some documentation here!                                               XXX
        """
        G_Gv = pd.get_reciprocal_vectors(add_q=False)

        dG_G = []
        for G_v in G_Gv:
            diff_dG = np.linalg.norm(dG_dGv - G_v, axis=1)
            dG = np.argmin(diff_dG)
            assert diff_dG[dG] < 1e-6, f'{diff_dG[dG]}, {G_v}'
            dG_G.append(dG)

        return np.array(dG_G)


class OldAdiabaticFXCCalculator:
    """Calculator for adiabatic local exchange-correlation kernels in pw mode.
    """

    def __init__(self, gs, context, rshelmax=-1, rshewmin=None):
        """Construct the AdiabaticFXCCalculator.

        Parameters
        ----------
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
        self.gs = gs
        self.context = context

        # Do not carry out the expansion in real spherical harmonics, if lmax
        # is chosen as None
        self.rshe = rshelmax is not None

        if self.rshe:
            # Perform rshe up to l<=lmax(<=5)
            if rshelmax == -1:
                self.rshelmax = 5
            else:
                assert isinstance(rshelmax, int)
                assert rshelmax in range(6)
                self.rshelmax = rshelmax

            self.rshewmin = rshewmin if rshewmin is not None else 0.
            self.dfmask_g = None

    @timer('Calculate XC kernel')
    def __call__(self, fxc, spincomponent, pd):
        """Calculate the fxc kernel."""
        self.set_up_calculation(fxc, spincomponent)

        self.context.print('Calculating fxc')
        # Get the spin density we need and allocate fxc
        n_sG = self.get_density_on_grid(pd.gd)
        fxc_G = np.zeros(np.shape(n_sG[0]))

        self.context.print('    Calculating fxc on real space grid')
        self._add_fxc(pd.gd, n_sG, fxc_G)

        # Fourier transform to reciprocal space
        Kxc_GG = self.ft_from_grid(fxc_G, pd)

        if self.rshe:  # Do PAW correction to Fourier transformed kernel
            KxcPAW_GG = self.calculate_kernel_paw_correction(pd)
            Kxc_GG += KxcPAW_GG

        self.context.print('Finished calculating fxc\n')

        return Kxc_GG / pd.gd.volume

    def get_density_on_grid(self, gd):
        """Get the spin density on coarse real-space grid.

        Returns
        -------
        nt_sG or n_sG : nd.array
            Spin density on coarse real-space grid. If not self.rshe, use
            the PAW corrected all-electron spin density.
        """
        if self.rshe:
            return self.gs.nt_sR  # smooth density

        self.context.print('    Calculating all-electron density')
        with self.context.timer('Calculating all-electron density'):
            n_sG, gd1 = self.gs.all_electron_density(gridrefinement=1)
            assert (gd1.n_c == gd.n_c).all()
            assert gd1.comm.size == 1
            return n_sG

    @timer('Fourier transform of kernel from real-space grid')
    def ft_from_grid(self, fxc_G, pd):
        self.context.print('    Fourier transforming kernel from real-space')
        nG = pd.gd.N_c
        nG0 = nG[0] * nG[1] * nG[2]

        tmp_g = np.fft.fftn(fxc_G) * pd.gd.volume / nG0

        # The unfolding procedure could use vectorization and parallelization.
        # This remains a slow step for now.
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
    def calculate_kernel_paw_correction(self, pd):
        self.context.print("    Calculating PAW corrections to the kernel\n")

        # Calculate (G-G') reciprocal space vectors
        dG_GGv = self._calculate_dG(pd)

        # Reshape to composite K = (G, G') index
        dG_Kv = dG_GGv.reshape(-1, dG_GGv.shape[-1])

        # Find unique dG-vectors
        dG_dGv, dG_K = np.unique(dG_Kv, return_inverse=True, axis=0)
        ndG = len(dG_dGv)

        # Allocate array and distribute plane waves
        KxcPAW_dG = np.zeros(ndG, dtype=complex)
        dG_mydG = self._distribute_correction(ndG)
        dG_mydGv = dG_dGv[dG_mydG]

        # Calculate my (G-G') reciprocal space vector lengths and directions
        dGl_mydG, dGn_mydGv = self._normalize_by_length(dG_mydGv)

        # Calculate PAW correction to each augmentation sphere (to each atom)
        R_av = self.gs.atoms.positions / Bohr
        for a, R_v in enumerate(R_av):
            # Calculate dfxc on Lebedev quadrature and radial grid
            # Please note: Using the radial grid descriptor with _add_fxc
            # might give problems beyond ALDA
            df_ng, Y_nL, rgd = self._calculate_dfxc(a)

            # Calculate the surface norm square of df
            dfSns_g = self._ang_int(df_ng ** 2)
            # Reduce radial grid by excluding points where dfSns_g = 0
            df_ng, r_g, dv_g = self._reduce_radial_grid(df_ng, rgd, dfSns_g)

            # Expand correction in real spherical harmonics
            df_gL = self._perform_rshe(df_ng, Y_nL)
            # Reduce expansion by removing coefficients that do not contribute
            df_gM, L_M, l_M = self._reduce_rsh_expansion(a, df_gL, dfSns_g)

            # Expand plane wave differences (G-G')
            (ii_MmydG,
             j_gMmydG,
             Y_MmydG) = self._expand_plane_waves(dGl_mydG, dGn_mydGv,
                                                 r_g, L_M, l_M)

            # Perform integration
            with self.context.timer('Integrate PAW correction'):
                coefatomR_dG = np.exp(-1j * np.inner(dG_mydGv, R_v))
                coefatomang_MdG = ii_MmydG * Y_MmydG
                coefatomrad_MdG = np.tensordot(j_gMmydG * df_gL[:, L_M,
                                                                np.newaxis],
                                               dv_g, axes=([0, 0]))
                coefatom_dG = np.sum(coefatomang_MdG * coefatomrad_MdG, axis=0)
                KxcPAW_dG[dG_mydG] += coefatom_dG * coefatomR_dG

        self.context.world.sum(KxcPAW_dG)

        # Unfold PAW correction
        KxcPAW_GG = KxcPAW_dG[dG_K].reshape(dG_GGv.shape[:2])

        return KxcPAW_GG

    def _calculate_dG(self, pd):
        """Calculate (G-G') reciprocal space vectors"""
        world = self.context.world
        npw = pd.ngmax
        G_Gv = pd.get_reciprocal_vectors(add_q=False)

        # Distribute dG to calculate
        nGpr = (npw + world.size - 1) // world.size
        Ga = min(world.rank * nGpr, npw)
        Gb = min(Ga + nGpr, npw)
        G_myG = range(Ga, Gb)

        # Calculate dG_v for every set of (G-G')
        dG_GGv = np.zeros((npw, npw, 3))
        for v in range(3):
            dG_GGv[Ga:Gb, :, v] = np.subtract.outer(G_Gv[G_myG, v], G_Gv[:, v])
        world.sum(dG_GGv)

        return dG_GGv

    def _distribute_correction(self, ndG):
        """Distribute correction"""
        world = self.context.world
        ndGpr = (ndG + world.size - 1) // world.size
        dGa = min(world.rank * ndGpr, ndG)
        dGb = min(dGa + ndGpr, ndG)

        return range(dGa, dGb)

    @staticmethod
    def _normalize_by_length(dG_mydGv):
        """Find the length and direction of reciprocal lattice vectors."""
        dGl_mydG = np.linalg.norm(dG_mydGv, axis=1)
        dGn_mydGv = np.zeros_like(dG_mydGv)
        mask0 = np.where(dGl_mydG != 0.)
        dGn_mydGv[mask0] = dG_mydGv[mask0] / dGl_mydG[mask0][:, np.newaxis]

        return dGl_mydG, dGn_mydGv

    def _get_densities_in_augmentation_sphere(self, a):
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
        setup = self.gs.setups[a]
        n_qg = setup.xc_correction.n_qg
        nt_qg = setup.xc_correction.nt_qg
        nc_g = setup.xc_correction.nc_g
        nct_g = setup.xc_correction.nct_g

        D_sp = self.gs.D_asp[a]
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
    def _calculate_dfxc(self, a):
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
        n_sLg, nt_sLg = self._get_densities_in_augmentation_sphere(a)

        setup = self.gs.setups[a]
        Y_nL = setup.xc_correction.Y_nL
        rgd = setup.xc_correction.rgd
        f_g = rgd.zeros()
        ft_g = rgd.zeros()
        df_ng = np.array([rgd.zeros() for n in range(len(R_nv))])
        for n, Y_L in enumerate(Y_nL):
            f_g[:] = 0.
            n_sg = np.dot(Y_L, n_sLg)
            self._add_fxc(rgd, n_sg, f_g)

            ft_g[:] = 0.
            nt_sg = np.dot(Y_L, nt_sLg)
            self._add_fxc(rgd, nt_sg, ft_g)

            df_ng[n, :] = f_g - ft_g

        return df_ng, Y_nL, rgd

    @staticmethod
    def _ang_int(f_nA):
        """ Make surface integral on a sphere using Lebedev quadrature """
        return 4. * np.pi * np.tensordot(weight_n, f_nA, axes=([0], [0]))

    def _reduce_radial_grid(self, df_ng, rgd, dfSns_g):
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
        self.dfmask_g = np.where(dfSns_g > 0.)
        ng = np.max(self.dfmask_g) + 1

        # Integrate only r-values inside augmentation sphere
        df_ng = df_ng[:, :ng]

        r_g = rgd.r_g[:ng]
        dv_g = rgd.dv_g[:ng]

        return df_ng, r_g, dv_g

    @timer('Expand PAW correction in real spherical harmonics')
    def _perform_rshe(self, df_ng, Y_nL):
        """Perform expansion of dfxc in real spherical harmonics. Note that the
        Lebedev quadrature, which is used to calculate the expansion
        coefficients, is exact to order l=11. This implies that functions
        containing angular components l<=5 can be expanded exactly.
        Assumes df_ng to be a real function.

        Returns
        -------
        df_gL : nd.array
            dfxc in g=radial grid index, L=(l,m) spherical harmonic index
        """
        lmax = min(int(np.sqrt(Y_nL.shape[1])) - 1, 36)
        nL = (lmax + 1)**2
        L_L = np.arange(nL)

        # Perform the real spherical harmonics expansion
        df_ngL = np.repeat(df_ng, nL, axis=1).reshape((*df_ng.shape, nL))
        Y_ngL = np.repeat(Y_nL[:, L_L], df_ng.shape[1],
                          axis=0).reshape((*df_ng.shape, nL))
        df_gL = self._ang_int(Y_ngL * df_ngL)

        return df_gL

    def _reduce_rsh_expansion(self, a, df_gL, dfSns_g):
        """Reduce the composite index L=(l,m) to M, which indexes coefficients
        contributing with a weight larger than rshewmin to the surface norm
        square on average.

        Returns
        -------
        df_gM : nd.array
            PAW correction in reduced rsh index
        L_M : nd.array
            L=(l,m) spherical harmonics indices in reduced rsh index
        l_M : list
            l spherical harmonics indices in reduced rsh index
        """
        # Create L_L and l_L array
        lmax = min(self.rshelmax, int(np.sqrt(df_gL.shape[1])) - 1)
        nL = (lmax + 1)**2
        L_L = np.arange(nL)
        l_L = []
        for l in range(int(np.sqrt(nL))):
            l_L += [l] * (2 * l + 1)

        # Filter away (l,m)-coefficients that do not contribute
        rshew_L = self._evaluate_rshe_coefficients(a, nL, df_gL, dfSns_g)
        L_M = np.where(rshew_L[L_L] > self.rshewmin)[0]
        l_M = [l_L[L] for L in L_M]
        df_gM = df_gL[:, L_M]

        return df_gM, L_M, l_M

    def _evaluate_rshe_coefficients(self, a, nL, df_gL, dfSns_g):
        """If some of the rshe coefficients are very small for all radii g,
        we may choose to exclude them from the kernel PAW correction.

        The "smallness" is evaluated from their average weight in
        evaluating the surface norm square for each radii g.
        """
        # Compute each coefficient's fraction of the surface norm square
        nallL = df_gL.shape[1]
        dfSns_gL = np.repeat(dfSns_g, nallL).reshape(dfSns_g.shape[0], nallL)
        dfSw_gL = df_gL[self.dfmask_g] ** 2 / dfSns_gL[self.dfmask_g]

        # The smallness is evaluated from the average
        rshew_L = np.average(dfSw_gL, axis=0)

        # Print information about the expansion
        p = partial(self.context.print, flush=False)
        p('    RSHE of atom', a)
        p('      {0:6}  {1:10}  {2:10}  {3:8}'.format('(l,m)',
                                                      'max weight',
                                                      'avg weight',
                                                      'included'))
        for L, (dfSw_g, rshew) in enumerate(zip(dfSw_gL.T, rshew_L)):
            self.print_rshe_info(L, nL, dfSw_g, rshew)

        tot_avg_cov = np.average(np.sum(dfSw_gL, axis=1))
        avg_cov = np.average(np.sum(dfSw_gL[:, :nL]
                                    [:, rshew_L[:nL] > self.rshewmin], axis=1))
        p(f'      In total: {avg_cov} of the dfSns is covered on average')
        self.context.print(f'      In total: {tot_avg_cov} of the dfSns could',
                           'be covered on average\n')

        return rshew_L

    def print_rshe_info(self, L, nL, dfSw_g, rshew):
        """Print information about the importance of the rshe coefficient"""
        l = int(np.sqrt(L))
        m = L - l**2 - l
        included = 'yes' if (rshew > self.rshewmin and L < nL) else 'no'
        p = partial(self.context.print, flush=False)
        p('      {0:6}  {1:1.8f}  {2:1.8f}  {3:8}'.format(f'({l},{m})',
                                                          np.max(dfSw_g),
                                                          rshew, included))

    @timer('Expand plane waves')
    def _expand_plane_waves(self, dG_mydG, dGn_mydGv, r_g, L_M, l_M):
        r"""Expand plane waves in spherical Bessel functions and real spherical
        harmonics:
                         l
                     __  __
         -iK.r       \   \      l             ^     ^
        e      = 4pi /   /  (-i)  j (|K|r) Y (K) Y (r)
                     ‾‾  ‾‾        l        lm    lm
                     l  m=-l

        Returns
        -------
        ii_MmydG : nd.array
            (-i)^l for used (l,m) coefficients M
        j_gMmydG : nd.array
            j_l(|dG|r) for used (l,m) coefficients M
        Y_MmydG : nd.array
                 ^
            Y_lm(K) for used (l,m) coefficients M
        """
        nmydG = len(dG_mydG)
        # Setup arrays to fully vectorize computations
        nM = len(L_M)
        (r_gMmydG, l_gMmydG,
         dG_gMmydG) = [a.reshape(len(r_g), nM, nmydG)
                       for a in np.meshgrid(r_g, l_M, dG_mydG, indexing='ij')]

        with self.context.timer('Compute spherical bessel functions'):
            # Slow step
            j_gMmydG = spherical_jn(l_gMmydG, dG_gMmydG * r_gMmydG)

        Y_MmydG = Yarr(L_M, dGn_mydGv)
        ii_X = (-1j) ** np.repeat(l_M, nmydG)
        ii_MmydG = ii_X.reshape((nM, nmydG))

        return ii_MmydG, j_gMmydG, Y_MmydG

    def set_up_calculation(self, fxc, spincomponent):
        """Creator component to set up the right calculation."""
        assert fxc in ['ALDA_x', 'ALDA_X', 'ALDA']

        if spincomponent in ['00', 'uu', 'dd']:
            assert len(self.gs.nt_sR) == 1  # nspins, see XXX below

            self._calculate_fxc = partial(self.calculate_dens_fxc, fxc=fxc)
        elif spincomponent in ['+-', '-+']:
            assert len(self.gs.nt_sR) == 2  # nspins

            self._calculate_fxc = partial(self.calculate_trans_fxc, fxc=fxc)
        else:
            raise ValueError(spincomponent)

    def _add_fxc(self, gd, n_sG, fxc_G):
        """Calculate fxc and add it to the output array."""
        fxc_G += self._calculate_fxc(gd, n_sG)

    def calculate_dens_fxc(self, gd, n_sG, *, fxc):
        if fxc == 'ALDA_x':
            n_G = np.sum(n_sG, axis=0)
            fx_G = -1. / 3. * (3. / np.pi)**(1. / 3.) * n_G**(-2. / 3.)
            return fx_G

        assert len(n_sG) == 1
        from gpaw.xc.libxc import LibXC
        kernel = LibXC(fxc[1:])
        fxc_sG = np.zeros_like(n_sG)
        kernel.xc.calculate_fxc_spinpaired(n_sG.ravel(), fxc_sG)

        return fxc_sG[0]  # not tested for spin-polarized calculations XXX

    def calculate_trans_fxc(self, gd, n_sG, *, fxc):
        """Calculate polarized fxc of spincomponents '+-', '-+'."""
        m_G = n_sG[0] - n_sG[1]

        if fxc == 'ALDA_x':
            fx_G = - (6. / np.pi)**(1. / 3.) \
                * (n_sG[0]**(1. / 3.) - n_sG[1]**(1. / 3.)) / m_G
            return fx_G
        else:
            v_sG = np.zeros(np.shape(n_sG))
            xc = XC(fxc[1:])
            xc.calculate(gd, n_sG, v_sg=v_sG)

            return (v_sG[0] - v_sG[1]) / m_G
