from __future__ import annotations

from abc import abstractmethod

import numpy as np

from gpaw.response import (ResponseGroundStateAdapter, ResponseContext,
                           GPWFilename, TXTFilename,
                           ensure_gs, ensure_gs_and_context)
from gpaw.response.symmetry import SymmetryAnalyzer
from gpaw.response.frequencies import ComplexFrequencyDescriptor
from gpaw.response.chiks import ChiKSCalculator, smat
from gpaw.response.localft import LocalFTCalculator, add_LSDA_Wxc
from gpaw.response.site_kernels import SiteKernels
from gpaw.response.site_data import AtomicSites, AtomicSiteData
from gpaw.response.pair_functions import SingleQPWDescriptor, PairFunction
from gpaw.response.pair_integrator import PairFunctionIntegrator
from gpaw.response.pair_transitions import PairTransitions
from gpaw.response.matrix_elements import (SitePairDensityCalculator,
                                           SiteZeemanPairEnergyCalculator)

from ase.units import Hartree


class IsotropicExchangeCalculator:
    r"""Calculator class for the Heisenberg exchange constants

    _           2
    J^ab(q) = - ‾‾ B^(xc†) K^(a†)(q) χ_KS^('+-)(q) K^b(q) B^(xc)            (1)
                V0

    calculated for an isotropic system in a plane wave representation using
    the magnetic force theorem within second order perturbation theory, see
    [J. Phys.: Condens. Matter 35 (2023) 105802].

    Entering the formula for the isotropic exchange constant at wave vector q
    between sublattice a and b is the unit cell volume V0, the functional
    derivative of the (LDA) exchange-correlation energy with respect to the
    magnitude of the magnetization B^(xc), the sublattice site kernels K^a(q)
    and K^b(q) as well as the reactive part of the static transverse magnetic
    susceptibility of the Kohn-Sham system χ_KS^('+-)(q).

    NB: To achieve numerical stability of the plane-wave implementation, we
    use instead the following expression to calculate exchange parameters:

    ˷           2
    J^ab(q) = - ‾‾ W_xc^(z†) K^(a†)(q) χ_KS^('+-)(q) K^b(q) W_xc^z          (2)
                V0

    We do this since B^(xc)(r) = -|W_xc^z(r)| is nonanalytic in points of space
    where the spin-polarization changes sign, why it is problematic to evaluate
    Eq. (1) numerically within a plane-wave representation.
    If the site partitionings only include spin-polarization of the same sign,
    Eqs. (1) and (2) should yield identical exchange parameters, but for
    antiferromagnetically aligned sites, the coupling constants differ by a
    sign.

    The site kernels encode the partitioning of real space into sites of the
    Heisenberg model. This is not a uniquely defined procedure, why the user
    has to define them externally through the SiteKernels interface."""

    def __init__(self,
                 chiks_calc: ChiKSCalculator,
                 localft_calc: LocalFTCalculator):
        """Construct the IsotropicExchangeCalculator object."""
        # Check that chiks has the assumed properties
        assumed_props = dict(
            gammacentered=True,
            nblocks=1
        )
        for key, item in assumed_props.items():
            assert getattr(chiks_calc, key) == item, \
                f'Expected chiks.{key} == {item}. '\
                f'Got: {getattr(chiks_calc, key)}'

        self.chiks_calc = chiks_calc
        self.context = chiks_calc.context

        # Check assumed properties of the LocalFTCalculator
        assert localft_calc.context is self.context
        assert localft_calc.gs is chiks_calc.gs
        self.localft_calc = localft_calc

        # W_xc^z buffer
        self._Wxc_G = None

        # χ_KS^('+-) buffer
        self._chiksr = None

    def __call__(self, q_c, site_kernels: SiteKernels, txt=None):
        """Calculate the isotropic exchange constants for a given wavevector.

        Parameters
        ----------
        q_c : nd.array
            Wave vector q in relative coordinates
        site_kernels : SiteKernels
            Site kernels instance defining the magnetic sites of the crystal
        txt : str
            Separate file to store the chiks calculation output in (optional).
            If not supplied, the output will be written to the standard text
            output location specified when initializing chiks.

        Returns
        -------
        J_abp : nd.array (dtype=complex)
            Isotropic Heisenberg exchange constants between magnetic sites a
            and b for all the site partitions p given by the site_kernels.
        """
        # Get ingredients
        Wxc_G = self.get_Wxc()
        chiksr = self.get_chiksr(q_c, txt=txt)
        qpd, chiksr_GG = chiksr.qpd, chiksr.array[0]  # array = chiksr_zGG
        V0 = qpd.gd.volume

        # Allocate an array for the exchange constants
        nsites = site_kernels.nsites
        J_pab = np.empty(site_kernels.shape + (nsites,), dtype=complex)

        # Compute exchange coupling
        for J_ab, K_aGG in zip(J_pab, site_kernels.calculate(qpd)):
            for a in range(nsites):
                for b in range(nsites):
                    J = np.conj(Wxc_G) @ np.conj(K_aGG[a]).T @ chiksr_GG \
                        @ K_aGG[b] @ Wxc_G
                    J_ab[a, b] = - 2. * J / V0

        # Transpose to have the partitions index last
        J_abp = np.transpose(J_pab, (1, 2, 0))

        return J_abp * Hartree  # Convert from Hartree to eV

    def get_Wxc(self):
        """Get B^(xc)_G from buffer."""
        if self._Wxc_G is None:  # Calculate if buffer is empty
            self._Wxc_G = self._calculate_Wxc()

        return self._Wxc_G

    def _calculate_Wxc(self):
        """Calculate the Fourier transform W_xc^z(G)."""
        # Create a plane wave descriptor encoding the plane wave basis. Input
        # q_c is arbitrary, since we are assuming that chiks.gammacentered == 1
        qpd0 = self.chiks_calc.get_pw_descriptor([0., 0., 0.])

        return self.localft_calc(qpd0, add_LSDA_Wxc)

    def get_chiksr(self, q_c, txt=None):
        """Get χ_KS^('+-)(q) from buffer."""
        q_c = np.asarray(q_c)

        # Calculate if buffer is empty or a new q-point is given
        if self._chiksr is None or not np.allclose(q_c, self._chiksr.q_c):
            self._chiksr = self._calculate_chiksr(q_c, txt=txt)

        return self._chiksr

    def _calculate_chiksr(self, q_c, txt=None):
        r"""Use the ChiKSCalculator to calculate the reactive part of the
        static Kohn-Sham susceptibility χ_KS^('+-)(q).

        First, the dynamic Kohn-Sham susceptibility

                                 __  __
                              1  \   \        f_nk↑ - f_mk+q↓
        χ_KS,GG'^+-(q,ω+iη) = ‾  /   /  ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
                              V  ‾‾  ‾‾ ħω - (ε_mk+q↓ - ε_nk↑) + iħη
                                 k  n,m
                                        x n_nk↑,mk+q↓(G+q) n_mk+q↓,nk↑(-G'-q)

        is calculated in the static limit ω=0 and without broadening η=0. Then,
        the reactive part (see [PRB 103, 245110 (2021)]) is extracted:

                              1
        χ_KS,GG'^(+-')(q,z) = ‾ [χ_KS,GG'^+-(q,z) + χ_KS,-G'-G^-+(-q,-z*)].
                              2
        """
        # Initiate new output file, if supplied
        if txt is not None:
            self.context.new_txt_and_timer(txt)

        zd = ComplexFrequencyDescriptor.from_array([0. + 0.j])
        chiks = self.chiks_calc.calculate('+-', q_c, zd)
        if np.allclose(q_c, 0.):
            chiks.symmetrize_reciprocity()

        # Take the reactive part
        chiksr = chiks.copy_reactive_part()

        return chiksr


def calculate_site_magnetization(
        gs: ResponseGroundStateAdapter | GPWFilename,
        sites: AtomicSites):
    """Calculate the site magnetization.

    Returns
    -------
    magmom_ap : np.ndarray
        Magnetic moment in μB of site a under partitioning p, calculated
        directly from the ground state density.
    """
    gs = ensure_gs(gs)
    site_data = AtomicSiteData(gs, sites)
    return site_data.calculate_magnetic_moments()


def calculate_site_zeeman_energy(
        gs: ResponseGroundStateAdapter | GPWFilename,
        sites: AtomicSites):
    """Calculate the site Zeeman energy.

    Returns
    -------
    EZ_ap : np.ndarray
        Local Zeeman energy in eV of site a under partitioning p, calculated
        directly from the ground state density.
    """
    gs = ensure_gs(gs)
    site_data = AtomicSiteData(gs, sites)
    return site_data.calculate_zeeman_energies() * Hartree  # Ha -> eV


def calculate_single_particle_site_magnetization(
        gs: ResponseGroundStateAdapter | GPWFilename,
        sites: AtomicSites,
        context: ResponseContext | TXTFilename = '-'):
    """Calculate the single-particle site magnetization.

    Returns
    -------
    sp_magmom_ap : np.ndarray
        Magnetic moment in μB of site a under partitioning p, calculated based
        on a single-particle sum rule.
    """
    gs, context = ensure_gs_and_context(gs, context=context)
    single_particle_calc = SingleParticleSiteMagnetizationCalculator(
        gs, sites, context=context)
    site_magnetization = single_particle_calc()
    return site_magnetization.array


def calculate_single_particle_site_zeeman_energy(
        gs: ResponseGroundStateAdapter | GPWFilename,
        sites: AtomicSites,
        context: ResponseContext | TXTFilename = '-'):
    """Calculate the single-particle site Zeeman energy.

    Returns
    -------
    sp_EZ_ap : np.ndarray
        Local Zeeman energy in eV of site a under partitioning p, calculated
        based on a single-particle sum rule.
    """
    gs, context = ensure_gs_and_context(gs, context=context)
    single_particle_calc = SingleParticleSiteZeemanEnergyCalculator(
        gs, sites, context=context)
    site_zeeman_energy = single_particle_calc()
    return site_zeeman_energy.array * Hartree  # Ha -> eV


def calculate_pair_site_magnetization(
        gs: ResponseGroundStateAdapter | GPWFilename,
        sites: AtomicSites,
        context: ResponseContext | TXTFilename = '-',
        q_c=[0., 0., 0.],
        nbands: int | None = None):
    """Calculate the pair site magnetization.

    Parameters
    ----------
    q_c : Vector
        q-vector to evaluate the pair site magnetization for.
    nbands : int or None
        Number of bands to include in the band summation of the pair site
        magnetization. If nbands is None, it includes all bands.

    Returns
    -------
    magmom_abp : np.ndarray
        Pair magnetization in μB of site a and b under partitioning p,
        calculated based on a two-particle sum rule.
    """
    gs, context = ensure_gs_and_context(gs, context=context)
    two_particle_calc = TwoParticleSiteMagnetizationCalculator(
        gs, sites, context=context, nbands=nbands)
    pair_site_magnetization = two_particle_calc(q_c)
    return pair_site_magnetization.array


def calculate_pair_site_zeeman_energy(
        gs: ResponseGroundStateAdapter | GPWFilename,
        sites: AtomicSites,
        context: ResponseContext | TXTFilename = '-',
        q_c=[0., 0., 0.],
        nbands: int | None = None):
    """Calculate the pair site Zeeman energy.

    Parameters
    ----------
    q_c : Vector
        q-vector to evaluate the pair site Zeeman energy for.
    nbands : int or None
        Number of bands to include in the band summation of the pair site
        Zeeman energy. If nbands is None, it includes all bands.

    Returns
    -------
    EZ_abp : np.ndarray
        Local pair Zeeman energy in eV of site a and b under partitioning p,
        calculated based on a two-particle sum rule.
    """
    gs, context = ensure_gs_and_context(gs, context=context)
    two_particle_calc = TwoParticleSiteZeemanEnergyCalculator(
        gs, sites, context=context, nbands=nbands)
    pair_site_zeeman_energy = two_particle_calc(q_c)
    return pair_site_zeeman_energy.array * Hartree  # Ha -> eV


class StaticSiteFunction(PairFunction):
    """Data object for static single-particle site functions."""
    def __init__(self,
                 qpd: SingleQPWDescriptor,
                 sites: AtomicSites):
        self.qpd = qpd
        self.q_c = qpd.q_c
        self.sites = sites
        self.array = self.zeros()

    @property
    def shape(self):
        return self.sites.shape

    def zeros(self):
        return np.zeros(self.shape, dtype=float)


class SingleParticleSiteSumRuleCalculator(PairFunctionIntegrator):
    r"""Calculator for single-particle site sum rules.

    For any site matrix element f^a_(nks,n'k's'), one may define a single-
    particle site sum rule by considering only the diagonal of the matrix
    element:
                 __  __
             1   \   \
    f_a^μ = ‾‾‾  /   /  σ^μ_ss f_nks f^a_(nks,nks)
            N_k  ‾‾  ‾‾
                 k   n,s

    where μ∊{0,z}.
    """

    def __init__(self, gs, sites, context):
        super().__init__(gs, context,
                         symmetry_analyzer=SymmetryAnalyzer(
                             point_group=False, time_reversal=False))

        # Set up calculator for the f^a matrix element
        self.sites = sites
        self.matrix_element_calc = self.create_matrix_element_calculator()

    @abstractmethod
    def create_matrix_element_calculator(self):
        """Create the desired site matrix element calculator."""

    def __call__(self):
        # Set up transitions
        # Loop over bands, which are fully or partially occupied
        nocc2 = self.kptpair_extractor.nocc2
        n_n = list(range(nocc2))
        n_t = np.array(n_n + n_n)
        s_t = np.array([0] * nocc2 + [1] * nocc2)
        transitions = PairTransitions(n1_t=n_t, n2_t=n_t, s1_t=s_t, s2_t=s_t)

        # Set up data object with q=0
        qpd = self.get_pw_descriptor([0., 0., 0.], ecut=1e-3)
        site_function = StaticSiteFunction(qpd, self.sites)

        # Perform actual calculation
        self._integrate(site_function, transitions)

        return site_function

    def add_integrand(self, kptpair, weight, site_function):
        r"""Add the integrand of the outer k-point integral.

        With
                   __
                1  \
        f_a^μ = ‾  /  (...)_k
                V  ‾‾
                   k

        the integrand has to be multiplied with the cell volume V0:
                     __
                     \
        (...)_k = V0 /  σ^μ_ss f_nks f^a_(nks,nks)
                     ‾‾
                     n,s
        """
        # Calculate matrix elements
        site_matrix_element = self.matrix_element_calc(
            kptpair, site_function.qpd)
        assert site_matrix_element.tblocks.blockcomm.size == 1
        f_tap = site_matrix_element.get_global_array()

        # Since we only use diagonal site matrix elements, corresponding
        # to the expectation value of the real functions Θ(r∊Ω_ap) and f(r),
        # f^a_(nks,nks) = <ψ_nks|Θ(r∊Ω_ap)f(r)|ψ_nks>,
        # the matrix elements are real
        assert np.allclose(f_tap.imag, 0.)
        f_tap = f_tap.real

        # Calculate Pauli matrix factors and multiply the occupations
        sigma = self.get_pauli_matrix()
        sigma_t = sigma[kptpair.transitions.s1_t, kptpair.transitions.s2_t]
        f_t = kptpair.get_all(kptpair.ikpt1.f_myt)
        sigmaf_t = sigma_t * f_t

        # Calculate and add integrand
        site_function.array[:] += self.gs.volume * weight * np.einsum(
            't, tap -> ap', sigmaf_t, f_tap)

    @abstractmethod
    def get_pauli_matrix(self):
        """Get the desired Pauli matrix σ^μ_ss."""


class SingleParticleSiteMagnetizationCalculator(
        SingleParticleSiteSumRuleCalculator):
    r"""Calculator for the single-particle site magnetization sum rule.

    The site magnetization is calculated from the site pair density:
                 __  __
             1   \   \
    n_a^z = ‾‾‾  /   /  σ^z_ss f_nks n^a_(nks,nks)
            N_k  ‾‾  ‾‾
                 k   n,s
    """
    def create_matrix_element_calculator(self):
        return SitePairDensityCalculator(self.gs, self.context, self.sites)

    def get_pauli_matrix(self):
        return smat('z')


class SingleParticleSiteZeemanEnergyCalculator(
        SingleParticleSiteMagnetizationCalculator):
    r"""Calculator for the single-particle site Zeeman energy sum rule.
                 __  __
             1   \   \
    E_a^Z = ‾‾‾  /   /  σ^z_ss f_nks E^(Z,a)_(nks,nks)
            N_k  ‾‾  ‾‾
                 k   n,s
    """
    def create_matrix_element_calculator(self):
        return SiteZeemanPairEnergyCalculator(
            self.gs, self.context, self.sites, rshewmin=1e-8)


class StaticSitePairFunction(StaticSiteFunction):
    """Data object for static pair site functions."""
    @property
    def shape(self):
        nsites = len(self.sites)
        npartitions = self.sites.npartitions
        return nsites, nsites, npartitions

    def zeros(self):
        return np.zeros(self.shape, dtype=complex)


class TwoParticleSiteSumRuleCalculator(PairFunctionIntegrator):
    r"""Calculator for two-particle site sum rules.

    For any set of site matrix elements f^a and g^b, one may define a two-
    particle site sum rule based on the lattice Fourier transformed quantity:
                     __  __   __
                 1   \   \    \   /
    ̄x_ab^z(q) = ‾‾‾  /   /    /   | σ^j_ss' (f_nks - f_n'k+qs')
                N_k  ‾‾  ‾‾   ‾‾  \                                       \
                     k  n,n' s,s'   × f^a_(nks,n'k+qs') g^b_(n'k+qs',nks) |
                                                                          /
    where σ^j is a Pauli matrix with j∊{0,+,-,z}.
    """

    def __init__(self,
                 gs: ResponseGroundStateAdapter,
                 sites: AtomicSites,
                 context: ResponseContext | None = None,
                 nbands: int | None = None):
        """Construct the two-particle site sum rule calculator."""
        if context is None:
            context = ResponseContext()
        super().__init__(gs, context,
                         symmetry_analyzer=SymmetryAnalyzer(
                             point_group=False, time_reversal=False))
        self.nbands = nbands

        # Set up calculators for the f^a and g^b matrix elements
        self.sites = sites
        mecalc1, mecalc2 = self.create_matrix_element_calculators()
        self.matrix_element_calc1 = mecalc1
        self.matrix_element_calc2 = mecalc2

    @abstractmethod
    def create_matrix_element_calculators(self):
        """Create the desired site matrix element calculators."""

    def __call__(self, q_c):
        """Calculate the site sum rule for a given wave vector q_c."""
        spincomponent = self.get_spincomponent()
        transitions = self.get_band_and_spin_transitions(
            spincomponent, nbands=self.nbands, bandsummation='double')
        self.context.print(self.get_info_string(
            q_c, self.nbands, len(transitions)))

        # Set up data object (without a plane-wave representation, which is
        # irrelevant in this case)
        qpd = self.get_pw_descriptor(q_c, ecut=1e-3)
        site_pair_function = StaticSitePairFunction(qpd, self.sites)

        # Perform actual calculation
        self._integrate(site_pair_function, transitions)

        return site_pair_function

    @abstractmethod
    def get_spincomponent(self):
        """Define how to rotate the spins via the spin component (μν)."""

    def add_integrand(self, kptpair, weight, site_pair_function):
        r"""Add the site sum rule integrand of the outer k-point integral.

        With
                       __
                    1  \
        ̄x_ab^z(q) = ‾  /  (...)_k
                    V  ‾‾
                       k

        the integrand is given by
                     __   __
                     \    \   /
        (...)_k = V0 /    /   | σ^j_ss' (f_nks - f_n'k+qs')
                     ‾‾   ‾‾  \                                       \
                    n,n' s,s'   × f^a_(nks,n'k+qs') g^b_(n'k+qs',nks) |
                                                                      /

        where V0 is the cell volume.
        """
        # Calculate site matrix elements
        qpd = site_pair_function.qpd
        matrix_element1 = self.matrix_element_calc1(kptpair, qpd)
        if self.matrix_element_calc2 is self.matrix_element_calc1:
            matrix_element2 = matrix_element1
        else:
            matrix_element2 = self.matrix_element_calc2(kptpair, qpd)

        # Calculate the product between the Pauli matrix and the occupational
        # differences
        sigma = self.get_pauli_matrix()
        s1_myt, s2_myt = kptpair.get_local_spin_indices()
        sigma_myt = sigma[s1_myt, s2_myt]
        df_myt = kptpair.ikpt1.f_myt - kptpair.ikpt2.f_myt
        sigmadf_myt = sigma_myt * df_myt

        # Calculate integrand
        f_mytap = matrix_element1.local_array_view
        g_mytap = matrix_element2.local_array_view
        fgcc_mytabp = f_mytap[:, :, np.newaxis] * g_mytap.conj()[:, np.newaxis]
        # Sum over local transitions
        integrand_abp = np.einsum('t, tabp -> abp', sigmadf_myt, fgcc_mytabp)
        # Sum over distributed transitions
        kptpair.tblocks.blockcomm.sum(integrand_abp)

        # Add integrand to output array
        site_pair_function.array[:] += self.gs.volume * weight * integrand_abp

    @abstractmethod
    def get_pauli_matrix(self):
        """Get the desired Pauli matrix σ^j_ss'."""

    def get_info_string(self, q_c, nbands, nt):
        """Get information about the calculation"""
        info_list = ['',
                     'Calculating two-particle site sum rule with:'
                     f'    q_c: [{q_c[0]}, {q_c[1]}, {q_c[2]}]',
                     self.get_band_and_transitions_info_string(nbands, nt),
                     '',
                     self.get_basic_info_string()]
        return '\n'.join(info_list)


class TwoParticleSiteMagnetizationCalculator(TwoParticleSiteSumRuleCalculator):
    r"""Calculator for the two-particle site magnetization sum rule.

    The site magnetization can be calculated from the site pair densities via
    the following sum rule [publication in preparation]:
                     __  __
                 1   \   \
    ̄n_ab^z(q) = ‾‾‾  /   /  (f_nk↑ - f_mk+q↓) n^a_(nk↑,mk+q↓) n^b_(mk+q↓,nk↑)
                N_k  ‾‾  ‾‾
                     k   n,m

              = δ_(a,b) n_a^z

    This is directly related to the sum rule of the χ^(+-) spin component of
    the four-component susceptibility tensor.
    """
    def create_matrix_element_calculators(self):
        site_pair_density_calc = SitePairDensityCalculator(
            self.gs, self.context, self.sites)
        return site_pair_density_calc, site_pair_density_calc

    def get_spincomponent(self):
        return '+-'

    def get_pauli_matrix(self):
        return smat('+')


class TwoParticleSiteZeemanEnergyCalculator(
        TwoParticleSiteMagnetizationCalculator):
    r"""Calculator for the two-particle site Zeeman energy sum rule.

    The site Zeeman energy can be calculated from the site pair density and
    site Zeeman pair energy via the following sum rule [publication in
    preparation]:
                     __  __
    ˍ            1   \   \  /
    E_ab^Z(q) = ‾‾‾  /   /  | (f_nk↑ - f_mk+q↓)
                N_k  ‾‾  ‾‾ \                                       \
                     k   n,m  × E^(Z,a)_(nk↑,mk+q↓) n^b_(mk+q↓,nk↑) |
                                                                    /
              = δ_(a,b) E_a^Z
    """
    def create_matrix_element_calculators(self):
        site_zeeman_pair_energy_calc = SiteZeemanPairEnergyCalculator(
            self.gs, self.context, self.sites, rshewmin=1e-8)
        site_pair_density_calc = SitePairDensityCalculator(
            self.gs, self.context, self.sites)
        return site_zeeman_pair_energy_calc, site_pair_density_calc
