# General modules
import numpy as np

# GPAW modules
from gpaw.sphere.integrate import default_spherical_drcut

from gpaw.response.frequencies import ComplexFrequencyDescriptor
from gpaw.response.chiks import ChiKSCalculator
from gpaw.response.localft import LocalFTCalculator, add_LSDA_Bxc
from gpaw.response.site_kernels import SiteKernels

# ASE modules
from ase.units import Hartree, Bohr

from ase.neighborlist import natural_cutoffs, build_neighbor_list


class IsotropicExchangeCalculator:
    r"""Calculator class for the Heisenberg exchange constants

    _           2
    J^ab(q) = - ‾‾ B^(xc†) K^(a†)(q) χ_KS^('+-)(q) K^b(q) B^(xc)
                V0

    calculated for an isotropic system in a plane wave representation using
    the magnetic force theorem within second order perturbation theory, see
    [arXiv:2204.04169].

    Entering the formula for the isotropic exchange constant at wave vector q
    between sublattice a and b is the unit cell volume V0, the functional
    derivative of the (LDA) exchange-correlation energy with respect to the
    magnitude of the magnetization B^(xc), the sublattice site kernels K^a(q)
    and K^b(q) as well as the reactive part of the static transverse magnetic
    susceptibility of the Kohn-Sham system χ_KS^('+-)(q).

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
            assert getattr(chiks_calc, key) == item,\
                f'Expected chiks.{key} == {item}. '\
                f'Got: {getattr(chiks_calc, key)}'

        self.chiks_calc = chiks_calc
        self.context = chiks_calc.context

        # Check assumed properties of the LocalFTCalculator
        assert localft_calc.context is self.context
        assert localft_calc.gs is chiks_calc.gs
        self.localft_calc = localft_calc

        # Bxc field buffer
        self._Bxc_G = None

        # chiksr buffer
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
        Bxc_G = self.get_Bxc()
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
                    J = np.conj(Bxc_G) @ np.conj(K_aGG[a]).T @ chiksr_GG \
                        @ K_aGG[b] @ Bxc_G
                    J_ab[a, b] = - 2. * J / V0

        # Transpose to have the partitions index last
        J_abp = np.transpose(J_pab, (1, 2, 0))

        return J_abp * Hartree  # Convert from Hartree to eV

    def get_Bxc(self):
        """Get B^(xc)_G from buffer."""
        if self._Bxc_G is None:  # Calculate if buffer is empty
            self._Bxc_G = self._calculate_Bxc()

        return self._Bxc_G

    def _calculate_Bxc(self):
        """Use the PlaneWaveBxc calculator to calculate the plane wave
        coefficients B^xc_G"""
        # Create a plane wave descriptor encoding the plane wave basis. Input
        # q_c is arbitrary, since we are assuming that chiks.gammacentered == 1
        qpd0 = self.chiks_calc.get_pw_descriptor([0., 0., 0.])

        return self.localft_calc(qpd0, add_LSDA_Bxc)

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


class AtomicSiteData:
    r"""Data object for spherical atomic sites."""

    def __init__(self, gs, indices, radii):
        """
        Some documentation here! XXX
        """
        self.A_a = np.asarray(indices)

        # Parse the input atomic radii
        rc_ap = np.asarray(radii)
        assert rc_ap.ndim == 2
        assert rc_ap.shape[0] == len(self.A_a)
        # Convert radii to internal units (Å to Bohr)
        self.rc_ap = rc_ap / Bohr

        assert self._in_valid_site_radii_range(gs),\
            'Please provide site radii in the valid range, see '\
            'AtomicSiteData.valid_site_radii_range()'

    @staticmethod
    def _valid_site_radii_range(gs):
        """For each atom in gs, determine the valid site radii range in Bohr.

        The lower bound is determined by the spherical truncation width, when
        truncating integrals on the real-space grid.
        The upper bound is determined by the distance to the nearest
        augmentation sphere.
        """
        atoms = gs.atoms
        drcut = default_spherical_drcut(gs.gd)
        rmin_A = np.array([drcut / 2] * len(atoms))

        # Find neighbours based on covalent radii
        cutoffs = natural_cutoffs(atoms, mult=2)
        neighbourlist = build_neighbor_list(atoms, cutoffs)
        # Determine rmax for each atom
        augr_A = gs.get_aug_radii()
        rmax_A = []
        for A in range(len(atoms)):
            pos = atoms.positions[A]
            # Calculate the distance to the augmentation sphere of each
            # neighbour
            aug_distances = []
            for An, offset in zip(*neighbourlist.get_neighbors(A)):
                if An == A and np.all(offset == 0):
                    continue  # The atom itself is somehow a neighbour...
                posn = atoms.positions[An] + offset @ atoms.get_cell()
                dist = np.linalg.norm(posn - pos) / Bohr  # Å -> Bohr
                aug_dist = dist - augr_A[An]
                assert aug_dist > 0.
                aug_distances.append(aug_dist)
            # In order for PAW corrections to be valid, we need a sphere of
            # radius rcut not to overlap with any neighbouring augmentation
            # spheres
            rmax_A.append(min(aug_distances))
        rmax_A = np.array(rmax_A)

        return rmin_A, rmax_A

    @staticmethod
    def valid_site_radii_range(gs):
        """
        Some documentation here! XXX
        """
        rmin_A, rmax_A = AtomicSiteData._valid_site_radii_range(gs)
        # Convert to external units (Bohr to Å)
        return rmin_A * Bohr, rmax_A * Bohr

    def _in_valid_site_radii_range(self, gs):
        rmin_A, rmax_A = AtomicSiteData._valid_site_radii_range(gs)
        for a, A in enumerate(self.A_a):
            if not np.all(
                    np.logical_and(
                        self.rc_ap[a] > rmin_A[A] - 1e-8,
                        self.rc_ap[a] < rmax_A[A] + 1e-8)):
                return False
        return True
        
        
