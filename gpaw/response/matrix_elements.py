from functools import partial

import numpy as np

from gpaw.response import timer
from gpaw.response.kspair import KohnShamKPointPair


class PairDensity:
    """Data class for transition distributed pair density arrays."""

    def __init__(self, tblocks, n_mytG):
        self.tblocks = tblocks
        self.n_mytG = n_mytG

    @classmethod
    def from_qpd(cls, tblocks, qpd):
        n_mytG = qpd.zeros(tblocks.blocksize)
        return cls(tblocks, n_mytG)

    @property
    def local_array_view(self):
        return self.n_mytG[:self.tblocks.nlocal]

    def get_global_array(self):
        """Get the global (all gathered) pair density array n_tG."""
        n_tG = self.tblocks.collect(self.n_mytG)

        return n_tG


class NewPairDensityCalculator:
    r"""Class for calculating pair densities

    n_kt(G+q) = n_nks,n'k+qs'(G+q) = <nks| e^-i(G+q)r |n'k+qs'>_V0

                /
              = | dr e^-i(G+q)r ψ_nks^*(r) ψ_n'k+qs'(r)
                /V0

    for a single k-point pair (k, k + q) in the plane-wave mode."""
    def __init__(self, gs, context):
        self.gs = gs
        self.context = context

        # Save PAW correction for all calls with same q_c
        self._pawcorr = None
        self._currentq_c = None

    def initialize_paw_corrections(self, qpd):
        """Initialize the PAW corrections ahead of the actual calculation."""
        self.get_paw_corrections(qpd)

    def get_paw_corrections(self, qpd):
        """Get PAW corrections correcsponding to a specific q-vector."""
        if self._pawcorr is None \
           or not np.allclose(qpd.q_c - self._currentq_c, 0.):
            with self.context.timer('Initialize PAW corrections'):
                self._pawcorr = self.gs.pair_density_paw_corrections(qpd)
                self._currentq_c = qpd.q_c

        return self._pawcorr

    @timer('Calculate pair density')
    def __call__(self, kptpair: KohnShamKPointPair, qpd) -> PairDensity:
        r"""Calculate the pair density for all transitions t.

        In the PAW method, the all-electron pair density is calculated in
        two contributions, the pseudo pair density and a PAW correction,

        n_kt(G+q) = ñ_kt(G+q) + Δn_kt(G+q),

        see [PRB 103, 245110 (2021)] for details.
        """
        # Construct symmetrizers for the periodic part of the pseudo waves
        # and for the PAW projectors
        ut1_symmetrizer, Ph1_symmetrizer, shift1_c = \
            self.construct_symmetrizers(kptpair.kpt1)
        ut2_symmetrizer, Ph2_symmetrizer, shift2_c = \
            self.construct_symmetrizers(kptpair.kpt2)

        # Initialize a blank pair density object
        pair_density = PairDensity.from_qpd(kptpair.tblocks, qpd)
        n_mytG = pair_density.local_array_view

        self.add_pseudo_pair_density(kptpair, qpd, n_mytG,
                                     ut1_symmetrizer, ut2_symmetrizer,
                                     shift1_c, shift2_c)
        self.add_paw_correction(kptpair, qpd, n_mytG,
                                Ph1_symmetrizer, Ph2_symmetrizer)

        return pair_density

    @timer('Calculate the pseudo pair density')
    def add_pseudo_pair_density(self, kptpair, qpd, n_mytG,
                                ut1_symmetrizer, ut2_symmetrizer,
                                shift1_c, shift2_c):
        r"""Add the pseudo pair density to an output array.

        The pseudo pair density is first evaluated on the coarse real-space
        grid and then FFT'ed to reciprocal space:

                    /               ˷          ˷
        ñ_kt(G+q) = | dr e^-i(G+q)r ψ_nks^*(r) ψ_n'k+qs'(r)
                    /V0
                                 ˷          ˷
                  = FFT_G[e^-iqr ψ_nks^*(r) ψ_n'k+qs'(r)]
        """
        kpt1 = kptpair.kpt1
        kpt2 = kptpair.kpt2
        # Fourier transform the periodic part of the pseudo waves to the coarse
        # real-space grid and symmetrize them
        ut1_hR = self.get_periodic_pseudo_waves(kpt1, ut1_symmetrizer)
        ut2_hR = self.get_periodic_pseudo_waves(kpt2, ut2_symmetrizer)

        # Calculate the pseudo pair density in real space
        ut1cc_mytR = ut1_hR[kpt1.h_myt].conj()
        n_mytR = ut1cc_mytR * ut2_hR[kpt2.h_myt]

        # Get the plane-wave indices to Fourier transform products of
        # Kohn-Sham orbitals in k and k + q
        dshift_c = shift1_c - shift2_c
        Q_G = self.get_fft_indices(kpt1.K, kpt2.K, qpd, dshift_c)

        # Add FFT of the pseudo pair density to the output array
        nlocalt = kptpair.tblocks.nlocal
        assert len(n_mytG) == nlocalt and len(n_mytR) == nlocalt
        for n_G, n_R in zip(n_mytG, n_mytR):
            n_G[:] += qpd.fft(n_R, 0, Q_G) * qpd.gd.dv

    @timer('Calculate the pair density PAW corrections')
    def add_paw_correction(self, kptpair, qpd, n_mytG,
                           Ph1_symmetrizer, Ph2_symmetrizer):
        r"""Add the pair-density PAW correction to the output array.

        The correction is calculated as a sum over augmentation spheres a
        and projector indices i and j,
                     __  __
                     \   \   ˷     ˷     ˷    ˷
        Δn_kt(G+q) = /   /  <ψ_nks|p_ai><p_aj|ψ_n'k+qs'> Q_aij(G+q)
                     ‾‾  ‾‾
                     a   i,j

        where the pair-density PAW correction tensor is calculated from the
        smooth and all-electron partial waves:

                     /
        Q_aij(G+q) = | dr e^-i(G+q)r [φ_ai^*(r-R_a) φ_aj(r-R_a)
                     /V0                ˷             ˷
                                      - φ_ai^*(r-R_a) φ_aj(r-R_a)]
        """
        kpt1 = kptpair.kpt1
        kpt2 = kptpair.kpt2

        # Symmetrize the projectors
        P1h = Ph1_symmetrizer(kpt1.Ph)
        P2h = Ph2_symmetrizer(kpt2.Ph)

        # Calculate the actual PAW corrections
        Q_aGii = self.get_paw_corrections(qpd).Q_aGii
        P1 = kpt1.projectors_in_transition_index(P1h)
        P2 = kpt2.projectors_in_transition_index(P2h)
        for a, Q_Gii in enumerate(Q_aGii):  # Loop over augmentation spheres
            # NB: There does not seem to be any strict guarantee that the order
            # of the PAW corrections matches the projections keys.
            # This is super dangerous and should be rectified in the future XXX
            P1_myti = P1[a]
            P2_myti = P2[a]
            # Make outer product of the projectors in the projector index i,j
            P1ccP2_mytii = P1_myti.conj()[..., np.newaxis] \
                * P2_myti[:, np.newaxis]
            # Sum over projector indices and add correction to the output
            n_mytG[:] += np.einsum('tij, Gij -> tG', P1ccP2_mytii, Q_Gii)

    def get_periodic_pseudo_waves(self, kpt, ut_symmetrizer):
        """FFT the Kohn-Sham orbitals to real space and symmetrize them."""
        ik = self.gs.kd.bz2ibz_k[kpt.K]
        ut_hR = self.gs.gd.empty(kpt.nh, self.gs.dtype)
        for h, psit_G in enumerate(kpt.psit_hG):
            ut_hR[h] = ut_symmetrizer(self.gs.global_pd.ifft(psit_G, ik))

        return ut_hR

    def construct_symmetrizers(self, kpt):
        """Construct functions to symmetrize ut_hR and Ph."""
        _, T, a_a, U_aii, shift_c, time_reversal = \
            self.gs.construct_symmetry_operators(kpt.K, kpt.k_c)

        ut_symmetrizer = T
        Ph_symmetrizer = partial(symmetrize_projections,
                                 a1_a2=a_a, U_aii=U_aii,
                                 time_reversal=time_reversal)

        return ut_symmetrizer, Ph_symmetrizer, shift_c

    def get_fft_indices(self, K1, K2, qpd, dshift_c):
        from gpaw.response.pair import fft_indices
        return fft_indices(self.gs.kd, K1, K2, qpd, dshift_c)


def symmetrize_projections(Ph, a1_a2, U_aii, time_reversal):
    """Symmetrize the PAW projections.

    NB: The projections of atom a1 are mapped onto a *different* atom a2
    according to the input map of atomic indices a1_a2."""
    # First, we apply the symmetry operations to the projections one at a time
    P_a2hi = []
    for a1, U_ii in zip(a1_a2, U_aii):
        P_hi = Ph[a1].copy(order='C')
        np.dot(P_hi, U_ii, out=P_hi)
        if time_reversal:
            np.conj(P_hi, out=P_hi)
        P_a2hi.append(P_hi)

    # Then, we store the symmetry mapped projectors in the projections object
    for a2, P_hi in enumerate(P_a2hi):
        I1, I2 = Ph.map[a2]
        Ph.array[..., I1:I2] = P_hi

    return Ph