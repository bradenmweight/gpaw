import numpy as np

from gpaw.response import timer
from gpaw.response.paw import get_pair_density_paw_corrections


class PlaneWavePairDensity:
    """Class for calculating pair densities

    n_kt(G+q) = n_nks,n'k+qs'(G+q) = <nks| e^-i(G+q)r |n'k+qs'>_V0

    for a single k-point pair (k,k+q) in the plane wave mode"""
    def __init__(self, kspair):
        self.gs = kspair.gs
        self.context = kspair.context
        self.transitionblockscomm = kspair.transitionblockscomm

        # Save PAW correction for all calls with same q_c
        self.pawcorr = None
        self.currentq_c = None

    def initialize(self, qpd):
        """Initialize PAW corrections ahead in time of integration."""
        self.initialize_paw_corrections(qpd)

    @timer('Initialize PAW corrections')
    def initialize_paw_corrections(self, qpd):
        """Initialize PAW corrections, if not done already, for the given q"""
        q_c = qpd.q_c
        if self.pawcorr is None or not np.allclose(q_c - self.currentq_c, 0.):
            self.pawcorr = self._initialize_paw_corrections(qpd)
            self.currentq_c = q_c

    def _initialize_paw_corrections(self, qpd):
        pawdatasets = self.gs.pawdatasets
        spos_ac = self.gs.spos_ac
        return get_pair_density_paw_corrections(pawdatasets, qpd, spos_ac)

    @timer('Calculate pair density')
    def __call__(self, kskptpair, qpd):
        """Calculate the pair densities for all transitions t of the (k,k+q)
        k-point pair:

        n_kt(G+q) = <nks| e^-i(G+q)r |n'k+qs'>_V0

                    /
                  = | dr e^-i(G+q)r psi_nks^*(r) psi_n'k+qs'(r)
                    /V0
        """
        Q_aGii = self.get_paw_projectors(qpd)
        Q_G = self.get_fft_indices(kskptpair, qpd)
        mynt, nt, ta, tb = kskptpair.transition_distribution()

        n_mytG = qpd.empty(mynt)

        # Calculate smooth part of the pair densities:
        with self.context.timer('Calculate smooth part'):
            ut1cc_mytR = kskptpair.kpt1.ut_tR.conj()
            n_mytR = ut1cc_mytR * kskptpair.kpt2.ut_tR
            # Unvectorized
            for myt in range(tb - ta):
                n_mytG[myt] = qpd.fft(n_mytR[myt], 0, Q_G) * qpd.gd.dv

        # Calculate PAW corrections with numpy
        with self.context.timer('PAW corrections'):
            P1 = kskptpair.kpt1.projections
            P2 = kskptpair.kpt2.projections
            for (Q_Gii, (a1, P1_myti),
                 (a2, P2_myti)) in zip(Q_aGii, P1.items(), P2.items()):
                P1cc_myti = P1_myti[:tb - ta].conj()
                C1_Gimyt = np.tensordot(Q_Gii, P1cc_myti, axes=([1, 1]))
                P2_imyt = P2_myti.T[:, :tb - ta]
                n_mytG[:tb - ta] += np.sum(C1_Gimyt * P2_imyt[np.newaxis,
                                                              :, :], axis=1).T

        # Attach the calculated pair density to the KohnShamKPointPair object
        kskptpair.attach('n_mytG', 'n_tG', n_mytG)

    def get_paw_projectors(self, qpd):
        """Make sure PAW correction has been initialized properly
        and return projectors"""
        self.initialize_paw_corrections(qpd)
        return self.pawcorr.Q_aGii

    @timer('Get G-vector indices')
    def get_fft_indices(self, kskptpair, qpd):
        """Get indices for G-vectors inside cutoff sphere."""
        from gpaw.response.pair import fft_indices

        kpt1 = kskptpair.kpt1
        kpt2 = kskptpair.kpt2
        kd = self.gs.kd
        q_c = qpd.q_c

        return fft_indices(kd=kd, K1=kpt1.K, K2=kpt2.K, q_c=q_c, qpd=qpd,
                           shift0_c=kpt1.shift_c - kpt2.shift_c)
