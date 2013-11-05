"""Todo: optical limit, Hilbert transform"""

import sys
from math import pi

import numpy as np
from ase.units import Hartree

from gpaw.mpi import world
from gpaw.utilities.blas import gemm, rk
from gpaw.wavefunctions.pw import PWDescriptor
from gpaw.kpt_descriptor import KPointDescriptor
from gpaw.response.math_func import two_phi_planewave_integrals


class Chi0:
    def __init__(self, calc, omega_w, ecut=50 / Hartree, hilbert=False,
                 eta=0.2 / Hartree, blocksize=50, spin=0, ftol=1e-6,
                 comm=world, fd=None):
        self.calc = calc
        self.omega_w = omega_w
        self.ecut = ecut
        self.hilbert = hilbert
        self.eta = eta
        self.blocksize = blocksize
        self.spin = spin
        self.ftol = ftol
        self.comm = comm
        self.fd = fd or sys.stdout
            
        self.spos_ac = calc.atoms.get_scaled_positions()
        
        self.nocc1 = None  # number of completely filled bands
        self.nocc2 = None  # number of non-empty bands
        self.count_occupied_bands()
        
        self.Kn1n2 = None  # my occupied k-points and bands
        self.distribute_k_points_and_bands()
        
        self.ut_knR = None  # periodic part of wave functions in real-space
        self.eps_kn = None  # eigenvalues
        self.f_kn = None  # occupation numbers
        self.P_kani = None  # PAW projections
        self.shift_kc = None  # long story - see construct_symmetry_operators()
        self.initialize_occupied_states()
        
        vol = abs(np.linalg.det(calc.wfs.gd.cell_cv))
        self.prefactor = 2 / vol / calc.wfs.kd.nbzkpts

    def count_occupied_bands(self):
        self.nocc1 = 9999999
        self.nocc2 = 0
        for kpt in self.calc.wfs.kpt_u:
            f_n = kpt.f_n / kpt.weight
            self.nocc1 = min((f_n > 1 - self.ftol).sum(), self.nocc1)
            self.nocc2 = max((f_n > self.ftol).sum(), self.nocc2)
        self.fd.write('Number of completely filled bands: %d\n' % self.nocc1)
        self.fd.write('Number of non-empty bands: %d\n' % self.nocc2)
        self.fd.write('Total number of bands: %d\n' % self.calc.wfs.bd.nbands)
        
    def distribute_k_points_and_bands(self):
        comm = self.comm
        nk = self.calc.wfs.kd.nbzkpts
        n = (nk * self.nocc2 + comm.size - 1) // comm.size
        i1 = comm.rank * n
        i2 = min((comm.rank + 1) * n, nk * self.nocc2)
        K1, n1 = divmod(i1, self.nocc2)
        K2, n2 = divmod(i2, self.nocc2)
        self.Kn1n2_k = []
        for K in range(K1, K2):
            self.Kn1n2_k.append((K, n1, self.nocc2))
            n1 = 0
        if n2 > 0:
            self.Kn1n2_k.append((K2, n1, n2))
        self.fd.write('k-points: %s\n' % self.calc.wfs.kd.description)
        self.fd.write('Distributing %d x %d bands over %d process%s' %
                      (nk, self.nocc2, comm.size, ['es', ''][comm.size == 1]))
            
    def initialize_occupied_states(self):
        self.ut_knR = []
        self.eps_kn = []
        self.f_kn = []
        self.P_kani = []
        self.shift_kc = []
        for K, n1, n2 in self.Kn1n2_k:
            ut_nG, eps_n, f_n, P_ani, shift_c = self.get_k_point(K, n1, n2)
            self.ut_knR.append(ut_nG)
            self.eps_kn.append(eps_n)
            self.f_kn.append(f_n)
            self.P_kani.append(P_ani)
            self.shift_kc.append(shift_c)

    def get_k_point(self, K, n1, n2):
        wfs = self.calc.wfs
        
        T, T_a, shift_c = self.construct_symmetry_operators(K)
        ik = wfs.kd.bz2ibz_k[K]
        kpt = wfs.kpt_u[ik]
        
        psit_nG = kpt.psit_nG
        ut_nR = wfs.gd.empty(n2 - n1, complex)
        for n in range(n1, n2):
            ut_nR[n - n1] = T(wfs.pd.ifft(psit_nG[n], ik))

        eps_n = kpt.eps_n[n1:n2]
        f_n = kpt.f_n[n1:n2] / kpt.weight
        
        P_ani = []
        for (b, T_ii, time_reversal) in T_a:
            P_ni = np.dot(kpt.P_ani[b][n1:n2], T_ii)
            if time_reversal:
                P_ni = P_ni.conj()
            P_ani.append(P_ni)
        
        return ut_nR, eps_n, f_n, P_ani, shift_c
    
    def calculate(self, q_c, chi0_wGG=None, nocc1=None, nbands=None,
                  direction=None):
        wfs = self.calc.wfs

        if self.eta == 0.0:
            update = self.update_hermetian
        elif self.hilbert:
            update = self.update_hilbert
        else:
            update = self.update
            
        qd = KPointDescriptor([q_c])
        pd = PWDescriptor(self.ecut, wfs.gd, complex, qd)
        nG = pd.ngmax
        
        if chi0_wGG is None:
            # Start from scratch and do all empty bands:
            chi0_wGG = np.zeros((len(self.omega_w), nG, nG), complex)
            nocc1 = self.nocc1
            nbands = wfs.bd.nbands
            
        Q_aGii = self.calculate_paw_corrections(pd)
        
        for k, (K, n1, n2) in enumerate(self.Kn1n2_k):
            P_ani = self.P_kani[k]
            K2 = wfs.kd.find_k_plus_q(q_c, [K])[0]
            ut_mR, eps_m, f_m, P_ami, shift_c = self.get_k_point(
                K2, nocc1, nbands)
            Q_G = self.get_fft_indices(K, K2, q_c, pd,
                                       self.shift_kc[k] - shift_c)
            for n in range(n2 - n1):
                eps = self.eps_kn[k][n]
                f = self.f_kn[k][n]
                utcc_R = self.ut_knR[k][n].conj()
                C_aGi = [np.dot(Q_Gii, P_ni[n].conj())
                         for Q_Gii, P_ni in zip(Q_aGii, P_ani)]
                for m1 in range(0, nbands - nocc1, self.blocksize):
                    m2 = min(m1 + self.blocksize, nbands - nocc1)
                    n_mG = self.calculate_pair_densities(utcc_R, C_aGi,
                                                         ut_mR, P_ami,
                                                         m1, m2, pd, Q_G)
                    deps_m = eps - eps_m[m1:m2]
                    update(n_mG, deps_m, f - f_m, chi0_wGG)
                    
        world.sum(chi0_wGG)
        
        if self.eta == 0.0:
            # Set lower triangle also:
            il = np.tril_indices(nG, -1)
            iu = il[::-1]
            for chi0_GG in chi0_wGG:
                chi0_GG[il] = chi0_GG[iu]
                
        elif self.hilbert:
            for G in range(nG):
                chi0_wGG[:, :, G] = np.dot(A_ww, chi0_wGG[:, :, G])
        
        return chi0_wGG, pd
        
    def update(self, n_mG, deps_m, df_m, chi0_wGG):
        for w, omega in enumerate(self.omega_w):
            x_m = df_m * (1.0 / (omega + deps_m + 1j * self.eta) -
                          1.0 / (omega - deps_m + 1j * self.eta))
            x_mG = n_mG * x_m[:, np.newaxis]
            gemm(self.prefactor, n_mG.conj(), np.ascontiguousarray(x_mG.T),
                 1.0, chi0_wGG[w])

    def update_hermitian(self, n_mG, deps_m, df_m, chi0_wGG):
        for w, omega in enumerate(self.omega_w):
            x_m = (2 * df_m * deps_m / (omega.imag**2 + deps_m**2))**0.5
            x_mG = n_mG * x_m[:, np.newaxis]
            rk(self.prefactor, x_mG, 1.0, chi0_wGG[w])

    def update_hilbert(self, n_mG, deps_m, df_m, chi0_wGG):
        domega = self.omega_w[1]
        for omega, df, n_G in zip(deps_m, df_m, n_mG):
            w = omega / domega
            iw = int(w)
            weights = df * np.array([[1 - w + iw], [w - iw]])
            x_2G = n_G * weights**0.5
            rk(self.prefactor, x_2G, 1.0, chi0_wGG[iw:iw + 2])

    def calculate_pair_densities(self, utcc_R, C_aGi, ut_mR, P_ami,
                                 m1, m2, pd, Q_G):
        dv = pd.gd.dv
        n_mG = pd.empty(m2 - m1)
        for m in range(m1, m2):
            n_R = utcc_R * ut_mR[m]
            pd.tmp_R[:] = n_R
            pd.fftplan.execute()
            n_mG[m - m1] = pd.tmp_Q.ravel()[Q_G] * dv
        
        # PAW corrections:
        for C_Gi, P_mi in zip(C_aGi, P_ami):
            gemm(1.0, C_Gi, P_mi[m1:m2], 1.0, n_mG, 't')
            
        return n_mG

    def get_fft_indices(self, K, K2, q_c, pd, shift0_c):
        kd = self.calc.wfs.kd
        Q_G = pd.Q_qG[0]
        shift_c = (shift0_c +
                   (q_c - kd.bzk_kc[K2] + kd.bzk_kc[K]).round().astype(int))
        if shift_c.any():
            q_cG = np.unravel_index(Q_G, pd.gd.N_c)
            q_cG = [q_G + shift for q_G, shift in zip(q_cG, shift_c)]
            Q_G = np.ravel_multi_index(q_cG, pd.gd.N_c, 'wrap')
        return Q_G
        
    def construct_symmetry_operators(self, K):
        """Construct symmetry operators for wave function and PAW projections.
        
        We want to transform a k-point in the irreducible part of the BZ to
        the corresponding k-point with index K.
        
        Returns T, T_a and shift_c, where:
            
        * T() is a function that transforms the periodic part of the wave
          function.
        * T_a is a list of (b, U_ii, time_reversal) tuples (one for each
          atom a), where:
        
          * b is the symmetry related atom index
          * U_ii is a rotation matrix for the PAW projections
          * time_reversal is a flag - if True, projections should be complex
            conjugated
            
          See the get_k_point() method for how tu use these tuples.
        
        * shift_c is three integers: see code below.
        """
        
        wfs = self.calc.wfs
        kd = wfs.kd

        s = kd.sym_k[K]
        U_cc = kd.symmetry.op_scc[s]
        time_reversal = kd.time_reversal_k[K]
        ik = kd.bz2ibz_k[K]
        k_c = kd.bzk_kc[K]
        ik_c = kd.ibzk_kc[ik]
        
        sign = 1 - 2 * time_reversal
        shift_c = np.dot(U_cc, ik_c) - k_c * sign
        assert np.allclose(shift_c.round(), shift_c)
        shift_c = shift_c.round().astype(int)
        
        if (U_cc == np.eye(3)).all():
            T = lambda f_R: f_R
        else:
            N_c = self.calc.wfs.gd.N_c
            i_cr = np.dot(U_cc.T, np.indices(N_c).reshape((3, -1)))
            i = np.ravel_multi_index(i_cr, N_c, 'wrap')
            T = lambda f_R: f_R.ravel()[i].reshape(N_c)

        if time_reversal:
            T0 = T
            T = lambda f_R: T0(f_R).conj()
            shift_c *= -1
        
        T_a = []
        for a, id in enumerate(self.calc.wfs.setups.id_a):
            b = kd.symmetry.a_sa[s, a]
            S_c = np.dot(self.spos_ac[a], U_cc) - self.spos_ac[b]
            x = np.exp(2j * pi * np.dot(ik_c, S_c))
            U_ii = wfs.setups[a].R_sii[s].T * x
            T_a.append((b, U_ii, time_reversal))

        return T, T_a, shift_c

    def calculate_paw_corrections(self, pd):
        wfs = self.calc.wfs
        q_v = pd.K_qv[0]
        G_Gv = pd.G_Qv[pd.Q_qG[0]] + q_v
        pos_av = np.dot(self.spos_ac, pd.gd.cell_cv)
          
        # Collect integrals for all species:
        Q_xGii = {}
        for id, atomdata in wfs.setups.setups.items():
            Q_Gii = two_phi_planewave_integrals(G_Gv, atomdata)
            ni = atomdata.ni
            Q_xGii[id] = Q_Gii.reshape((-1, ni, ni))
        
        Q_aGii = []
        for a, atomdata in enumerate(wfs.setups):
            id = wfs.setups.id_a[a]
            Q_Gii = Q_xGii[id]
            x_G = np.exp(-1j * np.dot(G_Gv, pos_av[a]))
            Q_aGii.append(x_G[:, np.newaxis, np.newaxis] * Q_Gii)
        return Q_aGii
