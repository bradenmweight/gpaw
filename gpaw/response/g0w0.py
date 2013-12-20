from __future__ import division

import sys
from math import pi

import numpy as np
from ase.units import Hartree
from ase.utils import prnt, opencew
from ase.dft.kpoints import monkhorst_pack

import gpaw.mpi as mpi
from gpaw.response.chi0 import Chi0
from gpaw.response.pair import PairDensity
from gpaw.wavefunctions.pw import PWDescriptor
from gpaw.kpt_descriptor import KPointDescriptor
from gpaw.response.wstc import WignerSeitzTruncatedCoulomb


class G0W0(PairDensity):
    def __init__(self, calc, filename='gw',
                 kpts=None, bands=None, nbands=None, ppa=False, hilbert=False,
                 ecut=150.0, eta=0.1, E0=1.0 * Hartree,
                 domega=0.025, omegamax=100,
                 world=mpi.world, txt=sys.stdout):
    
        PairDensity.__init__(self, calc, ecut, world=world, txt=txt)
        
        self.filename = filename
        
        ecut /= Hartree
        
        self.ppa = ppa
        self.hilbert = hilbert
        self.eta = eta / Hartree
        self.E0 = E0 / Hartree
        self.domega = domega / Hartree
        self.omegamax = omegamax / Hartree
        
        if kpts is None:
            kpts = range(self.calc.wfs.kd.nibzkpts)
        
        if bands is None:
            bands = [0, self.nocc2]
            
        self.kpts = kpts
        self.bands = bands

        shape = (self.calc.wfs.nspins, len(kpts), bands[1] - bands[0])
        self.sigma_sin = np.zeros(shape)   # self-energies
        self.dsigma_sin = np.zeros(shape)  # derivatives of self-energies
        self.Z_sin = np.empty(shape)       # renormalization factors
        self.exx_sin = np.empty(shape)     # exact exchange contributions
        self.eps_sin = np.empty(shape)     # KS-eigenvalues

        if nbands is None:
            nbands = int(self.vol * ecut**1.5 * 2**0.5 / 3 / pi**2)
        self.nbands = nbands
        
        self.mysKn1n2 = None  # my (s, K, n1, n2) indices
        self.distribute_k_points_and_bands(nbands)
        
        self.omega_w = None  # frequencies
        self.initialize_frequencies()
        
        # Find q-vectors and weights in IBZ:
        kd = self.calc.wfs.kd
        assert -1 not in kd.bz2bz_ks
        offset_c = 0.5 * ((kd.N_c + 1) % 2) / kd.N_c
        bzq_qc = monkhorst_pack(kd.N_c) + offset_c
        #bzq_qc = np.array([(-0.5,0,0),(0,0,0),(0.5,0,0)])
        self.qd = KPointDescriptor(bzq_qc)
        #self.qd.set_symmetry(self.calc.atoms, self.calc.wfs.setups,
        #                     usesymm=True, N_c=self.calc.wfs.gd.N_c)
        
    def initialize_frequencies(self, domega=0.05):
        domega /= Hartree
        epsmin = 10000.0
        epsmax = -10000.0
        for kpt in self.calc.wfs.kpt_u:
            epsmin = min(epsmin, kpt.eps_n[0])
            epsmax = max(epsmax, kpt.eps_n[self.nbands - 1])
        
        self.omega_w = np.linspace(0, self.omegamax,
                                   round(self.omegamax / self.domega) + 1)
        self.domega_w = np.ones_like(self.omega_w) * self.domega
        self.domega_w[0] *= 0.5
            
        prnt('Minimum eigenvalue: %10.3f eV' % (epsmin * Hartree),
             file=self.fd)
        prnt('Maximum eigenvalue: %10.3f eV' % (epsmax * Hartree),
             file=self.fd)
        prnt('Maximum frequency: %10.3f eV' % (self.omegamax * Hartree),
             file=self.fd)
        prnt('Number of frequencies:', len(self.omega_w), file=self.fd)
    
    def calculate(self):
        kd = self.calc.wfs.kd

        self.calculate_screened_potential()
        
        mykpts = [self.get_k_point(s, K, n1, n2)
                  for s, K, n1, n2 in self.mysKn1n2]

        for s in range(self.calc.wfs.nspins):
            for i, k1 in enumerate(self.kpts):
                K1 = kd.ibz2bz_k[k1]
                kpt1 = self.get_k_point(s, K1, *self.bands)
                self.eps_sin[s, i] = kpt1.eps_n
                for kpt2 in mykpts:
                    if kpt2.s == s:
                        self.calculate_q(i, kpt1, kpt2)
        prnt(np.array_str(self.eps_sin * Hartree, precision=3), file=self.fd)
        prnt(np.array_str(self.sigma_sin * Hartree, precision=3), file=self.fd)
        prnt(np.array_str(self.dsigma_sin, precision=3), file=self.fd)

    def calculate_q(self, i, kpt1, kpt2):
        wfs = self.calc.wfs
        qd = self.qd
        q_c = wfs.kd.bzk_kc[kpt2.K] - wfs.kd.bzk_kc[kpt1.K]
        #Q = abs((qd.bzk_kc - q_c) % 1).sum(axis=1).argmin()
        Q = abs(qd.bzk_kc - q_c).sum(axis=1).argmin()
        s = qd.sym_k[Q]
        #U_cc = qd.symmetry.op_scc[s]
        #time_reversal = qd.time_reversal_k[Q]
        iq = qd.bz2ibz_k[Q]
        iq_c = qd.ibzk_kc[iq]
        
        #sign = 1 - 2 * time_reversal
        shift_c = iq_c - q_c
        #shift_c = np.dot(U_cc, iq_c) - q_c * sign
        assert np.allclose(shift_c.round(), shift_c)
        shift_c = shift_c.round().astype(int)
        
        #if (U_cc == np.eye(3)).all():
        #    pass
        #else:
        #    N_c = self.calc.wfs.gd.N_c
        #    i_cr = np.dot(U_cc.T, np.indices(N_c).reshape((3, -1)))
        #    i = np.ravel_multi_index(i_cr, N_c, 'wrap')
        #    sdfg
            
        qd = KPointDescriptor([q_c])
        pd = PWDescriptor(self.ecut, wfs.gd, complex, qd)
        Q_G = self.get_fft_indices(kpt1.K, kpt2.K, q_c, pd,
                                   kpt1.shift_c - kpt2.shift_c)

        Q_aGii = self.initialize_paw_corrections(pd)
        
        for n in range(kpt1.n2 - kpt1.n1):
            fd = open('W.q%d.%s.npy' % (iq, self.filename))
            ut1cc_R = kpt1.ut_nR[n].conj()
            eps1 = kpt1.eps_n[n]
            C1_aGi = [np.dot(Q_Gii, P1_ni[n].conj())
                     for Q_Gii, P1_ni in zip(Q_aGii, kpt1.P_ani)]
            n_mG = self.calculate_pair_densities(ut1cc_R, C1_aGi, kpt2,
                                                 pd, Q_G)
            f_m = kpt2.f_n
            deps_m = eps1 - kpt2.eps_n
            sigma, dsigma = self.calculate_sigma(fd, n_mG, deps_m, f_m)
            self.sigma_sin[s, i, n] += sigma
            self.dsigma_sin[s, i, n] += dsigma

    def calculate_sigma(self, fd, n_mG, deps_m, f_m):
        sigma = 0.0
        dsigma = 0.0
        
        for omegap, domegap in zip(self.omega_w, self.domega_w):
            W_GG = np.load(fd)
            x1_m = 1 / (deps_m + omegap + 2j * self.eta * (f_m - 0.5))
            x2_m = 1 / (deps_m - omegap + 2j * self.eta * (f_m - 0.5))
            x_m = x1_m + x2_m
            dx_m = x1_m**2 + x2_m**2
            nW_mG = np.dot(n_mG, W_GG)
            sigma -= domegap * np.vdot(n_mG * x_m[:, np.newaxis], nW_mG).imag
            dsigma += domegap * np.vdot(n_mG * dx_m[:, np.newaxis], nW_mG).imag

        x = 1 / (self.qd.nbzkpts * 2 * pi * self.vol)
        return x * sigma, x * dsigma
        
    def calculate_screened_potential(self):
        chi0 = None
        for iq, q_c in enumerate(self.qd.ibzk_kc):
            fd = opencew('W.q%d.%s.npy' % (iq, self.filename))
            if fd is None:
                continue
                
            if chi0 is None:
                prnt('Calulating screened Coulomb potential:', file=self.fd)
                # Chi_0 calculator:
                chi0 = Chi0(self.calc,
                            self.omega_w * Hartree,
                            ecut=self.ecut * Hartree,
                            eta=self.eta * Hartree,
                            timeordered=True,
                            hilbert=not True,
                            real_space_derivatives=True)
                #wstc = WignerSeitzTruncatedCoulomb(self.calc.wfs.gd.cell_cv,
                #                                   self.calc.wfs.kd.N_c,
                #                                   self.fd)
            
            prnt(q_c, file=self.fd)
            pd, chi0_wGG = chi0.calculate(q_c)[:2]
            prnt(chi0_wGG.shape, file=self.fd)

            #iG_G = (wstc.get_potential(pd) / (4 * pi))**0.5
            iG_G = pd.G2_qG[0]**-0.5
            
            if not q_c.any():
                #chi0_wGG[:, 0] = 0.0
                #chi0_wGG[:, :, 0] = 0.0
                dq3 = (2 * pi)**3 / (self.qd.nbzkpts * self.vol)
                qc = (dq3 / 4 / pi * 3)**(1 / 3)
                G0inv = 2 * pi * qc**2 / dq3
                G20inv = 4 * pi * qc / dq3
                iG_G[0] = 1
                
            delta_GG = np.eye(len(iG_G))
            
            for chi0_GG in chi0_wGG:
                e_GG = delta_GG - 4 * pi * chi0_GG * iG_G * iG_G[:, np.newaxis]
                W_GG = 4 * pi * (np.linalg.inv(e_GG) -
                                 delta_GG) * iG_G * iG_G[:, np.newaxis]
                if not q_c.any():
                    W_GG[0, 0] *= G20inv
                    W_GG[1:, 0] *= G0inv
                    W_GG[0, 1:] *= G0inv
                    
                np.save(fd, W_GG)
            fd.close()
