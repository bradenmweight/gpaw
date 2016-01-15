from __future__ import print_function

import functools
import pickle
from time import time, ctime
from datetime import timedelta

import sys

import numpy as np
from ase.io import write
from ase.units import Hartree
from ase.utils import devnull
from ase.dft import monkhorst_pack

from gpaw import GPAW
from gpaw.kpt_descriptor import KPointDescriptor
from gpaw.wavefunctions.pw import PWDescriptor

from gpaw.io.tar import Writer, Reader
from gpaw.utilities.memory import maxrss
from gpaw.blacs import BlacsGrid, Redistributor
from gpaw.mpi import world, serial_comm
from gpaw.response.chi0 import Chi0
from gpaw.response.kernels import get_coulomb_kernel
from gpaw.response.kernels import get_integrated_kernel
from gpaw.response.wstc import WignerSeitzTruncatedCoulomb
from gpaw.response.pair import PairDensity


class BSE():
    def __init__(self, 
                 calc=None,
                 ecut=10.,
                 nbands=None,
                 valence_bands=None,
                 conduction_bands=None,
                 eshift=None,
                 gw_skn=None,
                 truncation=None,
                 txt=sys.stdout,
                 mode='BSE',
                 wfile='W_qGG.pckl',
                 write_h=True,
                 write_v=True): 
        """Creates the BSE object
        
        calc: str or calculator object
            The string should refer to the .gpw file contaning KS orbitals
        ecut: float
            Plane wave cutoff energy (eV)
        nbands: int
            Number of bands used for the screened interaction
        valence_bands: list
            Valence bands used in the BSE Hamiltonian
        conduction_bands: list
            Conduction bands used in the BSE Hamiltonian
        eshift: float
            Scissors operator opening the gap (eV)
        gw_skn: list / array
            List or array defining the gw quasiparticle energies used in the BSE 
            Hamiltonian. Should mathc spin, k-points and valence/conduction bands
        truncation: str
            Coulomb truncation scheme. Can be either wigner-seitz, 
            2D, 1D, or 0D
        txt: str
            txt output
        mode: str
            Theory level used. can be RPA TDHF or BSE. Only BSE is screened.
        wfile: str
            File for saving screened interaction and some other stuff needed later
        write_h: bool
            If True, write the BSE Hamiltonian to H_SS.gpw.
        write_v: bool
            If True, write eigenvalues and eigenstates to v_TS.gpw
        """

        assert mode in ['RPA', 'TDHF', 'BSE']
        assert calc.wfs.kd.nbzkpts % world.size == 0

        # Calculator
        if isinstance(calc, str):
            calc = GPAW(calc, txt=None, communicator=serial_comm)
        self.calc = calc

        # txt file
        if world.rank != 0:
            txt = devnull
        elif isinstance(txt, str):
            txt = open(txt, 'w', 1)
        self.fd = txt

        self.ecut = ecut / Hartree
        self.nbands = nbands
        self.mode = mode
        self.truncation = truncation
        self.wfile = wfile
        self.write_h = write_h
        self.write_v = write_v

        # Find q-vectors and weights in the IBZ:
        self.kd = calc.wfs.kd
        assert -1 not in self.kd.bz2bz_ks
        offset_c = 0.5 * ((self.kd.N_c + 1) % 2) / self.kd.N_c
        bzq_qc = monkhorst_pack(self.kd.N_c) + offset_c
        self.qd = KPointDescriptor(bzq_qc)
        self.qd.set_symmetry(self.calc.atoms, self.kd.symmetry)
        self.vol = abs(np.linalg.det(calc.wfs.gd.cell_cv))

        # bands
        if valence_bands is None:
            valence_bands = [self.calc.wfs.setups.nvalence // 2]
        if conduction_bands is None:
            conduction_bands = [valence_bands[-1] + 1]
        self.val_n = valence_bands
        self.con_n = conduction_bands
        self.td = True
        for n in self.val_n:
            if n in self.con_n:
                self.td = False
        self.nv = len(valence_bands)
        self.nc = len(conduction_bands)
        if eshift is not None:
            eshift /= Hartree
        if gw_skn is not None:
            assert self.nv + self.nc == len(gw_skn[0, 0])
            assert self.kd.nbzkpts == len(gw_skn[0])
            gw_skn /= Hartree
        self.gw_skn = gw_skn
        self.eshift = eshift
        self.nS = self.kd.nbzkpts * self.nv * self.nc
        
        self.print_initialization(self.td, self.eshift, self.gw_skn)

    def calculate(self, optical=True, q_c=[0.0, 0.0, 0.0], ac=1.0):

        # Parallelization stuff
        nK = self.kd.nbzkpts
        myKsize = -(-nK // world.size)
        myKrange = range(world.rank * myKsize,
                         min((world.rank + 1) * myKsize, nK))
        myKsize = len(myKrange)

        # Calculate direct (screened) interaction
        self.Q_qaGii = []
        self.W_qGG = []
        self.pd_q = []
        self.get_screened_potential(ac=ac)

        # Calculate exchange interaction
        self.Pair = PairDensity(self.calc, self.ecut, world=serial_comm,
                                txt='pair.txt')
        iq0 = self.qd.bz2ibz_k[self.kd.where_is_q(self.q_c, self.qd.bzk_kc)]
        pd0 = self.pd_q[iq0]
        v_G = get_coulomb_kernel(pd0, self.kd.N_c, truncation=self.truncation,
                                 wstc=self.wstc)
        if optical:
            v_G[0] = 0.0

        # Calculate pair densities, eigenvalues and occupations
        rhoex_KmnG = np.zeros((nK, self.nv, self.nc, len(v_G)), complex)
        rhoG0_Kmn = np.zeros((nK, self.nv, self.nc), complex)
        df_Kmn = np.zeros((nK, self.nv, self.nc), float)
        deps_kmn = np.zeros((myKsize, self.nv, self.nc), float)
        if np.allclose(self.q_c, 0.0):
            optical_limit = True
        else:
            optical_limit = False
        get_pair = self.Pair.get_kpoint_pair
        get_rho = self.Pair.get_pair_density
        vi, vf = self.val_n[0], self.val_n[-1] + 1
        ci, cf = self.con_n[0], self.con_n[-1] + 1
        for ik, iK in enumerate(myKrange):
            pair = get_pair(self.pd_q[iq0], 0, iK, vi, vf, ci, cf)
            #deps_nm = (self.gw_skn[0, ik, :self.nv][:, np.newaxis] -
            #           self.gw_skn[0, ik, self.nv:])

            deps_kmn[ik] = -pair.get_transition_energies(range(self.nv),
                                                         range(self.nc))
            df_Kmn[iK] = pair.get_occupation_differences(range(self.nv),
                                                         range(self.nc))
            rhoex_KmnG[iK] = get_rho(self.pd_q[iq0], pair,
                                     range(self.nv), range(self.nc),
                                     optical_limit=optical_limit,
                                     Q_aGii=self.Q_qaGii[iq0])[0]
        if self.eshift is not None:
            deps_kmn[np.where(df_Kmn[myKrange] > 1.0e-3)] += self.eshift
            deps_kmn[np.where(df_Kmn[myKrange] < -1.0e-3)] -= self.eshift
        df_Kmn *= 2.0 / nK  # multiply by 2 for spin polarized calculation
        world.sum(rhoex_KmnG)
        world.sum(df_Kmn)

        # Calculate Hamiltonian
        t0 = time()
        print('Calculating %s matrix elements' % self.mode, file=self.fd)
        H_kKmnmn = np.zeros((myKsize, nK,
                             self.nv, self.nc, self.nv, self.nc),
                            complex)
        for ik1, iK1 in enumerate(myKrange):
            rho1_mnG = rhoex_KmnG[iK1]
            rho1ccV_mnG = rho1_mnG.conj()[:, :] * v_G
            rhoG0_Kmn[iK1] = rho1_mnG[:, :, 0]
            kptv1 = self.Pair.get_k_point(0, iK1, vi, vf)
            kptc1 = self.Pair.get_k_point(0, iK1, ci, cf)
            for iK2 in range(nK):
                rho2_mnG = rhoex_KmnG[iK2]
                H_kKmnmn[ik1, iK2] += np.dot(rho1ccV_mnG,
                                             np.swapaxes(rho2_mnG, 1, 2))
                if not self.mode == 'RPA':
                    kptv2 = self.Pair.get_k_point(0, iK2, vi, vf)
                    kptc2 = self.Pair.get_k_point(0, iK2, ci, cf)
                    rho3_mmG, iq = self.get_density_matrix(kptv1, kptv2)
                    rho4_nnG, iq = self.get_density_matrix(kptc1, kptc2)
                    rho3ccW_mmG = np.dot(rho3_mmG.conj(), self.W_qGG[iq])
                    W_mmnn = np.dot(rho3ccW_mmG,
                                    np.swapaxes(rho4_nnG, 1, 2))
                    H_kKmnmn[ik1, iK2] -= 0.5 * np.swapaxes(W_mmnn, 1, 2)
                    
            if iK1 % (myKsize // 5 + 1) == 0:
                dt = time() - t0
                tleft = dt * myKsize / (iK1 + 1) - dt
                print('  Finished %s pair orbitals in %s - Estimated %s left' % 
                      ((iK1 + 1) * self.nv * self.nc * world.size,
                       timedelta(seconds=round(dt)),
                       timedelta(seconds=round(tleft))), file=self.fd)
        H_kKmnmn /= self.vol

        # From here on s is a local pair-orbital index
        mySsize = myKsize * self.nv * self.nc
        if myKsize > 0:
            iS0 = myKrange[0] * self.nv * self.nc
            
        # Reshape and collect
        world.sum(rhoG0_Kmn)
        self.rhoG0_S = np.reshape(rhoG0_Kmn, -1)
        self.df_S = np.reshape(df_Kmn, -1)
        self.deps_s = np.reshape(deps_kmn, -1)
        H_sS = np.reshape(np.swapaxes(np.swapaxes(H_kKmnmn, 1, 2), 2, 3),
                          (mySsize, self.nS))

        for iS in range(mySsize):
            # Multiply by occupations and adiabatic coupling
            H_sS[iS] *= self.df_S[iS0 + iS] * ac
            # add bare transition energies
            H_sS[iS, iS0 + iS] += self.deps_s[iS]

        # Save H_sS matrix
        if self.write_h:
            self.par_save('H_SS','H_SS', H_sS)

        self.H_sS = H_sS

    def get_density_matrix(self, kpt1, kpt2):

        Q_c = self.kd.bzk_kc[kpt2.K] - self.kd.bzk_kc[kpt1.K]
        iQ = self.qd.where_is_q(Q_c, self.qd.bzk_kc)
        iq = self.qd.bz2ibz_k[iQ]
        q_c = self.qd.ibzk_kc[iq]
        
        if np.allclose(q_c, 0.0):
            optical_limit = True
        else:
            optical_limit = False

        # Find symmetry that transforms Q_c into q_c
        sym = self.qd.sym_k[iQ]
        U_cc = self.qd.symmetry.op_scc[sym]
        time_reversal = self.qd.time_reversal_k[iQ]
        sign = 1 - 2 * time_reversal
        d_c = sign * np.dot(U_cc, q_c) - Q_c
        assert np.allclose(d_c.round(), d_c)
        
        pd = self.pd_q[iq]
        N_c = pd.gd.N_c
        i_cG = sign * np.dot(U_cc, np.unravel_index(pd.Q_qG[0], N_c))

        shift0_c = Q_c - sign * np.dot(U_cc, pd.kd.bzk_kc[0])
        assert np.allclose(shift0_c.round(), shift0_c)
        shift0_c = shift0_c.round().astype(int)

        shift_c = kpt1.shift_c - kpt2.shift_c - shift0_c
        I_G = np.ravel_multi_index(i_cG + shift_c[:, None], N_c, 'wrap')
        
        G_Gv = pd.get_reciprocal_vectors()
        pos_ac = self.calc.atoms.get_scaled_positions()
        pos_av = np.dot(pos_ac, pd.gd.cell_cv)
        M_vv = np.dot(pd.gd.cell_cv.T, np.dot(U_cc.T,
                                              np.linalg.inv(pd.gd.cell_cv).T))

        Q_aGii = []
        for a, Q_Gii in enumerate(self.Q_qaGii[iq]):
            x_G = np.exp(1j * np.dot(G_Gv, (pos_av[a] -
                                            sign * np.dot(M_vv, pos_av[a]))))
            U_ii = self.calc.wfs.setups[a].R_sii[sym]
            Q_Gii = np.dot(np.dot(U_ii, Q_Gii * x_G[:, None, None]),
                           U_ii.T).transpose(1, 0, 2)
            Q_aGii.append(Q_Gii)
        
        rho_mnG = np.zeros((len(kpt1.eps_n), len(kpt2.eps_n), len(G_Gv)),
                           complex)
        for m in range(len(rho_mnG)):
            C1_aGi = [np.dot(Qa_Gii, P1_ni[m].conj())
                      for Qa_Gii, P1_ni in zip(Q_aGii, kpt1.P_ani)]
            ut1cc_R = kpt1.ut_nR[m].conj()
            rho_mnG[m] = self.Pair.calculate_pair_densities(ut1cc_R, C1_aGi, 
                                                            kpt2, pd, I_G)
        return rho_mnG, iq

    def get_screened_potential(self, ac=1.0):

        if self.truncation == 'wigner-seitz':
            self.wstc = WignerSeitzTruncatedCoulomb(self.calc.wfs.gd.cell_cv, 
                                                    self.kd.N_c, self.fd)
        else:
            self.wstc = None

        if self.wfile is not None:
            # Read screened potential from file
            try:
                f = open(self.wfile)
                print('Reading screened potential from % s' % self.wfile, 
                      file=self.fd)
                self.Q_qaGii, self.pd_q, self.W_qGG = pickle.load(f)
            except:
                self.calculate_screened_potential(ac)
                print('Saving screened potential to % s' % self.wfile,
                      file=self.fd)
                f = open(self.wfile, 'w')
                pickle.dump((self.Q_qaGii, self.pd_q, self.W_qGG),
                            f, pickle.HIGHEST_PROTOCOL)
        else:
            self.calculate_screened_potential(ac)

    def calculate_screened_potential(self, ac):
        """Calculate W_GG(q)"""

        chi0 = Chi0(self.calc,
                    frequencies=[0.0],
                    eta=0.001,
                    ecut=self.ecut,
                    intraband=False,
                    hilbert=False,
                    nbands=self.nbands,
                    txt='chi0.txt',
                    world=world,
                    )

        self.blockcomm = chi0.blockcomm
        wfs = self.calc.wfs

        t0 = time()
        print('Calculating screened potential', file=self.fd)
        for iq, q_c in enumerate(self.qd.ibzk_kc):
            thisqd = KPointDescriptor([q_c])
            pd = PWDescriptor(self.ecut, wfs.gd, complex, thisqd)
            nG = pd.ngmax

            chi0.Ga = self.blockcomm.rank * nG
            chi0.Gb = min(chi0.Ga + nG, nG)
            chi0_wGG = np.zeros((1, nG, nG), complex)
            if np.allclose(q_c, 0.0):
                chi0_wxvG = np.zeros((1, 2, 3, nG), complex)
                chi0_wvv = np.zeros((1, 3, 3), complex)
            else:
                chi0_wxvG = None
                chi0_wvv = None

            chi0._calculate(pd, chi0_wGG, chi0_wxvG, chi0_wvv,
                            0, self.nbands, [0, 1])
            chi0_GG = chi0_wGG[0]

            # Calculate eps^{-1}_GG
            if pd.kd.gamma:
                # Generate fine grid in vicinity of gamma
                kd = self.calc.wfs.kd
                N = 4
                N_c = [N, N, N]
                if not np.all(kd.N_c == [1, 1, 1]):
                    N_c[np.where(kd.N_c == 1)[0]] = 1
                qf_qc = monkhorst_pack(N_c)
                qf_qc *= 1.0e-6
                U_scc = kd.symmetry.op_scc
                qf_qc = kd.get_ibz_q_points(qf_qc, U_scc)[0]
                weight_q = kd.q_weights
                qf_qv = 2 * np.pi * np.dot(qf_qc, pd.gd.icell_cv)
                a_q = np.sum(np.dot(chi0_wvv[0], qf_qv.T)
                             * qf_qv.T, axis=0)
                a0_qG = np.dot(qf_qv, chi0_wxvG[0, 0])
                a1_qG = np.dot(qf_qv, chi0_wxvG[0, 1])
                einv_GG = np.zeros((nG, nG), complex)
                #W_GG = np.zeros((nG, nG), complex)
                for iqf in range(len(qf_qv)):
                    chi0_GG[0] = a0_qG[iqf]
                    chi0_GG[:, 0] = a1_qG[iqf]
                    chi0_GG[0, 0] = a_q[iqf]
                    sqrv_G = get_coulomb_kernel(pd,
                                                kd.N_c,
                                                truncation=self.truncation,
                                                wstc=self.wstc,
                                                q_v=qf_qv[iqf])**0.5
                    sqrv_G *= ac**0.5 # Multiply by adiabatic coupling
                    e_GG = np.eye(nG) - chi0_GG * sqrv_G * sqrv_G[:, np.newaxis]
                    einv_GG += np.linalg.inv(e_GG) * weight_q[iqf]
                    #einv_GG = np.linalg.inv(e_GG) * weight_q[iqf]
                    #W_GG += (einv_GG * sqrv_G * sqrv_G[:, np.newaxis] 
                    #         * weight_q[iqf])
            else:
                sqrv_G = get_coulomb_kernel(pd,
                                            self.kd.N_c,
                                            truncation=self.truncation,
                                            wstc=self.wstc)**0.5
                sqrv_G *= ac**0.5 # Multiply by adiabatic coupling
                e_GG = np.eye(nG) - chi0_GG * sqrv_G * sqrv_G[:, np.newaxis]
                einv_GG = np.linalg.inv(e_GG)
                #W_GG = einv_GG * sqrv_G * sqrv_G[:, np.newaxis]

            # Now calculate W_GG
            if pd.kd.gamma:
                sqrv_G = get_coulomb_kernel(pd,
                                            self.kd.N_c,
                                            truncation=self.truncation,
                                            wstc=self.wstc)**0.5 
                #bzvol = (2 * np.pi)**3 / self.vol / self.qd.nbzkpts
                #Rq0 = (3 * bzvol / (4 * np.pi))**(1. / 3.)
                #sqrv_G[0] = 4 * np.pi * (Rq0 / bzvol)**0.5
            sqrv_G[0] = get_integrated_kernel(pd,
                                              self.kd.N_c,
                                              truncation=self.truncation,
                                              N=100)**0.5 
            W_GG = einv_GG * sqrv_G * sqrv_G[:, np.newaxis]

            if pd.kd.gamma:
                e = 1 / einv_GG[0, 0].real
                print('RPA macroscopic dielectric constant is: %3.3f' %  e,
                      file=self.fd)
            self.Q_qaGii.append(chi0.Q_aGii)
            self.pd_q.append(pd)
            self.W_qGG.append(W_GG)

            if iq % (self.qd.nibzkpts // 5 + 1) == 0:
                dt = time() - t0
                tleft = dt * self.qd.nibzkpts / (iq + 1) - dt
                print('  Finished %s q-points in %s - Estimated %s left' % 
                      (iq, timedelta(seconds=round(dt)),
                       timedelta(seconds=round(tleft))), file=self.fd)

    def diagonalize(self):

        print('Diagonalizing Hamiltonian', file=self.fd)
        """The t and T represent local and global 
           eigenstates indices respectively
        """
 
        # Non-Hermitian matrix can only use linalg.eig
        if not self.td: 
            print('  Using numpy.linalg.eig...', file=self.fd)
            self.H_SS = np.zeros((self.nS, self.nS), dtype=complex)
            world.all_gather(self.H_sS, self.H_SS)
            self.w_T, self.v_St = np.linalg.eig(self.H_SS)
        # Here the eigenvectors are returned as complex conjugated rows
        else:
            if world.size == 1:
                print('  Using lapack...', file=self.fd)
                from gpaw.utilities.lapack import diagonalize
                self.w_T = np.zeros(self.nS)
                diagonalize(self.H_sS, self.w_T)
                self.v_St = self.H_sS.conj().T
            else:
                print('  Using scalapack...', file=self.fd)
                self.w_T, self.v_St = self.scalapack_diagonalize(self.H_sS)
                self.v_St = self.v_St.conj().T

        if self.write_v:
            self.par_save('v_TS', 'v_TS', self.v_St.T)
        return 

    def scalapack_diagonalize(self, H_sS):

        mb = 32
        N = self.nS
        
        g1 = BlacsGrid(world, world.size,    1)
        g2 = BlacsGrid(world, world.size // 2, 2)
        nndesc1 = g1.new_descriptor(N, N, len(H_sS),  N) 
        nndesc2 = g2.new_descriptor(N, N, mb, mb)
        
        A_ss = nndesc2.empty(dtype=H_sS.dtype)
        redistributor = Redistributor(world, nndesc1, nndesc2)
        redistributor.redistribute(H_sS, A_ss)
        
        # diagonalize
        v_ss = nndesc2.zeros(dtype=A_ss.dtype)
        w_S = np.zeros(N, dtype=float)
        nndesc2.diagonalize_dc(A_ss, v_ss, w_S, 'L')
        
        # distribute the eigenvectors to master
        v_sS = np.zeros_like(H_sS)
        redistributor = Redistributor(world, nndesc2, nndesc1)
        redistributor.redistribute(v_ss, v_sS)
        
        return w_S, v_sS

    def get_vchi(self, w_w=None, eta=0.1, q_c=[0.0, 0.0, 0.0],
                 ac=1.0, readfile=None):
        """Returns v * \chi where v is the bare Coulomb interaction"""

        self.q_c = q_c

        if readfile is None:
            self.calculate(optical=True, q_c=q_c, ac=ac)
            self.diagonalize()
        elif readfile == 'H_SS':
            print('Reading Hamiltonian from file', file=self.fd)
            self.par_load('H_SS', 'H_SS')
            self.diagonalize()
        elif readfile == 'v_TS':
            print('Reading eigenstates from file', file=self.fd)
            self.par_load('v_TS', 'v_TS') 
        else:
            raise ValueError('%s array not recognized' % readfile)

        w_w /= Hartree
        eta /= Hartree

        w_T= self.w_T
        v_St = self.v_St
        rhoG0_S = self.rhoG0_S
        df_S = self.df_S

        print('Calculating response function at %s frequency points' %
              len(w_w), file=self.fd)
        vchi_w = np.zeros(len(w_w), dtype=complex)

        A_t = np.dot(rhoG0_S, v_St)
        B_t = np.dot(rhoG0_S * df_S, v_St)
        if not self.td:
            # Indices are global in this case (t=T, s=S)
            tmp = np.dot(v_St.conj().T, v_St )
            overlap_tt = np.linalg.inv(tmp)
            C_T = np.dot(B_t.conj(), overlap_tt.T) * A_t
        else:
            C_T = np.zeros(self.nS, complex)
            world.all_gather(B_t.conj() * A_t, C_T)

        for iw, w in enumerate(w_w):
            tmp_T = 1. / (w - w_T + 1j * eta)
            vchi_w[iw] += np.dot(tmp_T, C_T)
        vchi_w *=  4 * np.pi / self.vol

        if not np.allclose(self.q_c, 0.0):
            cell_cv = self.calc.wfs.gd.cell_cv
            B_cv = 2 * np.pi * np.linalg.inv(cell_cv).T
            q_v = np.dot(q_c, B_cv)
            vchi_w /= np.dot(q_v, q_v)

        return vchi_w * ac

    def get_dielectric_function(self, w_w=None, eta=0.1, q_c=[0.0, 0.0, 0.0],
                                filename='df_bse.csv', readfile=None,
                                write_eig='bse_eig.dat'):
        """Returns and writes real and imaginary part of the dielectric function.

        w_w: list of frequencies (eV)
            Dielectric function is calculated at these frequencies
        eta: float
            Lorentzian broadening of the spectrum (eV)
        q_c: list of three floats
            Wavevector in reduced units on which the response is calculated 
        filename: str
            data file on which frequencies, real and imaginary part of 
            dielectric function is written
        readfile: str
            If H_SS is given, the method will load the BSE Hamiltonian from H_SS.gpw
            If v_TS is given, the method will load the eigenstates from v_TS.gpw
        write_eig: str
            File on which the BSE eigenvalues are written
        """

        epsilon_w = -self.get_vchi(w_w=w_w, eta=eta, q_c=q_c, readfile=readfile)
        epsilon_w += 1.0
    
        """Check f-sum rule."""
        nv = self.calc.wfs.setups.nvalence
        dw_w = w_w[1:] - w_w[:-1]
        weps_w = (w_w[1:] + w_w[:-1]) * (epsilon_w[1:] + epsilon_w[:-1]) / 4
        N = np.dot(dw_w, weps_w.imag) * self.vol / (2 * np.pi**2)
        print(file=self.fd)
        print('Checking f-sum rule:', file=self.fd)
        print('  N = %f, %f  %% error' % (N, (N - nv) / nv * 100), 
              file=self.fd)
        print(file=self.fd)

        w_w *= Hartree
        if world.rank == 0:
            f = open(filename, 'w')
            for iw, w in enumerate(w_w):
                print('%.6f, %.6f, %.6f' %  
                      (w, epsilon_w[iw].real, epsilon_w[iw].imag), file=f)
            f.close()
        world.barrier()

        self.w_T *= Hartree
        if write_eig is not None:
            if world.rank == 0:
                f = open(write_eig, 'w')
                print('# %s eigenvalues in eV' % self.mode, file=f)
                for iw, w in enumerate(self.w_T):
                    print('%8d %12.6f' % (iw, w.real), file=f)
                f.close()
            
        print('Calculation completed at:', ctime(), file=self.fd)

        return w_w, epsilon_w

    def par_save(self, filename, name, A_sS):
        from gpaw.io import open 

        nS = self.nS
        
        if world.rank == 0:
            w = open(filename, 'w', world)
            w.dimension('nS', nS)
            
            if name == 'v_TS':
                w.add('w_T', ('nS',), dtype=self.w_T.dtype)
                w.fill(self.w_T)
            w.add('rhoG0_S', ('nS',), dtype=complex)
            w.fill(self.rhoG0_S)
            w.add('df_S', ('nS',), dtype=complex)
            w.fill(self.df_S)

            w.add(name, ('nS', 'nS'), dtype=complex)
            tmp = np.zeros_like(A_sS)

        # Assumes that H_SS is written in order from rank 0 - rank N
        for irank in range(world.size):
            if irank == 0:
                if world.rank == 0:
                    w.fill(A_sS)
            else:
                if world.rank == irank:
                    world.send(A_sS, 0, irank+100)
                if world.rank == 0:
                    world.receive(tmp, irank, irank+100)
                    w.fill(tmp)
        if world.rank == 0:
            w.close()
        world.barrier()

    def par_load(self, filename, name):
        from gpaw.io import open 

        r = open(filename, 'r')
        nS = r.dimension('nS')
        mySsize = -(-nS // world.size)
        mySrange = range(world.rank * mySsize,
                         min((world.rank + 1) * mySsize, nS))
        mySsize = len(mySrange)

        if name == 'H_SS':
            self.H_sS = np.zeros((mySsize, nS), dtype=complex)
            for si, s in enumerate(mySrange):
                self.H_sS[si] = r.get('H_SS', s)

        if name == 'v_TS':
            self.w_T = r.get('w_T')
            v_tS = np.zeros((mySsize, nS), dtype=complex)
            for it, t in enumerate(mySrange):
                v_tS[it] = r.get('v_TS', t)
            self.v_St = v_tS.T

        self.rhoG0_S = r.get('rhoG0_S')
        self.df_S = r.get('df_S')

        r.close()
        
    def print_initialization(self, td, eshift, gw_skn):
        p = functools.partial(print, file=self.fd)
        p('----------------------------------------------------------')
        p('BSE Calculation')
        p('----------------------------------------------------------')
        p('Started at:  ', ctime())
        p()
        p('Atoms                          :',
          self.calc.atoms.get_chemical_formula(mode='hill'))
        p('Ground state XC functional     :', self.calc.hamiltonian.xc.name)
        p('Valence electrons              :', self.calc.wfs.setups.nvalence)
        p('Number of bands                :', self.calc.wfs.bd.nbands)
        p('Number of spins                :', self.calc.wfs.nspins)
        p('Number of k-points             :', self.kd.nbzkpts)
        p('Number of irreducible k-points :', self.kd.nibzkpts)
        p('Number of q-points             :', self.qd.nbzkpts)
        p('Number of irreducible q-points :', self.qd.nibzkpts)
        p()
        for q in self.qd.ibzk_kc:
            p('    q: [%1.4f %1.4f %1.4f]' % (q[0], q[1], q[2]))
        p()
        if gw_skn is not None:
            p('User specified BSE bands')
        p('Screening bands included       :', self.nbands)
        p('BSE valence bands              :', self.val_n)
        p('BSE conduction bands           :', self.con_n)
        if eshift is not None:
            p('Scissors operator              :', eshift * Hartree, 'eV')
        p('Tamm-Dancoff approximation     :', td)
        p('Number of pair orbitals        :', self.nS)
        p()
        p('----------------------------------------------------------')
        p('----------------------------------------------------------')
        p()
        p('Parallelization - Total number of CPUs   : % s' % world.size)
        p('  Screened potential')
        p('    K-point/band decomposition           : % s' % world.size)
        p('  Hamiltonian')
        p('    Pair orbital decomposition           : % s' % world.size)
        p()


