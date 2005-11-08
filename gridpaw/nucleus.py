# Copyright (C) 2003  CAMP
# Please see the accompanying LICENSE file for further information.

"""Nucleus class.

A Paw object has a list of nuclei. Each nucleus is described by a
``Setup`` object and a scaled position plus some extra stuff..."""

import Numeric as num

from gridpaw.utilities.complex import real, cc
from gridpaw.localized_functions import create_localized_functions
from gridpaw.utilities import unpack, pack
import gridpaw.utilities.mpi as mpi


NOT_INITIALIZED = -1
NOTHING = 0
COMPENSATION_CHARGE = 1
PROJECTOR_FUNCTION = 2
EVERYTHING = 3


class Nucleus:
    """Nucleus-class.

    The ``Nucleus`` object basically consists of a ``Setup`` object, a
    scaled position and some localized functions.  It takes care of
    adding localized functions to functions on extended grids and
    calculating integrals of functions on extended grids and localized
    functions.

     ============= ========================================================
     ``setup``     ``Setup`` object.
     ``spos_i``    Scaled position.
     ``a``         Index number for this nucleus.
     ``typecode``  Data type of wave functions (``Float`` or ``Complex``).
     ``neighbors`` List of overlapping neighbor nuclei.
     ``onohirose`` Number of grid points used for Ono-Hirose_ interpolation
                   - currently 1, 2, 3, and 5 is implemented (1 means no
                   interpolation).
     ============= ========================================================

    Localized functions:
     ========== ===========================================================
     ``nct``    Pseudo core electron density.
     ``ghat_L`` Shape functions for compensation charges.
     ``vhat_L`` Correction potentials for overlapping compensation charges.
     ``pt_i``   Projector functions.
     ``vbar``   Arbitrary localized potential.
     ``phit_i`` Pseudo partial waves used for initial wave function guess.
     ========== ===========================================================

    Arrays:
     ========= ===============================================================
     ``P_uni`` Integral of products of all wave functions and the projector
               functions of this atom (``P_{\sigma\vec{k}ni}^a``).
     ``D_sp``  Atomic density matrix (``D_{\sigma i_1i_2}^a``).
     ``dH_sp`` Atomic Hamiltonian correction (``\Delta H_{\sigma i_1i_2}^a``).
     ``Q_L``   Multipole moments  (``Q_{\ell m}^a``).
     ``F_i``   Force.
     ========= ===============================================================

    Parallel stuff: ``comm``, ``rank`` and ``domain_overlap``
    """
    def __init__(self, setup, a, typecode, onohirose=1):
        """Construct a ``Nucleus`` object."""
        self.setup = setup
        self.a = a
        self.typecode = typecode
        self.onohirose = onohirose
        self.spos_c = None
        self.rank = None
        lmax = self.setup.lmax
        self.domain_overlap = NOT_INITIALIZED
        self.Q_L = num.zeros((lmax + 1)**2, num.Float)
        self.neighbors = []
        self.comm = mpi.serial_comm
        
    def allocate(self, nspins, nkpts, nbands):
        ni = self.get_number_of_partial_waves()
        np = ni * (ni + 1) / 2
        self.D_sp = num.zeros((nspins, np), num.Float)
        self.H_sp = num.zeros((nspins, np), num.Float)
        self.P_uni = num.zeros((nspins * nkpts, nbands, ni), self.typecode)
        self.F_c = num.zeros(3, num.Float)

    def reallocate(self, nbands):
        nu, nao, ni = self.P_uni.shape
        if nbands < nao:
            self.P_uni = self.P_uni[:, :nbands, :].copy()
        else:
            P_uni = num.zeros((nu, nbands, ni), self.typecode)
            P_uni[:, :nao, :] = self.P_uni
            self.P_uni = P_uni

    def set_communicators(self, comm, root, g_comm, g_root, p_comm, p_root):
        self.comm = comm

        if g_comm is not None:
            self.ghat_L.set_communicator(g_comm, g_root)
            self.vhat_L.set_communicator(g_comm, g_root)

        if p_comm is not None:
            self.pt_i.set_communicator(p_comm, p_root)
            if self.vbar is not None:
                self.vbar.set_communicator(p_comm, p_root)

        if self.nct is not None:
            self.nct.set_communicator(comm, root)

    def make_localized_grids(self, gd, finegd, k_ki, lfbc):
        # Smooth core density:
        nct = self.setup.get_smooth_core_density()
        create = create_localized_functions
        self.nct = create([nct], gd, self.spos_c, cut=True, lfbc=lfbc,
                          onohirose=self.onohirose)

        # Shape function:
        ghat_l = self.setup.get_shape_function()
        self.ghat_L = create(ghat_l, finegd, self.spos_c, lfbc=lfbc,
                             onohirose=self.onohirose)
            
        # Potential:
        vhat_l = self.setup.get_potential()
        self.vhat_L = create(vhat_l, finegd, self.spos_c, lfbc=lfbc,
                             onohirose=self.onohirose)

        assert (self.ghat_L is None) == (self.vhat_L is None)
        
        # Projectors:
        pt_j = self.setup.get_projectors()
        self.pt_i = create(pt_j, gd, self.spos_c, typecode=self.typecode,
                           lfbc=lfbc, onohirose=self.onohirose)
        if self.typecode == num.Complex:
            self.pt_i.set_phase_factors(k_ki)
        
        # Localized potential:
        vbar = self.setup.get_localized_potential()
        self.vbar = create([vbar], finegd, self.spos_c, lfbc=lfbc,
                           onohirose=self.onohirose)

        if self.pt_i is None:
            assert self.vbar is None
        
    def initialize_atomic_orbitals(self, gd, k_ki, lfbc):
        phit_j = self.setup.get_atomic_orbitals()
        self.phit_i = create_localized_functions(
            phit_j, gd, self.spos_c, onohirose=1, typecode=self.typecode,
            cut=True, forces=False, lfbc=lfbc)
        if self.typecode == num.Complex:
            self.phit_i.set_phase_factors(k_ki)

    def get_number_of_atomic_orbitals(self):
        return self.setup.get_number_of_atomic_orbitals()

    def get_number_of_partial_waves(self):
        return self.setup.get_number_of_partial_waves()
    
    def create_atomic_orbitals(self, psit_iG, k):
        if self.phit_i is None:
            # Nothing to do in this domain:
            return

        coefs = num.identity(len(psit_iG), psit_iG.typecode())
        self.phit_i.add(psit_iG, coefs, k)

    def add_atomic_density(self, nt_sG, magmom, hund):
        if self.phit_i is None:
            # Nothing to do in this domain:
            print 'self.phitj is none'
            return

        ns = len(nt_sG)
        niao = self.get_number_of_atomic_orbitals()
        f_si = num.zeros((ns, niao), num.Float)
        
        i = 0
        for l, f in zip(self.setup.l_j, self.setup.f_j):
            f = int(f)
            if f == 0:
                continue
            degeneracy = 2 * l + 1
            if hund:
                # Use Hunds rules:
                f_si[0, i:i + min(f, degeneracy)] = 1.0      # spin up
                f_si[1, i:i + max(f - degeneracy, 0)] = 1.0  # spin down
                if f < degeneracy:
                    magmom -= f
                else:
                    magmom -= 2 * degeneracy - f
            else:
                if ns == 1:
                    f_si[0, i:i + degeneracy] = 1.0 * f / degeneracy
                else:
                    if f == 2 * degeneracy:
                        mag = 0.0
                    else:
                        mag = min(magmom, degeneracy)
                    f_si[0, i:i + degeneracy] = 0.5 * (f + mag) / degeneracy
                    f_si[1, i:i + degeneracy] = 0.5 * (f - mag) / degeneracy
                    magmom -= mag
            i += degeneracy
        assert i == niao
        assert magmom == 0.0

        for s in range(ns):
            self.phit_i.add_density(nt_sG[s], f_si[s])

        if self.domain_overlap == EVERYTHING:
            ni = self.get_number_of_partial_waves()
            p = 0
            i0 = 0
            i = 0
            for f, l in zip(self.setup.f_j, self.setup.l_j):
                deg = 2 * l + 1
                for m in range(deg):
                    if f > 0.0: # XXX move out of loop?
                        self.D_sp[:, p] = f_si[:, i]
                        i += 1
                    p += ni - i0
                    i0 += 1
            assert i == niao
            assert i0 == ni
            assert p == ni * (ni + 1) / 2

    def add_smooth_core_density(self, nct_G):
        if self.nct is not None:
            self.nct.add(nct_G, num.ones(1, num.Float))
            
    def add_compensation_charge(self, nt2):
        self.ghat_L.add(nt2, self.Q_L)

    def add_hat_potential(self, vt2):
        self.vhat_L.add(vt2, self.Q_L)

    def add_localized_potential(self, vt2):
        if self.vbar is not None:
            self.vbar.add(vt2, num.array([1.0]))
        
    def calculate_projections(self, kpt):
        assert self.domain_overlap >= PROJECTOR_FUNCTION
        if self.domain_overlap == EVERYTHING:
            P_ni = self.P_uni[kpt.u]
            P_ni[:] = 0.0 # why????
            self.pt_i.integrate(kpt.psit_nG, P_ni, kpt.k)
        else:
            self.pt_i.integrate(kpt.psit_nG, None, kpt.k)

    def calculate_multipole_moments(self):
        if self.domain_overlap == EVERYTHING:
            self.Q_L[:] = num.dot(num.sum(self.D_sp), self.setup.Delta_pL)
            self.Q_L[0] += self.setup.Delta0
        self.comm.broadcast(self.Q_L, self.rank)
            
    def calculate_hamiltonian(self, nt, vHt_g):
        assert self.domain_overlap >= COMPENSATION_CHARGE
        if self.domain_overlap == EVERYTHING:
            a = self.setup
            W_L = num.zeros((a.lmax + 1)**2, num.Float)
            for neighbor in self.neighbors:
                W_L += num.dot(neighbor.v, neighbor.nucleus.Q_L)
            U = 0.5 * num.dot(self.Q_L, W_L)
        
            self.vhat_L.integrate(nt, W_L)
            self.ghat_L.integrate(vHt_g, W_L)

            Exc = a.xc.calculate_energy_and_derivatives(self.D_sp, self.H_sp)

            D_p = num.sum(self.D_sp)
            dH_p = (a.K_p + a.M_p + a.MB_p + 2.0 * num.dot(a.M_pp, D_p) +
                    num.dot(a.Delta_pL, W_L))
            Ekin = num.dot(a.K_p, D_p) + a.Kc

            Ebar = a.MB + num.dot(a.MB_p, D_p)
            Epot = U + a.M + num.dot(D_p, (a.M_p + num.dot(a.M_pp, D_p)))
            for H_p in self.H_sp:
                H_p += dH_p

            # Move this kinetic energy contribution to Paw.py: ????!!!!
            Ekin -= num.dot(self.D_sp[0], self.H_sp[0])
            if len(self.D_sp) == 2:
                Ekin -= num.dot(self.D_sp[1], self.H_sp[1])

            return Ekin, Epot, Ebar, Exc
        
        else:
            self.vhat_L.integrate(nt, None)
            self.ghat_L.integrate(vHt_g, None)
            return 0.0, 0.0, 0.0, 0.0

    def adjust_residual(self, R_nG, eps_n, s, u, k):
        assert self.domain_overlap >= PROJECTOR_FUNCTION
        if self.domain_overlap == EVERYTHING:
            H_ii = unpack(self.H_sp[s])
            P_ni = self.P_uni[u]
            coefs_ni =  (num.dot(P_ni, H_ii) -
                         num.dot(P_ni * eps_n[:, None], self.setup.O_ii))
            self.pt_i.add(R_nG, coefs_ni, k, communicate=True)
        else:
            self.pt_i.add(R_nG, None, k, communicate=True)
            
    def adjust_residual2(self, pR_G, dR_G, eps, s, k):
        assert self.domain_overlap >= PROJECTOR_FUNCTION
        if self.domain_overlap == EVERYTHING:
            ni = self.get_number_of_partial_waves()
            dP_i = num.zeros(ni, self.typecode)
            self.pt_i.integrate(pR_G, dP_i, k)
        else:
            self.pt_i.integrate(pR_G, None, k)

        if self.domain_overlap == EVERYTHING:
            H_ii = unpack(self.H_sp[s])
            coefs_i = (num.dot(dP_i, H_ii) -
                       num.dot(dP_i * eps, self.setup.O_ii))
            self.pt_i.add(dR_G, coefs_i, k, communicate=True)
        else:
            self.pt_i.add(dR_G, None, k, communicate=True)

    def symmetrize(self, D_aii, map_sa, s):
        D_ii = self.setup.symmetrize(self.a, D_aii, map_sa)
        self.D_sp[s] = pack(D_ii)

    def calculate_force(self, vHt_g, nt_g, vt_G):
        if self.domain_overlap == EVERYTHING:
            lmax = self.setup.lmax
            nk = (3, 9, 22)[lmax]
            # ???? Optimization: do the sum over L before the sum over g and G.
            F_k = num.zeros(nk, num.Float)
            self.ghat_L.integrate(vHt_g, F_k, derivatives=True)
            self.vhat_L.integrate(nt_g, F_k, derivatives=True) 
            
            Q_L = self.Q_L
            F = self.F_c
            F[:] += Q_L[0] * F_k[:3]
            if lmax > 0:
                F += Q_L[1] * F_k[3:6]
                F += Q_L[2] * num.array([F_k[4], F_k[6], F_k[7]])
                F += Q_L[3] * num.array([F_k[5], F_k[7], F_k[8]])
            if lmax > 1:
                f_cm = num.zeros(15, num.Float)
                f_cm[:7] = F_k[9:16]
                f_cm[7] = F_k[10]
                f_cm[8:10] = F_k[16:18]
                f_cm[10] = F_k[10]
                f_cm[11:] = F_k[18:]
                f_cm.shape = (3, 5)
                F += num.dot(f_cm, Q_L[4:])

            # Force from smooth core charge:
            F_k = num.zeros(3, num.Float)
            self.nct.integrate(vt_G, F_k, derivatives=True)
            F += F_k

            # Force from localized potential:
            F_k = num.zeros(3, num.Float)
            self.vbar.integrate(nt_g, F_k, derivatives=True)
            F += F_k

            dF = num.zeros(((lmax + 1)**2, 3), num.Float)
            for neighbor in self.neighbors:
                for i in range(3):
                    dF[:, i] += num.dot(neighbor.dvdr[:, :, i],
                                        neighbor.nucleus.Q_L)
            F += num.dot(self.Q_L, dF)
        else:
            if self.domain_overlap >= COMPENSATION_CHARGE:
                self.ghat_L.integrate(vHt_g, None, derivatives=True)
                self.vhat_L.integrate(nt_g, None, derivatives=True) 
            if self.nct is None:
                self.comm.sum(num.zeros(3, num.Float), self.rank)
            else:
                self.nct.integrate(vt_G, None, derivatives=True)
            if self.domain_overlap >= PROJECTOR_FUNCTION:
                if self.vbar is None:
                    self.pt_i.comm.sum(num.zeros(3, num.Float), self.pt_i.root)
                else:
                    self.vbar.integrate(nt_g, None, derivatives=True)

    def calculate_force_kpoint(self, kpt):
        assert self.domain_overlap >= PROJECTOR_FUNCTION
        f_n = kpt.f_n
        eps_n = kpt.eps_n
        psit_nG = kpt.psit_nG
        s = kpt.s
        u = kpt.u
        k = kpt.k
        if self.domain_overlap == EVERYTHING:
            P_ni = cc(self.P_uni[u])
            nb = P_ni.shape[0]
            H_ii = unpack(self.H_sp[s])
            O_ii = self.setup.O_ii
            nk = self.setup.get_number_of_derivatives()
            F_nk = num.zeros((nb, nk), self.typecode)
            # ???? Optimization: Take the real value of F_nk * P_ni early.
            self.pt_i.integrate(psit_nG, F_nk, k, derivatives=True)
            F_nk *= f_n[:, None]
            F_ik = num.dot(H_ii, num.dot(num.transpose(P_ni), F_nk))
            F_nk *= eps_n[:, None]
            F_ik -= num.dot(O_ii, num.dot(num.transpose(P_ni), F_nk))
            F_ik *= 2.0
            i = 0
            k = 0
            F = self.F_c
            for l in self.setup.l_j:
                f = real(F_ik[i:, k:])
                if l == 0:
                    F += f[0][:3]
                elif l == 1:
                    F[0] += f[0, 0] + f[1, 1] + f[2, 2]
                    F[1] += f[0, 1] + f[1, 3] + f[2, 4]
                    F[2] += f[0, 2] + f[1, 4] + f[2, 5]
                else:
                    F[0] += f[0, 0] + f[1, 1] + f[2, 2] + f[3, 3] + f[4, 4]
                    F[1] += f[0, 5] + f[1, 6] + f[2, 1] + f[3, 7] + f[4, 8]
                    F[2] += f[0, 1] + f[1, 9] + f[2, 10] + f[3, 11] + f[4, 12]
                i += 2 * l + 1
                k += 3 + l * (1 + 2 * l)
        else:
            self.pt_i.integrate(psit_nG, None, k, derivatives=True)
