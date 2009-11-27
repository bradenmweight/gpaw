from math import pi, sqrt
from os.path import isfile

import numpy as np
from ase.units import Hartree

from gpaw.coulomb import CoulombNEW
from gpaw.utilities import pack, unpack, unpack2 
from gpaw.utilities.lapack import diagonalize
from gpaw.response import CHI

class FiniteSys(CHI):
    def __init__(self):
        CHI.__init__(self)


    def get_dipole_strength(self, calc, q, wcut, wmin, wmax, dw, eta=0.2, sigma=2*1e-5):
        """Obtain the dipole strength spectra for a finite system.

        Parameters: 

        n_S: ndarray 
            Pair-orbitals in real space, (1, nS)
        specfunc_wSS: ndarray
            Spectral function, (NwS, nS, nS, dtype = C_knM.type), can be complex128 or float64
        chi0_wSS: ndarray
            The non-interacting density response function, (Nw, nS, nS, dtype=complex)
        kernelRPA_SS (or kernelLDA_SS): ndarray
            Kernel for the finite sys, (nS, nS), it is float64, but can be complex for periodic sys.
        SNonInter_w (or SRPA_w, SLDA_w): ndarray
            Dipole strength function, (Nw)
        """
        e_kn, f_kn, C_knM, orb_MG, spos_ac, nt_G, tmp = (
                self.initialize(calc, q, wcut, wmin, wmax, dw, eta)) 
        if self.HilbertTrans:
            assert tmp.shape == (self.NW, self.nS, self.nS) and tmp.dtype == complex
            chi0_wSS = tmp
        else:
            assert tmp.shape == (self.nkpt, 3)
            bzkpt_kG = tmp

        # Get pair-orbitals in real space
        n_S = self.pair_orbital_Rspace(orb_MG, calc.gd.h_c, calc.wfs.setups, 
                                         calc.wfs.kpt_u[0])

        # Get kernel
        if isfile('kernel.npz'):
            foo = np.load('kernel.npz')
            kernelRPA_SS = foo['KRPA']
            kernelLDA_SS = foo['KLDA']

        else:
            kernelRPA_SS, kernelLDA_SS = self.kernel_finite_sys(nt_G, calc.density.D_asp, orb_MG, 
                        calc.wfs.kpt_u[0], calc.gd, calc.wfs.setups, spos_ac)
            np.savez('kernel.npz',KRPA=kernelRPA_SS,KLDA=kernelLDA_SS)

        # Solve Dyson's equation and Get dipole strength function
        SNonInter_w = np.zeros((self.Nw,3))
        SRPA_w = np.zeros((self.Nw,3))
        SLDA_w = np.zeros((self.Nw,3))
        for iw in range(self.Nw):
            if not self.HilbertTrans:
                chi0_SS = self.calculate_chi0(bzkpt_kG, e_kn, f_kn, C_knM, q, iw*self.dw, eta=eta/Hartree)
            else:
                chi0_SS = chi0_wSS[iw]

            SNonInter_w[iw,:] = self.calculate_dipole_strength(chi0_SS, n_S, iw*self.dw)
            chi_SS = self.solve_Dyson(chi0_SS, kernelRPA_SS)
            SRPA_w[iw,:] = self.calculate_dipole_strength(chi_SS, n_S, iw*self.dw)
            chi_SS = self.solve_Dyson(chi0_SS, kernelLDA_SS)
            SLDA_w[iw,:] = self.calculate_dipole_strength(chi_SS, n_S, iw*self.dw)

        # Solve Casida's equation to get the excitation energies in eV
        eCasidaRPA_s, sCasidaRPA_s = self.solve_casida(e_kn[0], f_kn[0], C_knM[0], kernelRPA_SS, n_S)
        eCasidaLDA_s, sCasidaLDA_s = self.solve_casida(e_kn[0], f_kn[0], C_knM[0], kernelLDA_SS, n_S)

        return SNonInter_w, SRPA_w, SLDA_w, eCasidaRPA_s, eCasidaLDA_s, sCasidaRPA_s, sCasidaLDA_s


    def kernel_finite_sys(self, nt_G, D_asp, orb_MG, kpt, gd, setups, spos_ac):
        """Calculate the Kernel for a finite system. 
    
        The kernel is expressed as, refer to report 4/11/2009, Eq. (18) - (22)::
                                                                            
                     //                 (      1                )
            K      = || dr1 dr2 n (r1 ) | --------  + f  (r1,r2)|  n (r2)
             S1,S2   //          S1     ( |r1 - r2|    xc       )   S2  
                                                                   
        while::

                     ~        ----  a       ~ a
            n (r)  = n (r)  + \    n (r)  - n (r)
             S        S       /___  S        S
                                a
            ~
            n (r)  = phi (r) phi (r)
             S          mu      nu
    
             a       ----          ~ a    ~ a
            n (r)  = \    < phi  | p  > < p  | phi  >  phi (r) phi (r)
             S       /___      mu   i      j      nu      i       j
                      ij
            ~a       ----          ~ a    ~ a           ~       ~  
            n (r)  = \    < phi  | p  > < p  | phi  >  phi (r) phi (r)
             S       /___      mu   i      j      nu      i       j
                      ij

        Note, phi_mu is LCAO orbital, while phi_i or phi_j are partial waves

        Coulomb Kernel: use coulomb.calculate (note it returns the E_coul in eV)

        XC kernel is obtained by::
 
             xc
            K       = < n   | f [n] | n   >  (note, n is the total density)
             S1,S2       S1    xc      S2
                        ~        ~    ~
                    = < n   | f [n] | n   > 
                         S1    xc      S2
                         ----     a        a     a         ~a       ~a    ~a
                      +  \     < n   | f [n ] | n   >  - < n   | f [n ] | n   >
                         /___     S1    xc       S2         S1    xc       S2
                           a

        The second term of the XC kernel can be further evaluated by::

            ---- ----           ~ a     ~ a                     ~ a     ~ a
            \    \     < phi  | p   > < p   | phi  >   < phi  | p   > < p   | phi  >
            /___ /___       mu   i1      i2      nu         mu   i3      i4      nu
              a  i1,i2        1                    1          2                    2
                 i3,i4

                    (  /      a       a         a     a       a 
                  * | | dr phi (r) phi (r)  f [n ] phi (r) phi (r) 
                    ( /       i1      i2     xc       i3      i4

                       /    ~ a     ~ a        ~a   ~ a     ~ a     )
                    - | dr phi (r) phi (r)  f [n ] phi (r) phi (r)  |
                      /       i1      i2     xc       i3      i4    )

        The method four_phi_integrals calculate the () term in the above equation
        """

        Kcoul_SS = np.zeros((self.nS, self.nS))
        Kxc_SS = np.zeros_like(Kcoul_SS)
        P1_ap = {}
        P2_ap = {}
        J_II = {}

        fxc_G = self.fxc(nt_G)  # nt_G contains core density
        for a, D_sp in D_asp.items():
            J_pp = setups[a].xc_correction.four_phi_integrals(D_sp, self.fxc)
            ni = setups[a].ni
            J_II[a] = np.zeros((ni*ni, ni*ni))   
            nii = J_pp.shape[0]
            J_pI = np.zeros((nii, ni*ni))
            for ip, J_p in enumerate(J_pp):
                J_pI[ip] = unpack2(J_p).ravel() # D_sp uses pack
            for ii in range(ni*ni):
                J_II[a][:, ii] = unpack2(J_pI[:, ii].copy()).ravel()

        coulomb = CoulombNEW(gd, setups, spos_ac)
        for n in range(self.nLCAO):
            for m in range(self.nLCAO):
                nt1_G = orb_MG[n] * orb_MG[m] 
                for a, P_Mi in kpt.P_aMi.items():
                    D_ii = np.outer(P_Mi[n].conj(), P_Mi[m])
                    P1_ap[a] = pack(D_ii, tolerance=1e30)
                for p in range(self.nLCAO):
                    for q in range(self.nLCAO):
                        nt2_G = orb_MG[p] * orb_MG[q]
                        # Coulomb Kernel
                        for a, P_Mi in kpt.P_aMi.items():
                            D_ii = np.outer(P_Mi[p].conj(), P_Mi[q])
                            P2_ap[a] = pack(D_ii, tolerance=1e30)
                        Kcoul_SS[self.nLCAO*n+m, self.nLCAO*p+q] = coulomb.calculate(
                                    nt1_G, nt2_G, P1_ap, P2_ap)
                        # XC Kernel
                        Kxc_SS[self.nLCAO*n+m, self.nLCAO*p+q] = gd.integrate(nt1_G*fxc_G*nt2_G)

                        for a, P_Mi in kpt.P_aMi.items():
                            P1_I = np.outer(P_Mi[n].conj(), P_Mi[m]).ravel()                            
                            P2_I = np.outer(P_Mi[p].conj(), P_Mi[q]).ravel() 
                            Kxc_SS[self.nLCAO*n+m, self.nLCAO*p+q] += (
                                    np.inner(np.inner(P1_I, J_II[a]), P2_I) )

            print 'finished', n, 'cycle', ' (max: nLCAO = ', self.nLCAO, ')'
        tmp = Kcoul_SS / Hartree

        return tmp, tmp + Kxc_SS


    def calculate_dipole_strength(self, chi_SS, n_S, omega):
        """Calculate dipole strength for a particular omega.

        The dipole strength is obtained by (atomic unit)::

                    2w
            S(w) = ---- Im alpha(w) , 
                    pi

        while alpha is the dynamical polarizability defined as::

                         //
            alpha(w) = - || dr dr' r chi(r,r',w) r'
                        //
                         //          ----               
                     = - || dr dr' r \    n (r) chi (w) n (r') r'
                        //           /___  S       SS'   S'
                                      SS'

        The pair density is defined as::

            n (r) = phi (r) * phi (r) 
             S         mu        nu
     
                      ----          ~ a    ~ a          (                    ~       ~      )
                    + \    < phi  | p  > < p  | phi  >  | phi (r) phi (r) - phi (r) phi (r) |
                      /___      mu   i      j      nu   (    i       j         i       j    )
                       ij       

        As a result::

                         ----   /                    /
            alpha(w) = - \     | dr r n (r) chi (w) | dr' r' n (r')
                         /___ /        S       SS'  /         S'
                          SS'
                         ----
                     = - \    n  chi (w) n  
                         /___  S    SS'   S'
                          SS'
 
        where n_S is defined in pair_orbital_Rspace
        """

        alpha = np.zeros(3, dtype=complex)
        for i in range(3):
            alpha[i] = - np.dot( np.dot(n_S[:,i], chi_SS), n_S[:,i]) 

        S = 2. * omega / pi * np.imag(alpha)

        return S


    def pair_orbital_Rspace(self, orb_MG, h_c, setups, kpt):
        """Calculate pair LCAO orbital in real space. 

        The pair density is defined as::
             
                   /
            n   =  | dr  r  phi (r) phi (r) 
             S    /           mu      nu
                      ----          ~ a    ~ a          /      (                    ~       ~      )
                    + \    < phi  | p  > < p  | phi  >  | dr r | phi (r) phi (r) - phi (r) phi (r) |
                      /___      mu   i      j      nu  /       (    i       j         i       j    )
                       ij       

        Parameters:

        orb: ndarray
            LCAO orbital on the grid, (nband, Nx, Ny, Nz)
        Delta_pL: ndarray
            L = 1, 2, 3 corresponds to y, z, x, refer to c/bmgs/sharmonic.c
        """

        N_gd = orb_MG.shape[1:4] # number of grid points
        r = np.zeros((N_gd))        
        n_MM = np.zeros((self.nLCAO, self.nLCAO))
        n_S = np.zeros((self.nS, 3))
        tmp =  sqrt(4. * pi / 3.)    
        Li = np.array([3, 1, 2])

        for ix in range(3): # loop over x, y, z axis
            phi_I={}
            for a in range(len(setups)):
                phi_p = setups[a].Delta_pL[:,Li[ix]].copy()
                phi_I[a] = unpack(phi_p).ravel()
       
            for i in range(N_gd[0]):
                for j in range(N_gd[1]):
                    for k in range(N_gd[2]):
                        if ix == 0:
                            r[i,j,k] = i*h_c[0]
                        elif ix == 1:
                            r[i,j,k] = j*h_c[1] 
                        else:
                            r[i,j,k] = k*h_c[2]
    
            for mu in range(self.nLCAO):
                for nu in range(self.nLCAO):  
                    n_MM[mu,nu] = np.sum(orb_MG[mu] * orb_MG[nu] * r)
                    for a, P_Mi in kpt.P_aMi.items():
                        P_I = np.outer(P_Mi[mu], P_Mi[nu]).ravel()
                        n_MM[mu,nu] += np.sum(P_I * phi_I[a]) * tmp
            n_S[:,ix] = np.reshape(n_MM, self.nS)

        return n_S


    def solve_casida(self, e_n, f_n, C_nM, kernel_SS, n_S):
        """Solve Casida's equation with input from LCAO calculations (nspin = 1).

        The Casida matrix equation is written as::

                            2 
            Omega F  = omega  F
                   I        I  I

        while the Omega matrix is defined as::

                                      2          ---------------  
            Omega   = delta  delta   e    + 2   / f   e  f   e   K  
                ss'        ik     jq  s       \/   s   s  s'  s'  ss'

        Note, s(s') is a combined index for ij (kq)

        The kernel is obtained from ::

                    ----
            K    =  \    C      C      C      C      K
             ss'    /___  i,mu   j,nu   k,mu   q,nu   S S
                     S S      1      1      2      2   1 2
                      1 2

        Parameters:

        i (or k): integer
            Index for occupied states
        j (or q): integer
            Index for unoccupied states
        s (or s'): integer
            Combined index for ij (or kq)
        S (or S'): integer
            Combined index for mu,nu (or mu',nu')
        C_nM: ndarray
            The LCAO coefficients at kpt=0, (nband, nLCAO)
        omega_I: ndarray
            Excitation energies
        """

        # Count number of occupied and unoccupied states pairs
        Nocc = 0
        Nunocc = 0
        for n in range(self.nband):
            if f_n[n] > 0:
                Nocc += 1
            else: 
                Nunocc +=1
        if Nocc + Nunocc != self.nband:
            raise ValueError('Nocc + Nunocc != nband')
        npair = Nocc * Nunocc

        # calculate the factor before the K_ij,kq matrix
        e_s = np.zeros(npair)
        f_s = np.zeros_like(e_s)

        ipair = 0
        for i in range(Nocc):
            for j in range(Nocc, self.nband):
                e_s[ipair] = e_n[j] - e_n[i] # s: ij pair
                f_s[ipair] = f_n[i] - f_n[j]
                ipair += 1

        fe_ss = np.outer(e_s * f_s, e_s * f_s).ravel()
        fe_ss = (np.array([2. * sqrt(fe_ss[i]) for i in range(npair**2)])).reshape(npair, npair)
 
        # calculate kernel K_ij,kq
        npair1 = 0
        npair2 = 0
        kernel_ss = np.zeros((npair, npair))

        for i in range(Nocc):
            for j in range(Nocc, self.nband):
                C1_S = np.outer(C_nM[i], C_nM[j]).ravel() # S: mu nu pair
                for k in range(Nocc):
                    for q in range(Nocc, self.nband):
                        C2_S = np.outer(C_nM[k], C_nM[q]).ravel() # S: mu nu pair
                        kernel_ss[npair1, npair2] = np.inner(np.inner(C1_S, kernel_SS), C2_S)
                        npair2 += 1
                npair1 += 1
                npair2 = 0

        kernel_ss *= fe_ss

        # add the delta matrix to obtain the Omega matrix
        delta_ss = np.eye(npair,npair)
        for i in range(npair):
            delta_ss[i,i] *= e_s[i]**2

        Omega_ss = kernel_ss + delta_ss

        # diagonalize the Omega matrix
        eExcitation_s = np.zeros(npair)
        diagonalize(Omega_ss, eExcitation_s)

        # get the excitation energies in Hartree
        eExcitation_s = np.array([sqrt(eExcitation_s[i]) for i in range(npair)])
        
        # get the dipole strength 
        ipair = 0
        mu_s = np.zeros((npair, 3))
        for i in range(Nocc):
            for j in range(Nocc, self.nband):
                for ix in range(3): # x,y,z three directions
                    mu_s[ipair, ix] = np.inner(np.outer(C_nM[i], C_nM[j]).ravel(), n_S[:, ix])
                ipair += 1

        fe_s = np.array([ sqrt(f_s[i] * e_s[i]) for i in range(npair)])

        DipoleStrength = np.zeros((npair, 3))
        for s1 in range(npair):
            FI_s = Omega_ss[s1]            
            DipoleStrength[s1] = np.array([ 2. * ((mu_s[:, ix] * fe_s * FI_s).sum())**2 for ix in range(3) ])

        return eExcitation_s * Hartree, DipoleStrength


