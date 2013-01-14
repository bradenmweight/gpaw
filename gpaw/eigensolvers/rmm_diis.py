"""Module defining  ``Eigensolver`` classes."""

import numpy as np

from gpaw.utilities.blas import axpy, dotc
from gpaw.utilities.mblas import multi_axpy, multi_scal, multi_dotc
from gpaw.eigensolvers.eigensolver import Eigensolver
from gpaw.utilities import unpack
from gpaw.mpi import run



class RMM_DIIS(Eigensolver):
    """RMM-DIIS eigensolver

    It is expected that the trial wave functions are orthonormal
    and the integrals of projector functions and wave functions
    ``nucleus.P_uni`` are already calculated

    Solution steps are:

    * Subspace diagonalization
    * Calculation of residuals
    * Improvement of wave functions:  psi' = psi + lambda PR + lambda PR'
    * Orthonormalization"""

    def __init__(self, keep_htpsit=True, blocksize=1, cuda=False):
        Eigensolver.__init__(self, keep_htpsit, blocksize, cuda)

    def iterate_one_k_point(self, hamiltonian, wfs, kpt):
        """Do a single RMM-DIIS iteration for the kpoint"""

        self.subspace_diagonalize(hamiltonian, wfs, kpt)

        self.timer.start('RMM-DIIS')

        if self.cuda:
            psit_nG = kpt.psit_nG_gpu
        else:
            psit_nG = kpt.psit_nG


        B = self.blocksize
        dR_nG = self.operator.suggest_temporary_buffer(wfs.dtype, self.cuda)
            
        P_axi = wfs.pt.dict(B)

        if self.keep_htpsit:
            R_nG = self.Htpsit_nG
        else:
            print "no keep htpsit"
            R_nG = self.gd.empty(wfs.bd.mynbands, wfs.dtype, cuda=self.cuda)
            wfs.apply_pseudo_hamiltonian(kpt, hamiltonian, psit_nG, R_nG)
            wfs.pt.integrate(psit_nG, kpt.P_ani, kpt.q)
            
        self.calculate_residuals(kpt, wfs, hamiltonian, psit_nG,
                                 kpt.P_ani, kpt.eps_n, R_nG)
            
        if kpt.f_n is None:
            weight = kpt.weight
        else:
            weight = kpt.f_n
            
        if self.nbands_converge != 'occupied':
            for n in range(0, wfs.bd.mynbands):
                if wfs.bd.global_index(n) < self.nbands_converge:
                    weight[n] = kpt.weight
                else:
                    weight[n] = 0.0
                            
        error = sum(weight * multi_dotc(R_nG, R_nG).real)
        
        for n1 in range(0, wfs.bd.mynbands, B):
            n2 = n1 + B
            if n2 > wfs.bd.mynbands:
                n2 = wfs.bd.mynbands
                B = n2 - n1
                P_axi = dict([(a, P_xi[:B]) for a, P_xi in P_axi.items()])
                dR_xG = dR_xG[:B]
                
            n_x = range(n1, n2)
            
            R_xG = R_nG[n1:n2]
            dR_xG = dR_nG[n1:n2]
                
            # Precondition the residual:
            self.timer.start('precondition')
            dpsit_xG = self.preconditioner(R_xG, kpt)
            self.timer.stop('precondition')

            # Calculate the residual of dpsit_G, dR_G = (H - e S) dpsit_G:
            wfs.apply_pseudo_hamiltonian(kpt, hamiltonian, dpsit_xG, dR_xG)
            wfs.pt.integrate(dpsit_xG, P_axi, kpt.q)
            self.calculate_residuals(kpt, wfs, hamiltonian, dpsit_xG,
                                     P_axi, kpt.eps_n[n_x], dR_xG, n_x,
                                     calculate_change=True)
            
        # Find lam that minimizes the norm of R'_G = R_G + lam dR_G
        RdR_n=np.array(multi_dotc(R_nG, dR_nG).real)
        dRdR_n=np.array(multi_dotc(dR_nG, dR_nG).real)

        self.gd.comm.sum(RdR_n)
        self.gd.comm.sum(dRdR_n)

        lam_n = -RdR_n / dRdR_n
        # Calculate new psi'_G = psi_G + lam pR_G + lam pR'_G
        #                      = psi_G + p(2 lam R_G + lam**2 dR_G)
        multi_scal(2.0*lam_n,R_nG)
        multi_axpy(lam_n**2,dR_nG,R_nG)
        #for lam, R_G, dR_G in zip(lam_x, R_xG, dR_xG):
        #    R_G *= 2.0 * lam
        #    axpy(lam**2, dR_G, R_G)  # R_G += lam**2 * dR_G
        self.timer.start('precondition')
        for n1 in range(0, wfs.bd.mynbands, self.blocksize):
            # XXX GPUarray does not support properly multi-d slicing
            n2 = min(n1+self.blocksize, wfs.bd.mynbands)
            psit_G = psit_nG[n1:n2]
            R_xG = R_nG[n1:n2]
            psit_G += self.preconditioner(R_xG, kpt)
        self.timer.stop('precondition')
            
        self.timer.stop('RMM-DIIS')
        error = self.gd.comm.sum(error)
        return error
