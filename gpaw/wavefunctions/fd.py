import numpy as np

from gpaw.overlap import Overlap
from gpaw.fd_operators import Laplace
from gpaw.lfc import LocalizedFunctionsCollection as LFC
from gpaw.utilities import unpack
from gpaw.io.tar import TarFileReference
from gpaw.lfc import BasisFunctions
from gpaw.utilities.blas import axpy
from gpaw.utilities.linalg import multi_ax2py
from gpaw.transformers import Transformer
from gpaw.fd_operators import Gradient
from gpaw.band_descriptor import BandDescriptor
from gpaw import extra_parameters
from gpaw.wavefunctions.fdpw import FDPWWaveFunctions
from gpaw.hs_operators import MatrixOperator
from gpaw.preconditioner import Preconditioner

import pycuda.gpuarray as gpuarray
import pycuda.driver as cuda

class FDWaveFunctions(FDPWWaveFunctions):
    def __init__(self, stencil, diagksl, orthoksl, initksl,
                 gd, nvalence, setups, bd,
                 dtype, world, kd, timer=None, cuda=False):
        self.cuda = cuda
        FDPWWaveFunctions.__init__(self, diagksl, orthoksl, initksl,
                                   gd, nvalence, setups, bd,
                                   dtype, world, kd, timer, cuda)

        self.wd = self.gd  # wave function descriptor

        # Kinetic energy operator:
        self.kin = Laplace(self.gd, -0.5, stencil, self.dtype, allocate=False,
                           cuda=self.cuda)

        self.matrixoperator = MatrixOperator(self.bd, self.gd, orthoksl, 
                                             cuda=self.cuda)

        if self.cuda:
            self.nt_G_gpu=None
            from pycuda.elementwise import ElementwiseKernel
            self.axpysquarez = ElementwiseKernel(
                "double f, const pycuda::complex<double> *a, double *b",
                "b[i] += f*(a[i]._M_re*a[i]._M_re+a[i]._M_im*a[i]._M_im)",
                "axpy_squarez",preamble="#include <pycuda-complex.hpp>")
            self.axpysquare = ElementwiseKernel(
                "double f, const double *a, double *b",
                "b[i] += f*a[i]*a[i]",
                "axpy_square")          
            

    def set_setups(self, setups):
        self.pt = LFC(self.gd, [setup.pt_j for setup in setups],
                      self.kpt_comm, dtype=self.dtype, forces=True, 
                      cuda=self.cuda)
        FDPWWaveFunctions.set_setups(self, setups)

    def set_positions(self, spos_ac):
        if not self.kin.is_allocated():
            self.kin.allocate()
        FDPWWaveFunctions.set_positions(self, spos_ac)

    def summary(self, fd):
        fd.write('Mode: Finite-difference\n')
        
    def make_preconditioner(self, block=1):
        return Preconditioner(self.gd, self.kin, self.dtype, block, self.cuda)
    
    def apply_pseudo_hamiltonian(self, kpt, hamiltonian, psit_xG, Htpsit_xG):
        self.kin.apply(psit_xG, Htpsit_xG, kpt.phase_cd)
        hamiltonian.apply_local_potential(psit_xG, Htpsit_xG, kpt.s)

    def add_orbital_density(self, nt_G, kpt, n):
        if self.dtype == float:
            axpy(1.0, kpt.psit_nG[n]**2, nt_G)
        else:
            axpy(1.0, kpt.psit_nG[n].real**2, nt_G)
            axpy(1.0, kpt.psit_nG[n].imag**2, nt_G)

    def add_to_density_from_k_point_with_occupation(self, nt_sG, kpt, f_n, cuda_psit_nG=False):
        # Used in calculation of response part of GLLB-potential

        
        if cuda_psit_nG:
            assert self.cuda
            if self.nt_G_gpu is None:
                self.nt_G_gpu = gpuarray.empty(nt_sG[kpt.s].shape,nt_sG[kpt.s].dtype)
            self.nt_G_gpu.set(nt_sG[kpt.s])
            psit_nG = kpt.psit_nG_gpu
            multi_ax2py(f_n,psit_nG,self.nt_G_gpu)
            #for f, psit_G in zip(f_n, psit_nG):
            #    ax2py(f,psit_G,nt_G)     #axpy(f, psit_G*psit_G, nt_G)
            #if self.dtype == float:
            #    for f, psit_G in zip(f_n, psit_nG):
            #        self.axpysquare(f,psit_G,nt_G)
            #else:
            #    for f, psit_G in zip(f_n, psit_nG):
            #        self.axpysquarez(f,psit_G,nt_G)


        else:
            nt_G = nt_sG[kpt.s]
            psit_nG = kpt.psit_nG
            
            if self.dtype == float:
                for f, psit_G in zip(f_n, psit_nG):
                    axpy(f, psit_G*psit_G, nt_G)  # PyCUDA does not support power?
                    #nt_G+=f*(psit_G**2)
            else:
                for f, psit_G in zip(f_n, psit_nG):
                    axpy(f, psit_G.real*psit_G.real, nt_G)
                    axpy(f, psit_G.imag*psit_G.imag, nt_G)
                    #axpy(f, psit_G.real**2, nt_G)
                    #axpy(f, psit_G.imag**2, nt_G)

        # Hack used in delta-scf calculations:
        if hasattr(kpt, 'c_on'):
            assert self.bd.comm.size == 1
            d_nn = np.zeros((self.bd.mynbands, self.bd.mynbands),
                            dtype=complex)
            for ne, c_n in zip(kpt.ne_o, kpt.c_on):
                d_nn += ne * np.outer(c_n.conj(), c_n)
            for d_n, psi0_G in zip(d_nn, psit_nG):
                for d, psi_G in zip(d_n, psit_nG):
                    if abs(d) > 1.e-12:
                        nt_G += (psi0_G.conj() * d * psi_G).real

        if cuda_psit_nG:
            self.nt_G_gpu.get(nt_sG[kpt.s])

    def calculate_kinetic_energy_density(self, tauct, grad_v):
        assert not hasattr(self.kpt_u[0], 'c_on')
        if isinstance(self.kpt_u[0].psit_nG, TarFileReference):
            raise RuntimeError('Wavefunctions have not been initialized.')

        taut_sG = self.gd.zeros(self.nspins)
        dpsit_G = self.gd.empty(dtype=self.dtype)
        for kpt in self.kpt_u:
            for f, psit_G in zip(kpt.f_n, kpt.psit_nG):
                for v in range(3):
                    grad_v[v](psit_G, dpsit_G)
                    axpy(0.5 * f, abs(dpsit_G)**2, taut_sG[kpt.s])

        self.kpt_comm.sum(taut_sG)
        self.band_comm.sum(taut_sG)
        return taut_sG
        
    def calculate_forces(self, hamiltonian, F_av):
        # Calculate force-contribution from k-points:
        F_av.fill(0.0)
        F_aniv = self.pt.dict(self.bd.mynbands, derivative=True)
        for kpt in self.kpt_u:
            self.pt.derivative(kpt.psit_nG, F_aniv, kpt.q)
            for a, F_niv in F_aniv.items():
                F_niv = F_niv.conj()
                F_niv *= kpt.f_n[:, np.newaxis, np.newaxis]
                dH_ii = unpack(hamiltonian.dH_asp[a][kpt.s])
                P_ni = kpt.P_ani[a]
                F_vii = np.dot(np.dot(F_niv.transpose(), P_ni), dH_ii)
                F_niv *= kpt.eps_n[:, np.newaxis, np.newaxis]
                dO_ii = hamiltonian.setups[a].dO_ii
                F_vii -= np.dot(np.dot(F_niv.transpose(), P_ni), dO_ii)
                F_av[a] += 2 * F_vii.real.trace(0, 1, 2)

            # Hack used in delta-scf calculations:
            if hasattr(kpt, 'c_on'):
                assert self.bd.comm.size == 1
                self.pt.derivative(kpt.psit_nG, F_aniv, kpt.q) #XXX again
                d_nn = np.zeros((self.bd.mynbands, self.bd.mynbands),
                                dtype=complex)
                for ne, c_n in zip(kpt.ne_o, kpt.c_on):
                    d_nn += ne * np.outer(c_n.conj(), c_n)
                for a, F_niv in F_aniv.items():
                    F_niv = F_niv.conj()
                    dH_ii = unpack(hamiltonian.dH_asp[a][kpt.s])
                    Q_ni = np.dot(d_nn, kpt.P_ani[a])
                    F_vii = np.dot(np.dot(F_niv.transpose(), Q_ni), dH_ii)
                    F_niv *= kpt.eps_n[:, np.newaxis, np.newaxis]
                    dO_ii = hamiltonian.setups[a].dO_ii
                    F_vii -= np.dot(np.dot(F_niv.transpose(), Q_ni), dO_ii)
                    F_av[a] += 2 * F_vii.real.trace(0, 1, 2)

        self.bd.comm.sum(F_av, 0)

        if self.bd.comm.rank == 0:
            self.kpt_comm.sum(F_av, 0)

    def cuda_psit_nG_htod(self):
        for kpt in self.kpt_u:
            kpt.cuda_psit_nG_htod()
        
    def cuda_psit_nG_dtoh(self):
        for kpt in self.kpt_u:
            kpt.cuda_psit_nG_dtoh()

    def estimate_memory(self, mem):
        FDPWWaveFunctions.estimate_memory(self, mem)
        self.kin.estimate_memory(mem.subnode('Kinetic operator'))
