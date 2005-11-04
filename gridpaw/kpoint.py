# Copyright (C) 2003  CAMP
# Please see the accompanying LICENSE file for further information.

"""This module defines a ``KPoint`` class."""

from math import pi
from cmath import exp

import Numeric as num
import LinearAlgebra as linalg

from gridpaw.utilities.blas import axpy, rk, r2k, gemm
from gridpaw.utilities.complex import cc, real
from gridpaw.utilities.lapack import diagonalize
from gridpaw.utilities import unpack
from gridpaw.utilities.timing import Timer
from gridpaw.operators import Gradient


class KPoint:
    """Class for a singel **k**-point.

    The ``KPoint`` class takes care of all wave functions for a
    certain **k**-point and a certain spin."""
    
    def __init__(self, gd, weight, s, k, u, k_c, typecode):
        """Construct **k**-point object.

        Parameters:
         ============ ==========================================
         ``gd``       Grid-descriptor.
         ``weight``   Weight of this **k**-point.
         ``s``        Spin-index.
         ``k``        **k**-point index.
         ``u``        Combined spin and **k**-point index.
         ``k_c``      Scaled **k**-point vector.
         ``typecode`` Data type for wave functions (``Float`` or
                      ``Complex``).
         ============ ==========================================
        """

        self.gd = gd
        self.weight = weight
        self.typecode = typecode
        
        self.phase_cd = num.ones((3, 2), num.Complex)
        if typecode == num.Float:
            # Gamma-point calculation:
            self.k_c = None
        else:
            sdisp_cd = self.gd.domain.sdisp_cd
            for c in range(3):
                for d in range(2):
                    self.phase_cd[c, d] = exp(2j * pi *
                                              sdisp_cd[c, d] * k_c[c])
            self.k_c = k_c

        self.s = s  # spin index
        self.k = k  # k-point index
        self.u = u  # combined spin and k-point index

        # Which CPU does overlap-matrix Cholesky-decomposition and
        # Hamiltonian-matrix diagonalization?
        self.comm = self.gd.comm
        self.root = u % self.comm.size
        
        self.psit_nG = None
        self.Htpsit_nG = None

        self.timer = Timer()
        
    def allocate(self, nbands):
        """Allocate arrays."""
        self.nbands = nbands
        self.eps_n = num.zeros(nbands, num.Float)
        self.f_n = num.ones(nbands, num.Float) * self.weight
        self.H_nn = num.zeros((nbands, nbands), self.typecode)
        self.S_nn = num.zeros((nbands, nbands), self.typecode)

    def diagonalize(self, kin, vt_sG, my_nuclei, nbands):
        """Subspace diagonalization of wave functions.

        First, the Hamiltonian (defined by ``kin``, ``vt_sG``, and
        ``my_nuclei``) is applied to the wave functions, then the
        ``H_nn`` matrix is calculated and diagonalized, and finally,
        the wave functions are rotated.  Also the projections
        ``P_uni`` (an attribute of the nuclei) are rotated.

        If this is the first iteration and we are starting from atomic
        orbitals, then the desired number of bands (``nbands``) will
        most likely differ from the number of current atomic orbitals
        (``self.nbands``).  If this is the case, then new arrays are
        allocated:

        * Too many bands: The bands with the lowest eigenvalues are
          used.
        * Too few bands: Extra wave functions calculated as the
          derivative of the wave functions with respect to the
          *x*-coordinate.
        """

        kin.apply(self.psit_nG, self.Htpsit_nG, self.phase_cd)
        self.Htpsit_nG += self.psit_nG * vt_sG[self.s]
        r2k(0.5 * self.gd.dv, self.psit_nG, self.Htpsit_nG, 0.0, self.H_nn)

        for nucleus in my_nuclei:
            P_ni = nucleus.P_uni[self.u]
            self.H_nn += num.dot(P_ni, num.dot(unpack(nucleus.H_sp[self.s]),
                                               cc(num.transpose(P_ni))))

        self.comm.sum(self.H_nn, self.root)

        yield None
        
        if self.comm.rank == self.root:
            info = diagonalize(self.H_nn, self.eps_n)
            if info != 0:
                raise RuntimeError, 'Very Bad!!'
        
        yield None
        
        self.comm.broadcast(self.H_nn, self.root)
        self.comm.broadcast(self.eps_n, self.root)

        # Rotate psit_nG:
        # We should block this so that we can use a smaller temp !!!!!
        temp = num.array(self.psit_nG)
        gemm(1.0, temp, self.H_nn, 0.0, self.psit_nG)
        
        # Rotate Htpsit_nG:
        temp[:] = self.Htpsit_nG
        gemm(1.0, temp, self.H_nn, 0.0, self.Htpsit_nG)
        
        # Rotate P_ani:
        for nucleus in my_nuclei:
            P_ni = nucleus.P_uni[self.u]
            temp_ni = P_ni.copy()
            gemm(1.0, temp_ni, self.H_nn, 0.0, P_ni)
        
        if nbands != self.nbands:
            nao = self.nbands  # number of atomic orbitals
            nmin = min(nao, nbands)
            
            tmp_nG = self.psit_nG
            self.psit_nG = self.gd.new_array(nbands, self.typecode)
            self.psit_nG[:nmin] = tmp_nG[:nmin]

            tmp_nG = self.Htpsit_nG
            self.Htpsit_nG = self.gd.new_array(nbands, self.typecode)
            self.Htpsit_nG[:nmin] = tmp_nG[:nmin]
            del tmp_nG

            tmp_n = self.eps_n
            self.allocate(nbands)
            self.eps_n[:nmin] = tmp_n[:nmin]

            extra = nbands - nao
            if extra > 0:
                self.eps_n[nao:] = self.eps_n[nao - 1] + 0.5
                slice_nG = self.psit_nG[nao:]
                ddx = Gradient(self.gd, 0, typecode=self.typecode).apply
                ddx(self.psit_nG[:extra], slice_nG, self.phase_cd)
        
        yield None
        
    def calculate_residuals(self, p_nuclei):
        """Calculate wave function residuals.

        On entry, ``Htpsit_nG`` contains the soft part of the
        Hamiltonian applied to the wave functions.  After this call,
        ``Htpsit_nG`` holds the residuals::

          ^  ~        ^  ~   
          H psi - eps S psi =
                                _ 
              ~  ~         ~   \   ~a    a           a     ~a   ~
              H psi - eps psi + )  p  (dH    - eps dS    )<p  |psi>
                               /_   i1   i1i2        i1i2   i2
                              ai1i2

                                
        The size of the residuals is returned."""
        
        R_nG = self.Htpsit_nG
        # optimize XXX 
        for R_G, eps, psit_G in zip(R_nG, self.eps_n, self.psit_nG):
            R_G -= eps * psit_G
        
        for nucleus in p_nuclei:
            nucleus.adjust_residual(R_nG, self.eps_n, self.s, self.u, self.k)

        error = 0.0
        for R_G, f in zip(R_nG, self.f_n):
            error += f * real(num.dot(cc(R_G).flat, R_G.flat))

        return error
        
    def orthonormalize(self, my_nuclei):
        """Orthogonalize wave functions."""
        S_nn = self.S_nn

        # Fill in the lower triangle:
        rk(self.gd.dv, self.psit_nG, 0.0, S_nn)
        
        for nucleus in my_nuclei:
            P_ni = nucleus.P_uni[self.u]
            S_nn += num.dot(P_ni,
                            cc(num.innerproduct(nucleus.setup.O_ii, P_ni)))
        
        self.comm.sum(S_nn, self.root)

        yield None

        if self.comm.rank == self.root:
            # inverse returns a non-contigous matrix - grrrr!  That is
            # why there is a copy.  Should be optimized with a
            # different lapack call to invert a triangular matrix XXXXX
            S_nn[:] = linalg.inverse(
                linalg.cholesky_decomposition(S_nn)).copy()

        yield None

        self.comm.broadcast(S_nn, self.root)

        # This step will overwrite the Htpsit_nG array!
        gemm(1.0, self.psit_nG, S_nn, 0.0, self.Htpsit_nG)
        self.psit_nG, self.Htpsit_nG = self.Htpsit_nG, self.psit_nG  # swap

        for nucleus in my_nuclei:
            P_ni = nucleus.P_uni[self.u]
            gemm(1.0, P_ni.copy(), S_nn, 0.0, P_ni)

        yield None

    def add_to_density(self, nt_G):
        """Add contribution to pseudo electron-density."""
        if self.typecode is num.Float:
            for psit_G, f in zip(self.psit_nG, self.f_n):
                nt_G += f * psit_G**2
        else:
            for psit_G, f in zip(self.psit_nG, self.f_n):
                nt_G += f * (psit_G * num.conjugate(psit_G)).real
                
    def rmm_diis(self, p_nuclei, preconditioner, kin, vt_sG):
        """Improve the wave functions.

        Take two steps along the preconditioned residuals.  Step
        lengths are optimized for the first step and reused for the
        seconf."""
        
        vt_G = vt_sG[self.s]
        for n in range(self.nbands):
            R_G = self.Htpsit_nG[n]

            dR_G = num.zeros(R_G.shape, self.typecode)

            pR_G = preconditioner(R_G, self.phase_cd, self.psit_nG[n],
                                  self.k_c)
            
            kin.apply(pR_G, dR_G, self.phase_cd)

            dR_G += vt_G * pR_G

            dR_G -= self.eps_n[n] * pR_G

            for nucleus in p_nuclei:
                nucleus.adjust_residual2(pR_G, dR_G, self.eps_n[n],
                                         self.s, self.k)
            
            RdR = self.comm.sum(real(num.dot(cc(R_G).flat, dR_G.flat)))
            dRdR = self.comm.sum(real(num.dot(cc(dR_G.flat), dR_G.flat)))
            lam = -RdR / dRdR

            R_G *= 2.0 * lam
            axpy(lam**2, dR_G, R_G)
##            R_G += lam**2 * dR_G
            self.psit_nG[n] += preconditioner(R_G, self.phase_cd,
                                              self.psit_nG[n], self.k_c)

    def create_atomic_orbitals(self, nao, nuclei):
        """Initialize the wave functions from atomic orbitals.

        Create ``nao`` atomic orbitals."""
        
        # Allocate space for wave functions, occupation numbers,
        # eigenvalues and projections:
        self.allocate(nao)
        self.psit_nG = self.gd.new_array(nao, self.typecode)
        self.Htpsit_nG = self.gd.new_array(nao, self.typecode)
        
        # fill in the atomic orbitals:
        nao0 = 0
        for nucleus in nuclei:
            nao1 = nao0 + nucleus.get_number_of_atomic_orbitals()
            nucleus.create_atomic_orbitals(self.psit_nG[nao0:nao1], self.k)
            nao0 = nao1
        assert nao0 == nao
