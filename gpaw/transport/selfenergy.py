import numpy as np
from gpaw import GPAW
from gpaw.transport.tools import k2r_hs, get_hs, dagger, dot
from gpaw.utilities.lapack import inverse_general, inverse_symmetric
from gpaw.lcao.tools import get_realspace_hs

class LeadSelfEnergy:
    conv = 1e-8 # Convergence criteria for surface Green function
    
    def __init__(self, hs_dii, hs_dij, hs_dim, eta=1e-4):
        self.h_ii, self.s_ii = hs_dii # onsite principal layer
        self.h_ij, self.s_ij = hs_dij # coupling between principal layers
        self.h_im, self.s_im = hs_dim # coupling to the central region
        self.nbf = self.h_im.shape[1] # nbf for the scattering region
        self.eta = eta
        self.energy = None
        self.bias = 0
        self.sigma_mm = np.empty((self.nbf, self.nbf), complex)

    def __call__(self, energy):
        """Return self-energy (sigma) evaluated at specified energy."""
        if self.h_im.dtype == float:
            inv = inverse_symmetric
        if self.h_im.dtype == complex:
            inv = inverse_general
        if energy != self.energy:
            self.energy = energy
            z = energy - self.bias + self.eta * 1.j           
            tau_im = z * self.s_im - self.h_im
            ginv = self.get_sgfinv(energy)
            inv(ginv)
            a_im = dot(ginv, tau_im)
            tau_mi = z * dagger(self.s_im) - dagger(self.h_im)
            self.sigma_mm[:] = dot(tau_mi, a_im)
        return self.sigma_mm

    def set_bias(self, bias):
        self.bias = bias
        
    def get_lambda(self, energy):
        """Return the lambda (aka Gamma) defined by i(S-S^d).

        Here S is the retarded selfenergy, and d denotes the hermitian
        conjugate.
        """
        sigma_mm = self(energy)
        return 1.j * (sigma_mm - dagger(sigma_mm))
        
    def get_sgfinv(self, energy):
        """The inverse of the retarded surface Green function"""
        if self.h_im.dtype == float:
            inv = inverse_symmetric
        if self.h_im.dtype == complex:
            inv = inverse_general
        z = energy - self.bias + self.eta * 1.0j
        
        v_00 = z * dagger(self.s_ii) - dagger(self.h_ii)
        v_11 = v_00.copy()
        v_10 = z * self.s_ij - self.h_ij
        v_01 = z * dagger(self.s_ij) - dagger(self.h_ij)
        delta = self.conv + 1
        n = 0
        while delta > self.conv:
            inv(v_11)
            a = dot(v_11, v_01)
            b = dot(v_11, v_10)
            v_01_dot_b = dot(v_01, b)
            v_00 -= v_01_dot_b
            v_11 -= dot(v_10, a)
            v_11 -= v_01_dot_b
            v_01 = -dot(v_01, a)
            v_10 = -dot(v_10, b)
            delta = np.abs(v_01).max()
            n += 1
        return v_00
    
class CellSelfEnergy:
    def __init__(self, hs_kii, hs_ii, kpts, weight, eta=1e-4):
        self.h_skmm, self.s_kmm = hs_kii
        self.h_smm, self.s_mm = hs_ii
        self.kpts = kpts
        self.weight = weight
        self.bias = 0
        self.energy = None
        self.sigma = None
        
    def set_bias(self, bias):
        self.bias = bias
        
    def __call__(self, energy):
        h_skmm, s_kmm = self.h_skmm, self.s_kmm
        h_mm, s_mm = self.h_smm, self.s_mm        
        ns = h_skmm.shape[0]
        nk = h_skmm.shape[1]
        nb = h_skmm.shape[-1]
        kpts = self.kpts
        weight = self.weight
        
        g_skmm = np.empty(h_skmm.shape, h_skmm.dtype)
        sigma = np.empty([ns, nb, nb], complex)
        inv = inverse_general
        for s in range(ns):
            for k in range(nk):
                g_skmm[s, k] = energy * s_kmm[k] - h_skmm[s, k]
                inv(g_skmm[s, k])
            g_mm = get_realspace_hs(g_skmm, None, kpts, weight)
            inv(g_mm[s])
            sigma[s] = energy * s_mm - h_mm - g_mm[s]
        self.sigma = sigma[0]
        return self.sigma
    
    def get_lambda(self, energy):
        return 1.j * (self.sigma - dagger(self.sigma))
            
        
        
