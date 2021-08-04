from math import pi

import numpy as np
from ase.utils import seterr
from gpaw.pw.descriptor import PWDescriptor
from gpaw.typing import Array1D
from scipy.special import erf


class ReciprocalSpacePoissonSolver:
    def __init__(self,
                 pd: PWDescriptor,
                 realpbc_c: Array1D):
        self.pd = pd
        self.realpbc_c = realpbc_c
        self.G2_q = pd.G2_qG[0].copy()
        if pd.gd.comm.rank == 0:
            # Avoid division by zero:
            self.G2_q[0] = 1.0

    def initialize(self):
        pass

    def get_stencil(self):
        return '????'

    def estimate_memory(self, mem):
        pass

    def todict(self):
        return {}

    def solve(self,
              vHt_q: Array1D,
              rhot_q: Array1D) -> float:
        vHt_q[:] = 4 * pi * rhot_q
        vHt_q /= self.G2_q
        epot = 0.5 * self.pd.integrate(vHt_q, rhot_q)
        return epot


class ChargedReciprocalSpacePoissonSolver(ReciprocalSpacePoissonSolver):
    def __init__(self,
                 pd: PWDescriptor,
                 realpbc_c: Array1D,
                 charge: float,
                 eps: float = 1e-5):
        assert not realpbc_c.any()
        ReciprocalSpacePoissonSolver.__init__(self, pd, realpbc_c)
        self.charge = charge
        # Shortest distance from center to edge of cell:
        self.rcut = 0.5 / (pd.gd.icell_cv**2).sum(axis=1).max()**0.5
        self.alpha = -self.rcut**-2 * np.log10(eps)
        center_v = pd.gd.cell_cv.sum(axis=0) / 2
        G2_q = pd.G2_qG[0]
        G_qv = pd.get_reciprocal_vectors()
        self.charge_q = np.exp(-1 / (4 * self.alpha) * G2_q +
                               1j * (G_qv @ center_v)) * charge
        R_Rv = pd.gd.get_grid_point_coordinates().transpose((1, 2, 3, 0))
        d_R = ((R_Rv - center_v)**2).sum(axis=3)**0.5

        with seterr(invalid='ignore'):
            potential_R = erf(self.alpha**0.5 * d_R) / d_R
        if ((pd.gd.N_c % 2) == 0).all():
            potential_R[tuple(pd.gd.N_c // 2)] = (4 * self.alpha / pi)**0.5

        c_R = pd.ifft(self.charge_q)
        c_R /= pd.gd.dv
        print(self.alpha, self.rcut)

        print(pd.gd.integrate(c_R), pd.gd.dv, c_R[20,20,20],
              (self.alpha / pi)**1.5)
        print(pd.fft(c_R)[0])
        print(pd.gd.integrate(c_R, potential_R))
        self.potential_q = charge / pd.gd.dv * pd.fft(potential_R)

    def solve(self,
              vHt_q: Array1D,
              rhot_q: Array1D) -> float:
        neutral_q = rhot_q + self.charge_q
        if self.pd.gd.comm.rank == 0:
            assert abs(neutral_q[0]) < 1e-12
        vHt_q[:] = 4 * pi * neutral_q
        vHt_q /= self.G2_q
        vHt_q -= self.potential_q
        epot = 0.5 * self.pd.integrate(vHt_q, neutral_q)
        print(epot)
        epot -= self.pd.integrate(self.potential_q, rhot_q)
        print(epot)
        epot -= (self.alpha / 2 / pi)**0.5
        print(epot)
        return epot
