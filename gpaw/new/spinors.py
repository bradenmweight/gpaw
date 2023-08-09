from __future__ import annotations

import numpy as np
from gpaw.core.domain import Domain
from gpaw.core.plane_waves import PlaneWaves
from gpaw.typing import Vector


class SpinorWaveFunctionDescriptor(Domain):
    def __init__(self,
                 pw: PlaneWaves,
                 qspiral_v: Vector = None):
        self.pw = pw
        self.qspiral_v = (np.asarray(qspiral_v) if qspiral_v is not None else
                          None)
        Domain.__init__(self, pw.cell_cv, pw.pbc_c, pw.kpt_c, pw.comm,
                        complex)
        self.myshape = (2,) + pw.myshape
        self.itemsize = pw.itemsize
        self.shape = (2,) + pw.shape
        self.dv = pw.dv

    def __repr__(self):
        return f'{self.__class__.__name__}({self.pw}, {self.qspiral_v})'

    def new(self, *, kpt):
        pw = self.pw.new(kpt=kpt)
        pw.qspiral_v = self.qspiral_v
        return SpinorWaveFunctionDescriptor(pw, self.qspiral_v)

    def empty(self, nbands, band_comm, xp=None):
        return self.pw.empty((nbands, 2), band_comm)

    def global_shape(self) -> tuple[int, ...]:
        return (2,) + self.pw.global_shape()

    def indices(self, size):
        return self.pw.indices(size)
