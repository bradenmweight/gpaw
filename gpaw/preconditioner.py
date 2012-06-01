# Copyright (C) 2003  CAMP
# Please see the accompanying LICENSE file for further information.

from math import pi

import numpy as np

from gpaw.transformers import Transformer
from gpaw.fd_operators import Laplace
from gpaw.fd_operators import FDOperator

from gpaw.utilities.blas import axpy,scal
from gpaw.utilities.linalg import change_sign

import pycuda.gpuarray as gpuarray
import pycuda.driver as cuda

class Preconditioner:
    def __init__(self, gd0, kin0, dtype=float, block=1, cuda=False):
        gd1 = gd0.coarsen()
        gd2 = gd1.coarsen()
        self.kin0 = kin0

        self.cuda=cuda

        self.kin1 = Laplace(gd1, -0.5, 1, dtype, cuda=self.cuda)
        self.kin2 = Laplace(gd2, -0.5, 1, dtype, cuda=self.cuda)

        self.scratch0 = gd0.zeros((2, block), dtype, False,cuda=self.cuda)
        self.scratch1 = gd1.zeros((3, block),dtype, False,cuda=self.cuda)
        self.scratch2 = gd2.zeros((3, block), dtype, False,cuda=self.cuda)

        self.step = 0.66666666 / kin0.get_diagonal_element()

        self.restrictor_object0 = Transformer(gd0, gd1, 1,dtype, False,cuda=self.cuda)
        self.restrictor_object1 = Transformer(gd1, gd2, 1, dtype, False,cuda=self.cuda)
        self.interpolator_object2 = Transformer(gd2, gd1, 1, dtype, False,cuda=self.cuda)
        self.interpolator_object1 = Transformer(gd1, gd0, 1, dtype, False,cuda=self.cuda)
        self.restrictor0 = self.restrictor_object0.apply
        self.restrictor1 = self.restrictor_object1.apply
        self.interpolator2 = self.interpolator_object2.apply
        self.interpolator1 = self.interpolator_object1.apply
        self.allocated = False

    def allocate(self):
        assert not self.allocated
        for transformer in [self.restrictor_object0,
                            self.restrictor_object1,
                            self.interpolator_object2,
                            self.interpolator_object1]:
            transformer.allocate()
        self.allocated = True
        
    def __call__(self, residuals, kpt):
        nb = len(residuals) # number of bands
        phases = kpt.phase_cd
        step = self.step
        if self.cuda and not isinstance(residuals, gpuarray.GPUArray):
            local_residuals=gpuarray.to_gpu(residuals)
        else:
            local_residuals=residuals

        if self.cuda:
            # XXX GPUarray does not support properly multi-d slicing
            d0, q0 = gpuarray.GPUArray(shape=(self.scratch0.shape[0],)+
                                       (nb,)+self.scratch0.shape[2:],
                                       dtype=self.scratch0.dtype,
                                       allocator=self.scratch0.allocator,
                                       base=self.scratch0,
                                       gpudata=self.scratch0.gpudata)
            r1, d1, q1 = gpuarray.GPUArray(shape=(self.scratch1.shape[0],)+
                                           (nb,)+self.scratch1.shape[2:],
                                           dtype=self.scratch1.dtype,
                                           allocator=self.scratch1.allocator,
                                           base=self.scratch1,
                                           gpudata=self.scratch1.gpudata)
            
            r2, d2, q2 = gpuarray.GPUArray(shape=(self.scratch2.shape[0],)+
                                           (nb,)+self.scratch2.shape[2:],
                                           dtype=self.scratch2.dtype,
                                           allocator=self.scratch2.allocator,
                                           base=self.scratch2,
                                           gpudata=self.scratch2.gpudata)
            self.restrictor0(local_residuals, r1, phases)
            cuda.memcpy_dtod(d1.gpudata, r1.gpudata, r1.nbytes)
            scal(-4 * step, d1)
            self.kin1.apply(d1, q1, phases)
            q1 += r1
            self.restrictor1(q1, r2, phases)
            cuda.memcpy_dtod(d2.gpudata, r2.gpudata, r2.nbytes)
            scal(16 * step, d2)
            self.kin2.apply(d2, q2, phases)
            q2 -= r2
            axpy(-16*step, q2, d2)
            self.interpolator2(d2, q1, phases)
            d1 -= q1
            self.kin1.apply(d1, q1, phases)
            q1 += r1
            axpy(-4*step, q1, d1)
            change_sign(d1)
            self.interpolator1(d1, d0, phases)
            self.kin0.apply(d0, q0, phases)
            q0 -= local_residuals
            axpy(-step, q0, d0)  # d0 -= step * q0
            change_sign(d0)
        else:
            d0, q0 = self.scratch0[:,:nb]
            r1, d1, q1 = self.scratch1[:, :nb]
            r2, d2, q2 = self.scratch2[:, :nb]
            self.restrictor0(-local_residuals, r1, phases)
            d1 = 4 * step * r1  
            self.kin1.apply(d1, q1, phases)
            q1 -= r1
            self.restrictor1(q1, r2, phases)
            d2 = 16 * step * r2 
            self.kin2.apply(d2, q2, phases)
            q2 -= r2
            # d2 -= 16 * step * q2
            axpy(-16*step, q2, d2)
            self.interpolator2(d2, q1, phases)
            d1 -= q1
            self.kin1.apply(d1, q1, phases)
            q1 -= r1
            # d1 -= 4 * step * q1
            axpy(-4*step, q1, d1)
            self.interpolator1(-d1, d0, phases)
            self.kin0.apply(d0, q0, phases)
            q0 -= local_residuals
            axpy(-step, q0, d0)  # d0 -= step * q0
            scal(-1.0, d0)


        if self.cuda and not isinstance(residuals, gpuarray.GPUArray):
            return d0.get() 
        else:
            return d0


