# Copyright (C) 2003  CAMP
# Please see the accompanying LICENSE file for further information.

"""Grid transformers.

This module defines tools for doing interpolations/restrictions between
differentt uniform 3D grids.
"""

from __future__ import division
from math import pi
import numpy as np

from gpaw import debug
from gpaw.utilities import is_contiguous
import _gpaw

import pycuda.driver as cuda
import pycuda.gpuarray as gpuarray
from gpaw import debug_cuda,debug_cuda_reltol

class _Transformer:
    def __init__(self, gdin, gdout, nn=1, dtype=float, allocate=True, cuda=False):
        self.gdin = gdin
        self.gdout = gdout
        self.nn = nn
        assert 1 <= nn <= 4
        self.dtype = dtype

        pad_cd = np.empty((3, 2), int)
        neighborpad_cd = np.empty((3, 2), int)
        skip_cd = np.empty((3, 2), int)

        if (gdin.N_c == 2 * gdout.N_c).all():
            # Restriction:
            pad_cd[:, 0] = 2 * nn - 1 - 2 * gdout.beg_c + gdin.beg_c
            pad_cd[:, 1] = 2 * nn - 2 + 2 * gdout.end_c - gdin.end_c
            neighborpad_cd[:, 0] = 2 * nn - 2 + 2 * gdout.beg_c - gdin.beg_c
            neighborpad_cd[:, 1] = 2 * nn - 1 - 2 * gdout.end_c + gdin.end_c
            self.interpolate = False
        else:
            assert (gdout.N_c == 2 * gdin.N_c).all()
            # Interpolation:
            pad_cd[:, 0] = nn - 1 - gdout.beg_c // 2 + gdin.beg_c
            pad_cd[:, 1] = nn + gdout.end_c // 2 - gdin.end_c
            neighborpad_cd[:, 0] = nn + gdout.beg_c // 2 - gdin.beg_c
            neighborpad_cd[:, 1] = nn - 1 - gdout.end_c // 2 + gdin.end_c
            skip_cd[:, 0] = gdout.beg_c % 2
            skip_cd[:, 1] = gdout.end_c % 2
            self.interpolate = True
            assert ((gdin.n_c[0] + 2 * nn - 1) * (gdin.n_c[1] + 2 * nn - 1) <=
                    gdout.n_c[0] * gdout.n_c[1])

        assert np.alltrue(pad_cd.ravel() >= 0)
        self.ngpin = tuple(gdin.n_c)
        self.ngpout = tuple(gdout.n_c)
        assert dtype in [float, complex]
        self.transformer = None

        self.pad_cd = pad_cd
        self.neighborpad_cd = neighborpad_cd
        self.skip_cd = skip_cd

        self.cuda=cuda
        
        self.allocated = False
        if allocate:
            self.allocate()

    def allocate(self):
        assert not self.allocated
        gdin = self.gdin
        gdout = self.gdout
        if gdin.comm.size > 1:
            comm = gdin.comm.get_c_object()
        else:
            comm = None
        
        self.transformer = _gpaw.Transformer(gdin.n_c, gdout.n_c, 
                                             2 * self.nn, self.pad_cd, 
                                             self.neighborpad_cd, self.skip_cd,
                                             gdin.neighbor_cd,
                                             self.dtype == float, comm,
                                             self.interpolate, self.cuda)
        self.allocated = True
        
    def apply(self, input, output, phases=None):

        assert (type(input) == type(output))
        
        if isinstance(input,gpuarray.GPUArray) and  isinstance(output,gpuarray.GPUArray):
            #print "fd_transformer_apply_cuda_gpu"
            if debug_cuda:
                input_cpu = input.get()
                output_cpu = output.get()
                self.transformer.apply(input_cpu, output_cpu, phases)
            
            self.transformer.apply_cuda_gpu(input.gpudata, output.gpudata,
                                            input.shape, input.dtype,phases)
            if debug_cuda:
                diff = abs(output_cpu - output.get())
                error_i = np.unravel_index(np.argmax(diff), diff.shape)
                error = diff[error_i]
                if error > np.finfo(type(error)).eps and \
                       error > debug_cuda_reltol * abs(output_cpu[error_i]):
                    print "Debug cuda: transformer apply max error: ", error
                    #print input.shape,output.shape
                    #print input_cpu.shape,output_cpu.shape
                    #print phase_cd
                    #print output_cpu[1]-output.get()[1]
                    #assert 0
                    
        else:    
            self.transformer.apply(input, output, phases)

    def get_async_sizes(self):
        return self.transformer.get_async_sizes()

    def estimate_memory(self, mem):
        # Read transformers.c and bc.c for details
        # Notes: estimates are somewhat off, mostly around 100%-110%, but
        # below 100% for some small transformer objects.
        # How exactly it is possible for the estimate to be too high will
        # forever be a mystery.
        #
        # Accuracy not tested with OpenMP
        nbase_c = np.array(self.gdin.n_c)
        # Actual array size is gd array size plus pad_cd
        nactual_c = nbase_c + np.array(self.pad_cd).sum(1)
        # In the C code, the mysterious bc->ndouble is just 1 if real, else 2
        itemsize = mem.itemsize[self.dtype]
        inbytes = np.prod(nactual_c) * itemsize
        mem.subnode('buf', inbytes)
        if self.interpolate:
            mem.subnode('buf2 interp', 16 * inbytes)
        else:
            nactual_z = nactual_c[2] - 4 * self.nn + 3
            N = 1 * nactual_c[0] * nactual_c[1] * nactual_z // 2
            mem.subnode('buf2 restrict', N * itemsize)


class TransformerWrapper:
    def __init__(self, transformer):
        self.transformer = transformer
        self.allocated = transformer.allocated
        self.dtype = transformer.dtype
        self.ngpin = transformer.ngpin
        self.ngpout = transformer.ngpout

    def allocate(self):
        assert not self.allocated
        self.transformer.allocate()
        self.allocated = True

    def apply(self, input, output, phases=None):
        assert is_contiguous(input, self.dtype)
        assert is_contiguous(output, self.dtype)
        assert input.shape[-3:] == self.ngpin
        assert output.shape[-3:] == self.ngpout
        assert (self.dtype == float or
                (phases.dtype == complex and
                 phases.shape == (3, 2)))

        assert self.allocated
        self.transformer.apply(input, output, phases)

    def get_async_sizes(self):
        return self.transformer.get_async_sizes()

    def estimate_memory(self, mem):
        self.transformer.estimate_memory(mem)


def Transformer(gdin, gdout, nn=1, dtype=float, allocate=True,cuda=False):
    if nn != 9:
        t = _Transformer(gdin, gdout, nn, dtype, allocate,cuda)
        if debug:
            t = TransformerWrapper(t)
        return t
    class T:
        def allocate(self):
            pass
        def apply(self, input, output, phases=None):
            output[:] = input
        def estimate_memory(self, mem):
            mem.set('Nulltransformer', 0)
    return T()


def multiple_transform_apply(transformerlist, inputs, outputs, phases=None):
    return _gpaw.multiple_transform_apply(transformerlist, inputs, outputs, 
                                          phases)


def coefs(k, p):
    for i in range(0, k * p, p):
        print '%2d' % i,
        for x in range((k // 2 - 1) * p, k // 2 * p + 1):
            n = 1
            d = 1
            for j in range(0, k * p, p):
                if j == i:
                    continue
                n *= x - j
                d *= i - j
            print '%14.16f' % (n / d),
        print
