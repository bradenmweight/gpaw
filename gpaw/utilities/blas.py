# Copyright (C) 2003  CAMP
# Please see the accompanying LICENSE file for further information.

"""
Python wrapper functions for the ``C`` package:
Basic Linear Algebra Subroutines (BLAS)

See also:
http://en.wikipedia.org/wiki/Basic_Linear_Algebra_Subprograms
and
http://www.netlib.org/lapack/lug/node145.html
"""

import numpy as np

from gpaw.utilities import is_contiguous
from gpaw import debug

import _gpaw

from gpaw import debug_cuda,debug_cuda_tol
import pycuda.driver as cuda
import pycuda.gpuarray as gpuarray

def scal(alpha, x):
    """alpha x

    Performs the operation::

      x <- alpha * x 
      
    """
    if isinstance(alpha, complex):
        assert is_contiguous(x, complex)
    else:
        assert isinstance(alpha, float)
        assert x.dtype in [float, complex]
        assert x.flags.contiguous
    _gpaw.scal(alpha, x)
    

def gemm(alpha, a, b, beta, c, transa='n', cuda=False):
    """General Matrix Multiply.

    Performs the operation::
    
      c <- alpha * b.a + beta * c

    If transa is "n", ``b.a`` denotes the matrix multiplication defined by::
    
                      _
                     \  
      (b.a)        =  ) b  * a
           ijkl...   /_  ip   pjkl...
                      p
    
    If transa is "t" or "c", ``b.a`` denotes the matrix multiplication
    defined by::
    
                      _
                     \  
      (b.a)        =  ) b    *    a
           ij        /_  iklm...   jklm...
                     klm... 

    where in case of "c" also complex conjugate of a is taken.
    """

    if debug:
        assert np.isfinite(c).all()
    
        assert (a.dtype == float and b.dtype == float and c.dtype == float and
                isinstance(alpha, float) and isinstance(beta, float) or
                a.dtype == complex and b.dtype == complex and c.dtype == complex)
        assert a.flags.contiguous
        if transa == 'n':
            assert c.flags.contiguous or c.ndim == 2 and c.strides[1] == c.itemsize
            assert b.ndim == 2
            assert b.strides[1] == b.itemsize
            assert a.shape[0] == b.shape[1]
            assert c.shape == b.shape[0:1] + a.shape[1:]
        else:
            assert b.flags.contiguous
            assert c.strides[1] == c.itemsize
            assert a.shape[1:] == b.shape[1:]
            assert c.shape == (b.shape[0], a.shape[0])
            
    assert (type(a) == type(b) and type(b) == type(c))


    if isinstance(a,gpuarray.GPUArray) and  isinstance(b,gpuarray.GPUArray) and isinstance(c,gpuarray.GPUArray):
        #print "gemm_cuda_gpu",a.shape,b.shape,c.shape,a.dtype,transa
        if debug_cuda:
            a_cpu=a.get()
            b_cpu=b.get()
            c_cpu=c.get()
            _gpaw.gemm(alpha, a_cpu, b_cpu, beta, c_cpu, transa)
        _gpaw.gemm_cuda_gpu(alpha, a.gpudata,a.shape, b.gpudata, 
                            b.shape,beta, c.gpudata, c.shape, a.dtype, transa)
        if debug_cuda:
            error=np.max(abs(c_cpu-c.get()))
            if error>debug_cuda_tol:
                print "Debug cuda: gemm max error: ", error
    elif cuda:
        print "gemm_cuda_cpu"
        a_gpu = gpuarray.to_gpu(a)
        b_gpu = gpuarray.to_gpu(b)
        c_gpu = gpuarray.to_gpu(c)
        _gpaw.gemm_cuda_gpu(alpha, a_gpu.gpudata, a.shape, b_gpu.gpudata, 
                            b.shape,beta, c_gpu.gpudata, c.shape, a.dtype, transa)
        c_gpu.get(c)
        #_gpaw.gemm(alpha, a, b, beta, c, transa)
    else:
        _gpaw.gemm(alpha, a, b, beta, c, transa)

def gemv(alpha, a, x, beta, y, trans='t'):
    """General Matrix Vector product.

    Performs the operation::

      y <- alpha * a.x + beta * y

    ``a.x`` denotes matrix multiplication, where the product-sum is
    over the entire length of the vector x and 
    the first dimension of a (for trans='n'), or
    the last dimension of a (for trans='t' or 'c').

    If trans='c', the complex conjugate of a is used. The default is
    trans='t', i.e. behaviour like np.dot with a 2D matrix and a vector.

    Example::

      >>> y_m = np.dot(A_mn, x_n)
      >>> # or better yet
      >>> y_m = np.zeros(A_mn.shape[0], A_mn.dtype)
      >>> gemv(1.0, A_mn, x_n, 0.0, y_m)

    """
    assert (a.dtype == float and x.dtype == float and y.dtype == float and
            isinstance(alpha, float) and isinstance(beta, float) or
            a.dtype == complex and x.dtype == complex and y.dtype == complex)
    assert a.flags.contiguous
    assert y.flags.contiguous
    assert x.ndim == 1
    assert y.ndim == a.ndim-1
    if trans == 'n':
        assert a.shape[0] == x.shape[0]
        assert a.shape[1:] == y.shape
    else:
        assert a.shape[-1] == x.shape[0]
        assert a.shape[:-1] == y.shape
    _gpaw.gemv(alpha, a, x, beta, y, trans)


def axpy(alpha, x, y, cuda=False):
    """alpha x plus y.

    Performs the operation::

      y <- alpha * x + y
      
    """

    if debug:
        if isinstance(alpha, complex):
            assert is_contiguous(x, complex) and is_contiguous(y, complex)
        else:
            assert isinstance(alpha, float)
            assert x.dtype in [float, complex]
            assert x.dtype == y.dtype
            assert x.flags.contiguous and y.flags.contiguous
        assert x.shape == y.shape
    
    assert type(x) == type(y)
    
    if isinstance(x,gpuarray.GPUArray) and isinstance(y,gpuarray.GPUArray):
        #print "axpy_gpu"
        if debug_cuda:
            x_cpu=x.get()
            y_cpu=y.get()
            _gpaw.axpy(alpha, x_cpu, y_cpu)

        _gpaw.axpy_cuda_gpu(alpha, x.gpudata, x.shape, y.gpudata, y.shape, x.dtype)
        if debug_cuda:
            error=np.max(abs(y_cpu-y.get()))
            if error>debug_cuda_tol:
                print "Debug cuda: axpy max error: ", error
    elif cuda:
        print "axpy_cuda_cpu"
        x_gpu=gpuarray.to_gpu(x)
        y_gpu=gpuarray.to_gpu(y)

        _gpaw.axpy_cuda_gpu(alpha, x_gpu.gpudata, x.shape, y_gpu.gpudata, 
                            y.shape, x.dtype)
        y_gpu.get(y)
    else:
        _gpaw.axpy(alpha, x, y)

def rk(alpha, a, beta, c, cuda=False):
    """Rank-k update of a matrix.

    Performs the operation::
    
                        dag
      c <- alpha * a . a    + beta * c

    where ``a.b`` denotes the matrix multiplication defined by::

                 _
                \  
      (a.b)   =  ) a         * b
           ij   /_  ipklm...     pjklm...
               pklm...

    ``dag`` denotes the hermitian conjugate (complex conjugation plus a
    swap of axis 0 and 1).
    
    Only the lower triangle of ``c`` will contain sensible numbers.
    """

    if debug:
        assert np.isfinite(c).all()

        assert (a.dtype == float and c.dtype == float or
                a.dtype == complex and c.dtype == complex)
        assert a.flags.contiguous
        assert a.ndim > 1
        assert c.shape == (a.shape[0], a.shape[0])
        assert c.strides[1] == c.itemsize

    assert type(a) == type(c)
    
    if isinstance(a,gpuarray.GPUArray) and isinstance(c,gpuarray.GPUArray):
        # print "rk_gpu"
        if debug_cuda:
            a_cpu=a.get()
            c_cpu=c.get()
            _gpaw.rk(alpha, a_cpu, beta, c_cpu)

        _gpaw.rk_cuda_gpu(alpha, a.gpudata, a.shape,beta, c.gpudata, c.shape, a.dtype)
        if debug_cuda:
            error=np.max(abs(c_cpu-c.get()))
            if error>debug_cuda_tol:
                print "Debug cuda: rk max error: ", error
    elif cuda:
        # print "rk_cuda_cpu"
        a_gpu=gpuarray.to_gpu(a)
        c_gpu=gpuarray.to_gpu(c)

        _gpaw.rk_cuda_gpu(alpha, a_gpu.gpudata, a.shape,beta, 
                          c_gpu.gpudata, c.shape, a.dtype)
        c_gpu.get(c)
        #_gpaw.rk(alpha, a, beta, c)
    else:
        _gpaw.rk(alpha, a, beta, c)

def r2k(alpha, a, b, beta, c, cuda=False):
    """Rank-2k update of a matrix.

    Performs the operation::

                        dag        cc       dag
      c <- alpha * a . b    + alpha  * b . a    + beta * c

    where ``a.b`` denotes the matrix multiplication defined by::

                 _
                \ 
      (a.b)   =  ) a         * b
           ij   /_  ipklm...     pjklm...
               pklm...

    ``cc`` denotes complex conjugation.
    
    ``dag`` denotes the hermitian conjugate (complex conjugation plus a
    swap of axis 0 and 1).

    Only the lower triangle of ``c`` will contain sensible numbers.
    """

    if debug:
        assert np.isfinite(c).all()
        
        assert (a.dtype == float and b.dtype == float and c.dtype == float or
                a.dtype == complex and b.dtype == complex and c.dtype == complex)
        assert a.flags.contiguous and b.flags.contiguous
        assert np.rank(a) > 1
        assert a.shape == b.shape
        assert c.shape == (a.shape[0], a.shape[0])
        assert c.strides[1] == c.itemsize


    assert (type(a) == type(b) and type(b) == type(c))


    if isinstance(a,gpuarray.GPUArray) and  isinstance(b,gpuarray.GPUArray) and isinstance(c,gpuarray.GPUArray):
        # print "r2k_gpu"
        if debug_cuda:
            a_cpu=a.get()
            b_cpu=b.get()
            c_cpu=c.get()            
            _gpaw.r2k(alpha, a_cpu, b_cpu, beta, c_cpu)
            
        _gpaw.r2k_cuda_gpu(alpha, a.gpudata, a.shape, b.gpudata,
                           b.shape, beta, c.gpudata, c.shape, a.dtype)
        if debug_cuda:
            error=np.max(abs(c_cpu-c.get()))
            if error>debug_cuda_tol:
                print "Debug cuda: rk2 max error: ", error
    elif cuda:
        print "r2k_cuda_cpu"
        a_gpu=gpuarray.to_gpu(a)
        b_gpu=gpuarray.to_gpu(b)
        c_gpu=gpuarray.to_gpu(c)

        _gpaw.r2k_cuda_gpu(alpha, a_gpu.gpudata, a.shape,b_gpu.gpudata,
                           b.shape, beta, c_gpu.gpudata, c.shape, a.dtype)
        bc_gpu.get(c)
    else:
        _gpaw.r2k(alpha, a, b, beta, c)


def dotc(a, b):
    """Dot product, conjugating the first vector with complex arguments.

    Returns the value of the operation::

        _
       \   cc   
        ) a       * b
       /_  ijk...    ijk...       
       ijk...

    ``cc`` denotes complex conjugation.
    """

    if debug:
        assert ((is_contiguous(a, float) and is_contiguous(b, float)) or
                (is_contiguous(a, complex) and is_contiguous(b,complex)))
        assert a.shape == b.shape


    assert (type(a) == type(b))
    if isinstance(a,gpuarray.GPUArray) and  isinstance(b,gpuarray.GPUArray):
        if debug_cuda:
            gpu=_gpaw.dotc_cuda_gpu(a.gpudata, a.shape, b.gpudata, a.dtype)
            a_cpu=a.get()
            b_cpu=b.get()
            cpu=_gpaw.dotc(a_cpu, b_cpu)
            error=abs(gpu-cpu)
            if error>debug_cuda_tol:
                print "Debug cuda: dotc max error: ", error
            return gpu
        else:
            return _gpaw.dotc_cuda_gpu(a.gpudata, a.shape, b.gpudata, a.dtype)
    else:
        return _gpaw.dotc(a, b)

def dotu(a, b):
    """Dot product, NOT conjugating the first vector with complex arguments.

    Returns the value of the operation::

        _
       \ 
        ) a       * b
       /_  ijk...    ijk...
       ijk...


    """

    if debug:
        assert ((is_contiguous(a, float) and is_contiguous(b, float)) or
                (is_contiguous(a, complex) and is_contiguous(b,complex)))
        assert a.shape == b.shape

    assert (type(a) == type(b))
    if isinstance(a,gpuarray.GPUArray) and  isinstance(b,gpuarray.GPUArray):
        if debug_cuda:
            a_cpu=a.get()
            b_cpu=b.get()
            cpu=_gpaw.dotu(a_cpu, b_cpu)
            gpu=_gpaw.dotu_cuda_gpu(a.gpudata, a.shape, b.gpudata, a.dtype)
            error=abs(gpu-cpu)
            if error>debug_cuda_tol:
                print "Debug cuda: dotu max error: ", error
            return gpu
        else:
            return _gpaw.dotu_cuda_gpu(a.gpudata, a.shape, b.gpudata, a.dtype)
    else:
        return _gpaw.dotu(a, b)

def _gemmdot(a, b, alpha=1.0, beta=1.0, out=None, trans='n'):
    """Matrix multiplication using gemm.

    return reference to out, where::

      out <- alpha * a . b + beta * out

    If out is None, a suitably sized zero array will be created.

    ``a.b`` denotes matrix multiplication, where the product-sum is
    over the last dimension of a, and either
    the first dimension of b (for trans='n'), or
    the last dimension of b (for trans='t' or 'c').

    If trans='c', the complex conjugate of b is used.
    """
    # Store original shapes
    ashape = a.shape
    bshape = b.shape

    # Vector-vector multiplication is handled by dotu
    if a.ndim == 1 and b.ndim == 1:
        assert out is None
        if trans is 'c':
            return beta * _gpaw.dotc(b, a) # dotc conjugates *first* argument
        else:
            return beta * _gpaw.dotu(a, b)

##     # Use gemv if a or b is a vector, and the other is a matrix??
##     if a.ndim == 1 and trans == 'n':
##         gemv(alpha, b, a, beta, out, trans='n')
##     if b.ndim == 1 and trans == 'n':
##         gemv(alpha, a, b, beta, out, trans='t')

    # Map all arrays to 2D arrays
    a = a.reshape(-1, a.shape[-1])
    if trans == 'n':
        b = b.reshape(b.shape[0], -1)
        outshape = a.shape[0], b.shape[1]
    else: # 't' or 'c'
        b = b.reshape(-1, b.shape[-1])
    
    # Apply BLAS gemm routine
    outshape = a.shape[0], b.shape[trans == 'n']
    if out is None:
        # (ATLAS can't handle uninitialized output array)
        out = np.zeros(outshape, a.dtype)
    else:
        out = out.reshape(outshape)
    gemm(alpha, b, a, beta, out, trans)

    # Determine actual shape of result array
    if trans == 'n':
        outshape = ashape[:-1] + bshape[1:]
    else: # 't' or 'c'
        outshape = ashape[:-1] + bshape[:-1]
    return out.reshape(outshape)


def _rotate(in_jj, U_ij, a=1., b=0., out_ii=None, work_ij=None):
    """Perform matrix rotation using gemm

    For the 2D input matrices in, U, do the rotation::

      out <- a * U . in . U^d + b * out

    where '.' denotes matrix multiplication and '^d' the hermitian conjugate.

    work_ij is a temporary work array for storing the intermediate product.
    out_ii, and work_ij are created if not given.

    The method returns a reference to out.
    """
    if work_ij is None:
        work_ij = np.zeros_like(U_ij)
    if out_ii is None:
        out_ii = np.zeros(U_ij.shape[:1] * 2, U_ij.dtype)
    if in_jj.dtype == float:
        trans = 't'
    else:
        trans = 'c'
    gemm(1., in_jj, U_ij, 0., work_ij, 'n')
    gemm(a, U_ij, work_ij, b, out_ii, trans)
    return out_ii

if not debug:
    scal = _gpaw.scal
    #gemm = _gpaw.gemm
    gemv = _gpaw.gemv
    #axpy = _gpaw.axpy
    #rk = _gpaw.rk
    #r2k = _gpaw.r2k
    #dotc = _gpaw.dotc
    #dotu = _gpaw.dotu
    gemmdot = _gemmdot
    rotate = _rotate
else:
    def gemmdot(a, b, alpha=1.0, beta=1.0, out=None, trans='n'):
        assert a.flags.contiguous
        assert b.flags.contiguous
        assert a.dtype == b.dtype
        if trans == 'n':
            assert a.shape[-1] == b.shape[0]
        else:
            assert a.shape[-1] == b.shape[-1]
        if out is not None:
            assert out.flags.contiguous
            assert a.dtype == out.dtype
            assert a.ndim > 1 or b.ndim > 1
            if trans == 'n':
                assert out.shape == a.shape[:-1] + b.shape[1:]
            else:
                assert out.shape == a.shape[:-1] + b.shape[:-1]
        return _gemmdot(a, b, alpha, beta, out, trans)

    def rotate(in_jj, U_ij, a=1., b=0., out_ii=None, work_ij=None):
        assert in_jj.dtype == U_ij.dtype
        assert in_jj.flags.contiguous
        assert U_ij.flags.contiguous
        assert in_jj.shape == U_ij.shape[1:] * 2
        if out_ii is not None:
            assert out_ii.dtype == in_jj.dtype
            assert out_ii.flags.contiguous
            assert out_ii.shape == U_ij.shape[:1] * 2
        if work_ij is not None:
            assert work_ij.dtype == in_jj.dtype
            assert work_ij.flags.contiguous
            assert work_ij.shape == U_ij.shape
        return _rotate(in_jj, U_ij, a, b, out_ii, work_ij)
