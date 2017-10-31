"""BLACS distributed matrix object."""
import numpy as np
import scipy.linalg as linalg

import _gpaw
from gpaw import debug
from gpaw.mpi import serial_comm
import gpaw.utilities.blas as blas


_global_blacs_context_store = {}


def matrix_matrix_multiply(alpha, a, opa, b, opb, beta=0.0, c=None,
                           symmetric=False):
    """BLAS-style matrix-matrix multiplication.

    Will use dgemm/zgemm/dsyrk/zherk/dsyr2k/zher2k as apropriate or the
    equivalent PBLAS functions for distributed matrices.

    The coefficients alpha and beta are of type float.  Matrices a, b and c
    must have same type (float or complex).  The strings apa and opb must be
    'N', 'T', or 'C' .  For opa='N' and opb='N', the operation performed is
    equivalent to::

        c.array[:] =  alpha * np.dot(a.array, b.array) + beta * c.array

    Replace a.array with a.array.T or a.array.T.conj() for opa='T' and 'C'
    resprctively (similarly for opb).

    Use symmetric=True if the result matrix is symmetric/hermetian
    (only lower half of c will be evaluated).
    """
    return _matrix(a).multiply(alpha, opa, _matrix(b), opb,
                               beta, c if c is None else _matrix(c),
                               symmetric)


def suggest_blocking(N, ncpus):
    """Suggest blocking of NxN matrix.

    Returns rows, columns, blocksize tuple."""

    nprow = ncpus
    npcol = 1

    # Get a sort of reasonable number of columns/rows
    while npcol < nprow and nprow % 2 == 0:
        npcol *= 2
        nprow //= 2

    assert npcol * nprow == ncpus

    # ScaLAPACK creates trouble if there aren't at least a few
    # whole blocks; choose block size so there will always be
    # several blocks.  This will crash for small test systems,
    # but so will ScaLAPACK in any case
    blocksize = min(-(-N // 4), 64)

    return nprow, npcol, blocksize


class Matrix:
    def __init__(self, M, N, dtype=None, data=None, dist=None):
        """Matrix object.

        M: int
            Rows.
        N: int
            Columns.
        dtype: type
            Data type (float or complex).
        dist: tuple or None
            BLACS distribution given as (communicator, rows, colums, blocksize)
            tuple.  Default is None meaning no distribution.
        data: ndarray or None.
            Numpy ndarray to use for starage.  By default, a new ndarray
            will be allocated.
            """
        self.shape = (M, N)

        if dtype is None:
            if data is None:
                dtype = float
            else:
                dtype = data.dtype
        self.dtype = np.dtype(dtype)

        dist = dist or ()
        if isinstance(dist, tuple):
            dist = create_distribution(M, N, *dist)
        self.dist = dist

        if data is None:
            self.array = np.empty(dist.shape, self.dtype)
        else:
            self.array = data.reshape(dist.shape)

        self.comm = serial_comm
        self.state = 'everything is fine'

    def __len__(self):
        return self.shape[0]

    def __repr__(self):
        dist = str(self.dist).split('(')[1]
        return 'Matrix({}: {}'.format(self.dtype.name, dist)

    def new(self, dist='inherit'):
        """Create new matrix of same shape and dtype.

        Default is to use same BLACS distribution.  Use dist to use another
        distribution.
        """
        return Matrix(*self.shape, dtype=self.dtype,
                      dist=self.dist if dist == 'inherit' else dist)

    def __setitem__(self, i, x):
        # assert i == slice(None)
        if isinstance(x, np.ndarray):
            1 / 0  # sssssself.array[:] = x
        else:
            x.eval(self)

    def __iadd__(self, x):
        x.eval(self, 1.0)
        return self

    def multiply(self, alpha, opa, b, opb, beta=0.0, out=None,
                 symmetric=False):
        """BLAS-style Matrix-matrix multiplication.

        See matrix_matrix_multipliction() for details.
        """
        dist = self.dist
        if out is None:
            assert beta == 0.0
            if opa == 'N':
                M = self.shape[0]
            else:
                M = self.shape[1]
            if opb == 'N':
                N = b.shape[1]
            else:
                N = b.shape[0]
            out = Matrix(M, N, self.dtype,
                         dist=(dist.comm, dist.rows, dist.columns))
        if alpha == 1.0 and beta == 0.0 and opa == 'N' and opb == 'N':
            if dist.comm.size > 1 and len(self) % dist.comm.size == 0:
                return fastmmm(self, b, out)

        dist.multiply(alpha, self, opa, b, opb, beta, out, symmetric)
        return out

    def redist(self, other):
        """Redistribute to other BLACS layout."""
        if self is other:
            return
        d1 = self.dist
        d2 = other.dist
        n1 = d1.rows * d1.columns
        n2 = d2.rows * d2.columns
        if n1 == n2 == 1:
            other.array[:] = self.array
            return
        c = d1.comm if d1.comm.size > d2.comm.size else d2.comm
        n = max(n1, n2)
        if n < c.size:
            c = c.new_communicator(np.arange(n))
        if c is not None:
            M, N = self.shape
            d1 = create_distribution(M, N, c,
                                     d1.rows, d1.columns, d1.blocksize)
            d2 = create_distribution(M, N, c,
                                     d2.rows, d2.columns, d2.blocksize)
            if n1 == n:
                ctx = d1.desc[1]
            else:
                ctx = d2.desc[1]
            redist(d1, self.array, d2, other.array, ctx)

    def invcholesky(self):
        """Inverse of Cholesky decomposition.

        Only the lower part is used.
        """
        if self.state == 'a sum is needed':
            self.comm.sum(self.array, 0)

        if self.comm.rank == 0:
            if self.dist.comm.size > 1:
                S = self.new(dist=(self.dist.comm, 1, 1))
                self.redist(S)
            else:
                S = self
            if self.dist.comm.rank == 0:
                if debug:
                    S.array[np.triu_indices(S.shape[0], 1)] = 42.0
                L_nn = linalg.cholesky(S.array,
                                       lower=True,
                                       overwrite_a=True,
                                       check_finite=debug)
                S.array[:] = linalg.inv(L_nn,
                                        overwrite_a=True,
                                        check_finite=debug)
            if S is not self:
                S.redist(self)

        if self.comm.size > 1:
            self.comm.broadcast(self.array, 0)
            self.state == 'everything is fine'

    def eigh(self, cc=False, scalapack=(None, 1, 1, None)):
        """Calculate eigenvectors and eigenvalues.

        Matrix must be symmetric/hermitian and stored in lower half.

        cc: bool
            Complex conjugate matrix before finding eigenvalues.
        scalapack: tuple
            BLACS distribution for ScaLapack to use.  Default is to do serial
            diagonalization.
        """
        slcomm, rows, columns, blocksize = scalapack

        if self.state == 'a sum is needed':
            self.comm.sum(self.array, 0)

        slcomm = slcomm or self.dist.comm
        dist = (slcomm, rows, columns, blocksize)

        redist = (rows != self.dist.rows or
                  columns != self.dist.columns or
                  blocksize != self.dist.blocksize)

        if redist:
            H = self.new(dist=dist)
            self.redist(H)
        else:
            assert self.dist.comm.size == slcomm.size
            H = self

        eps = np.empty(H.shape[0])

        if rows * columns == 1:
            if self.comm.rank == 0 and self.dist.comm.rank == 0:
                if cc and H.dtype == complex:
                    np.negative(H.array.imag, H.array.imag)
                if debug:
                    H.array[np.triu_indices(H.shape[0], 1)] = 42.0
                eps[:], H.array.T[:] = linalg.eigh(H.array,
                                                   lower=True,  # ???
                                                   overwrite_a=True,
                                                   check_finite=debug)
            self.dist.comm.broadcast(eps, 0)
        elif slcomm.rank < rows * columns:
            assert cc
            array = H.array.copy()
            info = _gpaw.scalapack_diagonalize_dc(array, H.dist.desc, 'U',
                                                  H.array, eps)
            assert info == 0, info

        if redist:
            H.redist(self)

        assert (self.state == 'a sum is needed') == (
            self.comm.size > 1)
        if self.comm.size > 1:
            self.comm.broadcast(self.array, 0)
            self.comm.broadcast(eps, 0)
            self.state == 'everything is fine'

        return eps

    def complex_conjugate(self):
        """Inplace complex conjugation."""
        if self.dtype == complex:
            np.negative(self.array.imag, self.array.imag)


def _matrix(M):
    """Dig out Matrix object from wrapper(s)."""
    if isinstance(M, Matrix):
        return M
    return _matrix(M.matrix)


class NoDistribution:
    comm = serial_comm
    rows = 1
    columns = 1
    blocksize = None

    def __init__(self, M, N):
        self.shape = (M, N)

    def __str__(self):
        return 'NoDistribution({}x{})'.format(*self.shape)

    def global_index(self, n):
        return n

    def multiply(self, alpha, a, opa, b, opb, beta, c, symmetric):
        if symmetric:
            assert opa == 'N'
            assert opb == 'C' or opb == 'T' and a.dtype == float
            if a is b:
                blas.rk(alpha, a.array, beta, c.array)
            else:
                if beta == 1.0 and a.shape[1] == 0:
                    return
                blas.r2k(0.5 * alpha, a.array, b.array, beta, c.array)
        else:
            blas.mmm(alpha, a.array, opa, b.array, opb, beta, c.array)


class BLACSDistribution:
    serial = False

    def __init__(self, M, N, comm, r, c, b):
        self.comm = comm
        self.rows = r
        self.columns = c
        self.blocksize = b

        key = (comm, r, c)
        context = _global_blacs_context_store.get(key)
        if context is None:
            context = _gpaw.new_blacs_context(comm.get_c_object(), c, r, 'R')
            _global_blacs_context_store[key] = context

        if b is None:
            if c == 1:
                br = (M + r - 1) // r
                bc = max(1, N)
            elif r == 1:
                br = M
                bc = (N + c - 1) // c
            else:
                raise ValueError('Please specify block size!')
        else:
            br = bc = b

        n, m = _gpaw.get_blacs_local_shape(context, N, M, bc, br, 0, 0)
        if n < 0 or m < 0:
            n = m = 0
        self.shape = (m, n)
        lld = max(1, n)
        self.desc = np.array([1, context, N, M, bc, br, 0, 0, lld], np.intc)

    def __str__(self):
        return ('BLACSDistribution(global={}, local={}, blocksize={})'
                .format(*('{}x{}'.format(*shape)
                          for shape in [self.desc[3:1:-1],
                                        self.shape,
                                        self.desc[5:3:-1]])))

    def global_index(self, myi):
        return self.comm.rank * int(self.desc[5]) + myi

    def multiply(self, alpha, a, opa, b, opb, beta, c, symmetric):
        if symmetric:
            assert opa == 'N'
            assert opb == 'C' or opb == 'T' and a.dtype == float
            N, K = a.shape
            if a is b:
                _gpaw.pblas_rk(N, K, alpha, a.array,
                               beta, c.array,
                               a.dist.desc, c.dist.desc,
                               'U')
            else:
                _gpaw.pblas_r2k(N, K, 0.5 * alpha, b.array, a.array,
                                beta, c.array,
                                b.dist.desc, a.dist.desc, c.dist.desc,
                                'U')
        else:
            Ka, M = a.shape
            N, Kb = b.shape
            if opa == 'N':
                Ka, M = M, Ka
            if opb == 'N':
                N, Kb = Kb, N
            _gpaw.pblas_gemm(N, M, Ka, alpha, b.array, a.array,
                             beta, c.array,
                             b.dist.desc, a.dist.desc, c.dist.desc,
                             opb, opa)


def redist(dist1, M1, dist2, M2, context):
    _gpaw.scalapack_redist(dist1.desc, dist2.desc,
                           M1, M2,
                           dist1.desc[2], dist1.desc[3],
                           1, 1, 1, 1,  # 1-indexing
                           context, 'G')


def create_distribution(M, N, comm=None, r=1, c=1, b=None):
    if comm is None or comm.size == 1:
        assert r == 1 and abs(c) == 1 or c == 1 and abs(r) == 1
        return NoDistribution(M, N)

    return BLACSDistribution(M, N, comm,
                             r if r != -1 else comm.size,
                             c if c != -1 else comm.size,
                             b)


def fastmmm(m1, m2, m3):
    comm = m1.dist.comm

    n = len(m1.array)
    buf1 = m2.array
    buf2 = np.empty_like(buf1)

    beta = 0.0

    for r in range(comm.size - 1):
        rrank = (comm.rank + r + 1) % comm.size
        srank = (comm.rank - r - 1) % comm.size
        rrequest = comm.receive(buf2, rrank, 21, False)
        srequest = comm.send(m2.array, srank, 21, False)

        n1 = (comm.rank + r) % comm.size * n
        n2 = n1 + n
        blas.mmm(1.0, m1.array[:, n1:n2], 'N', buf1, 'N', beta, m3.array)

        beta = 1.0

        if r == 0:
            buf1 = np.empty_like(buf2)

        buf1, buf2 = buf2, buf1

        comm.wait(rrequest)
        comm.wait(srequest)

    n1 = rrank * n
    n2 = n1 + n
    blas.mmm(1.0, m1.array[:, n1:n2], 'N', buf1, 'N', beta, m3.array)

    return m3
