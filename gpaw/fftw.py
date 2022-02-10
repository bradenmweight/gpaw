"""Python wrapper for FFTW3 library."""

import numpy as np

import _gpaw
from gpaw.typing import Array3D

ESTIMATE = 64
MEASURE = 0
PATIENT = 32
EXHAUSTIVE = 8


def have_fftw() -> bool:
    return hasattr(_gpaw, 'FFTWPlan')


def check_fft_size(n: int, factors=[2, 3, 5, 7]) -> bool:
    """Check if n is an efficient fft size.

    Efficient means that n can be factored into small primes (2, 3, 5, 7).

    >>> check_fft_size(17)
    False
    >>> check_fft_size(18)
    True
    """

    if n == 1:
        return True
    for x in factors:
        if n % x == 0:
            return check_fft_size(n // x, factors)
    return False


def get_efficient_fft_size(N: int, n=1, factors=[2, 3, 5, 7]) -> int:
    """Return smallest efficient fft size.

    Must be greater than or equal to N and divisible by n.

    >>> get_efficient_fft_size(17)
    18
    """
    N = -(-N // n) * n
    while not check_fft_size(N, factors):
        N += n
    return N


def check_fftw_inputs(in_R, out_R):
    for arr in in_R, out_R:
        # Note: Arrays not necessarily contiguous due to 16-byte alignment
        assert arr.ndim == 3  # We can perhaps relax this requirement
        assert arr.dtype == float or arr.dtype == complex

    if in_R.dtype == out_R.dtype == complex:
        assert in_R.shape == out_R.shape
    else:
        # One real and one complex:
        R, C = (in_R, out_R) if in_R.dtype == float else (out_R, in_R)
        assert C.dtype == complex
        assert R.shape[:2] == C.shape[:2]
        assert C.shape[2] == 1 + R.shape[2] // 2


class FFTPlan:
    """FFT 3d transform."""
    def __init__(self,
                 in_R: Array3D,
                 out_R: Array3D,
                 sign: int,
                 flags: int = MEASURE):
        check_fftw_inputs(in_R, out_R)
        self.in_R = in_R
        self.out_R = out_R
        self.sign = sign
        self.flags = flags

    def execute(self) -> None:
        raise NotImplementedError


class FFTWPlan(FFTPlan):
    """FFTW3 3d transform."""
    def __init__(self, in_R, out_R, sign, flags=MEASURE):
        if not have_fftw():
            raise ImportError('Not compiled with FFTW.')
        self._ptr = _gpaw.FFTWPlan(in_R, out_R, sign, flags)
        FFTPlan.__init__(self, in_R, out_R, sign, flags)

    def execute(self):
        _gpaw.FFTWExecute(self._ptr)

    def __del__(self):
        if getattr(self, '_ptr', None) and _gpaw is not None:
            _gpaw.FFTWDestroy(self._ptr)
        self._ptr = None


class NumpyFFTPlan(FFTPlan):
    """Numpy fallback."""
    def execute(self):
        if self.in_R.dtype == float:
            self.out_R[:] = np.fft.rfftn(self.in_R)
        elif self.out_R.dtype == float:
            self.out_R[:] = np.fft.irfftn(self.in_R, self.out_R.shape)
            self.out_R *= self.out_R.size
        elif self.sign == 1:
            self.out_R[:] = np.fft.ifftn(self.in_R, self.out_R.shape)
            self.out_R *= self.out_R.size
        else:
            self.out_R[:] = np.fft.fftn(self.in_R)


def empty(shape, dtype=float):
    """numpy.empty() equivalent with 16 byte alignment."""
    assert dtype == complex
    N = np.prod(shape)
    a = np.empty(2 * N + 1)
    offset = (a.ctypes.data % 16) // 8
    a = a[offset:2 * N + offset].view(complex)
    a.shape = shape
    return a


def create_plan(in_R: Array3D,
                out_R: Array3D,
                sign: int,
                flags: int = MEASURE) -> FFTPlan:
    if have_fftw():
        return FFTWPlan(in_R, out_R, sign, flags)
    return NumpyFFTPlan(in_R, out_R, sign, flags)
