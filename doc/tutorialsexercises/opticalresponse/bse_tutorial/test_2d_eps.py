"""Test the calculated 2D dielectric function."""

import numpy as np
import pytest

from gpaw.test import findpeak


def test():
    """Test data in 2d_eps.png figure"""
    # Load data
    dat = np.loadtxt('2d_eps.dat')
    qid_q = dat[:, 0]
    eps_bare_q = dat[:, 1]
    eps_trunc_q = dat[:, 2]

    # Test the long wavelength limit
    assert eps_bare_q[0] == pytest.approx(6.547, abs=0.01)
    assert eps_trunc_q[0] == pytest.approx(1.)

    # Test maxima
    qmax_bare, epsmax_bare = findpeak(qid_q, eps_bare_q)
    qmax_trunc, epsmax_trunc = findpeak(qid_q, eps_trunc_q)
    assert qmax_bare == pytest.approx(4.563, abs=0.01)
    assert epsmax_bare == pytest.approx(9.206, abs=0.01)
    assert qmax_trunc == pytest.approx(4.947, abs=0.01)
    assert epsmax_trunc == pytest.approx(9.067, abs=0.01)

    # Test that the two kernels yield identical results at short wave lengths
    assert eps_trunc_q[10:] == pytest.approx(eps_bare_q[10:], abs=0.01)

    # Test values in the tail
    assert eps_bare_q[10::5] == pytest.approx([5.861, 3.143, 2.041], abs=0.01)


if __name__ == '__main__':
    test()
