import pytest
from gpaw.setup import create_setup


@pytest.mark.parametrize('symbol, M, n',
                         [('Pd', 1, 6),  # 4d -> 5s
                          ('Fe', 3, 5),
                          ('V', 3, 5),
                          ('Ti', 2, 5)])
def test_initial_occs(symbol, M, n):
    s = create_setup(symbol)
    f_si = s.calculate_initial_occupation_numbers(magmom=M,
                                                  hund=False,
                                                  charge=0,
                                                  nspins=2)
    print(f_si)
    magmom = (f_si[0] - f_si[1]).sum()
    assert abs(magmom - M) < 1e-10
    N = ((f_si[0] - f_si[1]) != 0).sum()
    assert n == N, 'Wrong # of values have changed'
