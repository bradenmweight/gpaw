from math import pi

import numpy as np
import pytest
from gpaw.core import UGDesc, PWDesc
from gpaw.fd_operators import Laplace
from gpaw.mpi import world


def test_real_to_complex_fft():
    L = 1.0
    n = 8
    a = UGDesc(cell=[L, L, L], size=(n, n, n), comm=world).empty()
    a.randomize()
    pw = PWDesc(ecut=50.0, cell=a.desc.cell, dtype=complex, comm=world)
    # real -> complex:
    with pytest.raises(TypeError):
        a.fft(pw=pw)
    # complex -> complex:
    b = a.to_complex().fft(pw=pw).gather(broadcast=True)
    # c(G) = c(-G)^*:
    coefs = {tuple(i_c): coef for i_c, coef in zip(pw.indices_cG.T, b.data)}
    for (i, j, k), coef in coefs.items():
        assert coef == coefs[(-i, -j, -k)].conj()


@pytest.mark.ci
def test_redist():
    a = 2.5
    n = 2
    grid1 = UGDesc(cell=[a, a, a], size=(n, n, n), comm=world)
    f1 = grid1.empty()
    f1.data[:] = world.rank + 1
    f2 = f1.gather()
    f3 = f1.gather(broadcast=True)
    if world.rank == 0:
        assert (f2.data == f3.data).all()
    print(f2)
    f4 = f1.new()
    f4.scatter_from(f2)
    assert (f4.data == f1.data).all()


def test_complex_laplace():
    a = 2.5
    grid = UGDesc(cell=[a, a, a],
                  size=(24, 8, 8),
                  kpt=[1 / 3, 0, 0],
                  comm=world)
    f = grid.empty()
    f.data[:] = 1.0
    f.multiply_by_eikr()
    lap = Laplace(grid._gd, n=2, dtype=complex)
    g = grid.empty()
    lap(f, g)
    k = 2 * pi / a / 3
    error = g.data.conj() * f.data + k**2
    assert abs(error).max() == pytest.approx(0.0, abs=1e-6)


def test_moment():
    L = 5
    n = 20
    grid = UGDesc(cell=[L, L, L], size=(n, n, n), comm=world)
    f = grid.zeros()

    # P-type Gaussian:
    l = 1
    a = 4.0
    rcut = 3.0
    p = (l, rcut, lambda r: np.exp(-a * r**2))

    if 0:  # Analytic result
        from sympy import integrate, exp, oo, var, Symbol, sqrt, pi
        x = var('x')
        a = Symbol('a', positive=True)
        m = 8 * (integrate(exp(-a * x**2) * x**2, (x, 0, oo)) *
                 integrate(exp(-a * x**2), (x, 0, oo))**2 *
                 sqrt(3 / (4 * pi)))
        print(m)  # sqrt(3)*pi/(4*a**(5/2))

    moment = 3**0.5 * np.pi / (4 * a**(5 / 2))

    # Add P-y function to grid:
    basis = grid.atom_centered_functions(
        [[p]],
        positions=[[0.5, 0.5, 0.5]])
    coefs = basis.layout.empty()
    if 0 in coefs:
        coefs[0] = [1.0, 0, 0]  # y, z, x
    basis.add_to(f, coefs)

    assert abs(f.integrate()) < 1e-14
    assert f.moment() == pytest.approx([0, moment, 0])

    pw = PWDesc(cell=grid.cell, ecut=700, comm=world)
    f2 = f.fft(pw=pw)

    assert abs(f2.integrate()) < 1e-14
    assert f2.moment() == pytest.approx([0, moment, 0], abs=1e-5)
