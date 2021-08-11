r"""
Real-valued spherical harmonics


=== === === =======
 L   l   m
=== === === =======
 0   0   0   1
 1   1  -1   y
 2   1   0   z
 3   1   1   x
 4   2  -2   xy
 5   2  -1   yz
 6   2   0   3z2-r2
 7   2   1   zx
 8   2   2   x2-y2
=== === === =======

For a more complete list, see c/bmgs/sharmonic.py


Gaunt coefficients::

                  __
     ^      ^    \   L      ^
  Y (r)  Y (r) =  ) G    Y (r)
   L      L      /__ L L  L
    1      2      L   1 2

"""

import numpy as np

from math import pi
from collections import defaultdict
from _gpaw import spherical_harmonics as Yl

__all__ = ['Y', 'YL', 'nablarlYL', 'Yl']

names = [['1'],
         ['y', 'z', 'x'],
         ['xy', 'yz', '3z2-r2', 'zx', 'x2-y2'],
         ['3x2y-y3', 'xyz', '4yz2-y3-x2y', '2z3-3x2z-3y2z', '4xz2-x3-xy2',
          'x2z-y2z', 'x3-3xy2']]


def Y(L, x, y, z):
    result = 0.0
    for c, n in YL[L]:
        result += c * x**n[0] * y**n[1] * z**n[2]
    return result


def Yarr(L_M, R_Av):
    """
    Calculate spherical harmonics L_M at positions R_Av, where
    A is some array like index.
    """
    Y_MA = np.zeros((len(L_M), *R_Av.shape[:-1]))
    for M, L in enumerate(L_M):
        for c, n in YL[L]:  # could be vectorized further
            Y_MA[M] += c * np.prod(np.power(R_Av, n), axis=-1)
    return Y_MA


def nablarlYL(L, R):
    """Calculate the gradient of a real solid spherical harmonic."""
    x, y, z = R
    dYdx = dYdy = dYdz = 0.0
    terms = YL[L]
    # The 'abs' avoids error in case powx == 0
    for N, (powx, powy, powz) in terms:
        dYdx += N * powx * x**abs(powx - 1) * y**powy * z**powz
        dYdy += N * powy * x**powx * y**abs(powy - 1) * z**powz
        dYdz += N * powz * x**powx * y**powy * z**abs(powz - 1)
    return dYdx, dYdy, dYdz


def YYY(l, m):
    from fractions import Fraction
    from sympy import assoc_legendre, sqrt, simplify, factorial as fac, I, pi
    from sympy.abc import x, y, z
    c = sqrt((2 * l + 1) * fac(l - m) / fac(l + m) / 4 / pi)
    if m > 0:
        return simplify(c * (x + I * y)**m / (1 - z**2)**Fraction(m, 2) *
                        assoc_legendre(l, m, z))
    return simplify(c * (x - I * y)**(-m) / (1 - z**2)**Fraction(-m, 2) *
                    assoc_legendre(l, m, z))


def YYYY(l, m):
    from sympy import I, Number, sqrt
    if m > 0:
        return (YYY(l, m) + (-1)**m * YYY(l, -m)) / sqrt(2) * (-1)**m
    if m < 0:
        return -(YYY(l, m) - Number(-1)**m * YYY(l, -m)) / (sqrt(2) * I)
    return YYY(l, m)


def f(l, m):
    from sympy import Poly
    from sympy.abc import x, y, z
    Y = YYYY(l, m)
    coeffs = {}
    for nx, coef in enumerate(reversed(Poly(Y, x).all_coeffs())):
        for ny, coef in enumerate(reversed(Poly(coef, y).all_coeffs())):
            for nz, coef in enumerate(reversed(Poly(coef, z).all_coeffs())):
                if coef:
                    coeffs[(nx, ny, nz)] = coef
    return coeffs


def fix(coeffs, l):
    from sympy import Number
    new = defaultdict(lambda: Number(0))
    for (nx, ny, nz), coef in coeffs.items():
        if nx + ny + nz == l:
            new[(nx, ny, nz)] += coef
        else:
            new[(nx + 2, ny, nz)] += coef
            new[(nx, ny + 2, nz)] += coef
            new[(nx, ny, nz + 2)] += coef

    new = {nxyz: coef for nxyz, coef in new.items() if coef}

    if not all(sum(nxyz) == l for nxyz in new):
        new = fix(new, l)

    return new


def print_YL_table_code():
    print('# Computer generated table - do not touch!')
    print('YL = [')
    print('    # s, l=0:')
    print(f'    [({(4 * pi)**-0.5}, (0, 0, 0))],')
    for l in range(1, 8):
        s = 'spdfghijk'[l]
        print(f'    # {s}, l={l}:')
        for m in range(-l, l + 1):
            y1 = f(l, m)
            e = fix(y1, l)
            if l**2 + m + l < len(YL):
                assert len(e) == len(YL[l**2 + m + l])
                for c0, n in YL[l**2 + m + l]:
                    c = e[n].evalf()
                    assert abs(c0 - c) < 1e-10
            terms = []
            for n, en in e.items():
                c = en.evalf()
                terms.append(f'({c!r}, {n})')
            print('    [' + ',\n     '.join(terms) + '],')
    print(']')


def write_c_code(l: int) -> None:
    print(f'          else if (l == {l})')
    print('            {')
    for m in range(2 * l + 1):
        terms = []
        for c, n in YL[l**2 + m]:
            terms.append(f'{c!r} * ' + '*'.join('x' * n[0] +
                                                'y' * n[1] +
                                                'z' * n[2]))
        print(f'              Y_m[{m}] = ' + ' + '.join(terms) + ';')
    print('            }')


# Computer generated table - do not touch!
YL = [
    # s, l=0:
    [(0.28209479177387814, (0, 0, 0))],
    # p, l=1:
    [(0.488602511902920, (0, 1, 0))],
    [(0.488602511902920, (0, 0, 1))],
    [(0.488602511902920, (1, 0, 0))],
    # d, l=2:
    [(1.09254843059208, (1, 1, 0))],
    [(1.09254843059208, (0, 1, 1))],
    [(-0.315391565252520, (2, 0, 0)),
     (-0.315391565252520, (0, 2, 0)),
     (0.630783130505040, (0, 0, 2))],
    [(1.09254843059208, (1, 0, 1))],
    [(-0.546274215296040, (0, 2, 0)),
     (0.546274215296040, (2, 0, 0))],
    # f, l=3:
    [(-0.590043589926644, (0, 3, 0)),
     (1.77013076977993, (2, 1, 0))],
    [(2.89061144264055, (1, 1, 1))],
    [(-0.457045799464466, (2, 1, 0)),
     (-0.457045799464466, (0, 3, 0)),
     (1.82818319785786, (0, 1, 2))],
    [(-1.11952899777035, (2, 0, 1)),
     (-1.11952899777035, (0, 2, 1)),
     (0.746352665180231, (0, 0, 3))],
    [(-0.457045799464466, (3, 0, 0)),
     (-0.457045799464466, (1, 2, 0)),
     (1.82818319785786, (1, 0, 2))],
    [(-1.44530572132028, (0, 2, 1)),
     (1.44530572132028, (2, 0, 1))],
    [(-1.77013076977993, (1, 2, 0)),
     (0.590043589926644, (3, 0, 0))],
    # g, l=4:
    [(-2.50334294179670, (1, 3, 0)),
     (2.50334294179670, (3, 1, 0))],
    [(-1.77013076977993, (0, 3, 1)),
     (5.31039230933979, (2, 1, 1))],
    [(-0.946174695757560, (3, 1, 0)),
     (-0.946174695757560, (1, 3, 0)),
     (5.67704817454536, (1, 1, 2))],
    [(-2.00713963067187, (2, 1, 1)),
     (-2.00713963067187, (0, 3, 1)),
     (2.67618617422916, (0, 1, 3))],
    [(0.317356640745613, (4, 0, 0)),
     (0.634713281491226, (2, 2, 0)),
     (-2.53885312596490, (2, 0, 2)),
     (0.317356640745613, (0, 4, 0)),
     (-2.53885312596490, (0, 2, 2)),
     (0.846284375321634, (0, 0, 4))],
    [(-2.00713963067187, (3, 0, 1)),
     (-2.00713963067187, (1, 2, 1)),
     (2.67618617422916, (1, 0, 3))],
    [(0.473087347878780, (0, 4, 0)),
     (-2.83852408727268, (0, 2, 2)),
     (-0.473087347878780, (4, 0, 0)),
     (2.83852408727268, (2, 0, 2))],
    [(-5.31039230933979, (1, 2, 1)),
     (1.77013076977993, (3, 0, 1))],
    [(0.625835735449176, (0, 4, 0)),
     (-3.75501441269506, (2, 2, 0)),
     (0.625835735449176, (4, 0, 0))],
    # h, l=5:
    [(0.656382056840170, (0, 5, 0)),
     (-6.56382056840170, (2, 3, 0)),
     (3.28191028420085, (4, 1, 0))],
    [(-8.30264925952416, (1, 3, 1)),
     (8.30264925952416, (3, 1, 1))],
    [(-0.978476598870501, (2, 3, 0)),
     (0.489238299435250, (0, 5, 0)),
     (-3.91390639548200, (0, 3, 2)),
     (-1.46771489830575, (4, 1, 0)),
     (11.7417191864460, (2, 1, 2))],
    [(-4.79353678497332, (3, 1, 1)),
     (-4.79353678497332, (1, 3, 1)),
     (9.58707356994665, (1, 1, 3))],
    [(0.452946651195697, (4, 1, 0)),
     (0.905893302391394, (2, 3, 0)),
     (-5.43535981434836, (2, 1, 2)),
     (0.452946651195697, (0, 5, 0)),
     (-5.43535981434836, (0, 3, 2)),
     (3.62357320956558, (0, 1, 4))],
    [(1.75425483680135, (4, 0, 1)),
     (3.50850967360271, (2, 2, 1)),
     (-4.67801289813694, (2, 0, 3)),
     (1.75425483680135, (0, 4, 1)),
     (-4.67801289813694, (0, 2, 3)),
     (0.935602579627389, (0, 0, 5))],
    [(0.452946651195697, (5, 0, 0)),
     (0.905893302391394, (3, 2, 0)),
     (-5.43535981434836, (3, 0, 2)),
     (0.452946651195697, (1, 4, 0)),
     (-5.43535981434836, (1, 2, 2)),
     (3.62357320956558, (1, 0, 4))],
    [(2.39676839248666, (0, 4, 1)),
     (-4.79353678497332, (0, 2, 3)),
     (-2.39676839248666, (4, 0, 1)),
     (4.79353678497332, (2, 0, 3))],
    [(0.978476598870501, (3, 2, 0)),
     (1.46771489830575, (1, 4, 0)),
     (-11.7417191864460, (1, 2, 2)),
     (-0.489238299435250, (5, 0, 0)),
     (3.91390639548200, (3, 0, 2))],
    [(2.07566231488104, (0, 4, 1)),
     (-12.4539738892862, (2, 2, 1)),
     (2.07566231488104, (4, 0, 1))],
    [(3.28191028420085, (1, 4, 0)),
     (-6.56382056840170, (3, 2, 0)),
     (0.656382056840170, (5, 0, 0))],
    # i, l=6:
    [(4.09910463115149, (1, 5, 0)),
     (-13.6636821038383, (3, 3, 0)),
     (4.09910463115149, (5, 1, 0))],
    [(2.36661916223175, (0, 5, 1)),
     (-23.6661916223175, (2, 3, 1)),
     (11.8330958111588, (4, 1, 1))],
    [(2.01825960291490, (1, 5, 0)),
     (-20.1825960291490, (1, 3, 2)),
     (-2.01825960291490, (5, 1, 0)),
     (20.1825960291490, (3, 1, 2))],
    [(-5.52723155708954, (2, 3, 1)),
     (2.76361577854477, (0, 5, 1)),
     (-7.36964207611939, (0, 3, 3)),
     (-8.29084733563431, (4, 1, 1)),
     (22.1089262283582, (2, 1, 3))],
    [(0.921205259514923, (5, 1, 0)),
     (1.84241051902985, (3, 3, 0)),
     (-14.7392841522388, (3, 1, 2)),
     (0.921205259514923, (1, 5, 0)),
     (-14.7392841522388, (1, 3, 2)),
     (14.7392841522388, (1, 1, 4))],
    [(2.91310681259366, (4, 1, 1)),
     (5.82621362518731, (2, 3, 1)),
     (-11.6524272503746, (2, 1, 3)),
     (2.91310681259366, (0, 5, 1)),
     (-11.6524272503746, (0, 3, 3)),
     (4.66097090014985, (0, 1, 5))],
    [(-0.317846011338142, (6, 0, 0)),
     (-0.953538034014426, (4, 2, 0)),
     (5.72122820408656, (4, 0, 2)),
     (-0.953538034014426, (2, 4, 0)),
     (11.4424564081731, (2, 2, 2)),
     (-7.62830427211541, (2, 0, 4)),
     (-0.317846011338142, (0, 6, 0)),
     (5.72122820408656, (0, 4, 2)),
     (-7.62830427211541, (0, 2, 4)),
     (1.01710723628205, (0, 0, 6))],
    [(2.91310681259366, (5, 0, 1)),
     (5.82621362518731, (3, 2, 1)),
     (-11.6524272503746, (3, 0, 3)),
     (2.91310681259366, (1, 4, 1)),
     (-11.6524272503746, (1, 2, 3)),
     (4.66097090014985, (1, 0, 5))],
    [(-0.460602629757462, (2, 4, 0)),
     (-0.460602629757462, (0, 6, 0)),
     (7.36964207611939, (0, 4, 2)),
     (-7.36964207611939, (0, 2, 4)),
     (0.460602629757462, (6, 0, 0)),
     (0.460602629757462, (4, 2, 0)),
     (-7.36964207611939, (4, 0, 2)),
     (7.36964207611939, (2, 0, 4))],
    [(5.52723155708954, (3, 2, 1)),
     (8.29084733563431, (1, 4, 1)),
     (-22.1089262283582, (1, 2, 3)),
     (-2.76361577854477, (5, 0, 1)),
     (7.36964207611939, (3, 0, 3))],
    [(2.52282450364362, (2, 4, 0)),
     (-0.504564900728724, (0, 6, 0)),
     (5.04564900728724, (0, 4, 2)),
     (2.52282450364362, (4, 2, 0)),
     (-30.2738940437234, (2, 2, 2)),
     (-0.504564900728724, (6, 0, 0)),
     (5.04564900728724, (4, 0, 2))],
    [(11.8330958111588, (1, 4, 1)),
     (-23.6661916223175, (3, 2, 1)),
     (2.36661916223175, (5, 0, 1))],
    [(-0.683184105191914, (0, 6, 0)),
     (10.2477615778787, (2, 4, 0)),
     (-10.2477615778787, (4, 2, 0)),
     (0.683184105191914, (6, 0, 0))],
    # j, l=7:
    [(-0.707162732524596, (0, 7, 0)),
     (14.8504173830165, (2, 5, 0)),
     (-24.7506956383609, (4, 3, 0)),
     (4.95013912767217, (6, 1, 0))],
    [(15.8757639708114, (1, 5, 1)),
     (-52.9192132360380, (3, 3, 1)),
     (15.8757639708114, (5, 1, 1))],
    [(4.67024020848234, (2, 5, 0)),
     (-0.518915578720260, (0, 7, 0)),
     (6.22698694464312, (0, 5, 2)),
     (2.59457789360130, (4, 3, 0)),
     (-62.2698694464312, (2, 3, 2)),
     (-2.59457789360130, (6, 1, 0)),
     (31.1349347232156, (4, 1, 2))],
    [(12.4539738892862, (1, 5, 1)),
     (-41.5132462976208, (1, 3, 3)),
     (-12.4539738892862, (5, 1, 1)),
     (41.5132462976208, (3, 1, 3))],
    [(2.34688400793441, (4, 3, 0)),
     (0.469376801586882, (2, 5, 0)),
     (-18.7750720634753, (2, 3, 2)),
     (-0.469376801586882, (0, 7, 0)),
     (9.38753603173764, (0, 5, 2)),
     (-12.5167147089835, (0, 3, 4)),
     (1.40813040476065, (6, 1, 0)),
     (-28.1626080952129, (4, 1, 2)),
     (37.5501441269506, (2, 1, 4))],
    [(6.63799038667474, (5, 1, 1)),
     (13.2759807733495, (3, 3, 1)),
     (-35.4026153955986, (3, 1, 3)),
     (6.63799038667474, (1, 5, 1)),
     (-35.4026153955986, (1, 3, 3)),
     (21.2415692373592, (1, 1, 5))],
    [(-0.451658037912587, (6, 1, 0)),
     (-1.35497411373776, (4, 3, 0)),
     (10.8397929099021, (4, 1, 2)),
     (-1.35497411373776, (2, 5, 0)),
     (21.6795858198042, (2, 3, 2)),
     (-21.6795858198042, (2, 1, 4)),
     (-0.451658037912587, (0, 7, 0)),
     (10.8397929099021, (0, 5, 2)),
     (-21.6795858198042, (0, 3, 4)),
     (5.78122288528111, (0, 1, 6))],
    [(-2.38994969192017, (6, 0, 1)),
     (-7.16984907576052, (4, 2, 1)),
     (14.3396981515210, (4, 0, 3)),
     (-7.16984907576052, (2, 4, 1)),
     (28.6793963030421, (2, 2, 3)),
     (-11.4717585212168, (2, 0, 5)),
     (-2.38994969192017, (0, 6, 1)),
     (14.3396981515210, (0, 4, 3)),
     (-11.4717585212168, (0, 2, 5)),
     (1.09254843059208, (0, 0, 7))],
    [(-0.451658037912587, (7, 0, 0)),
     (-1.35497411373776, (5, 2, 0)),
     (10.8397929099021, (5, 0, 2)),
     (-1.35497411373776, (3, 4, 0)),
     (21.6795858198042, (3, 2, 2)),
     (-21.6795858198042, (3, 0, 4)),
     (-0.451658037912587, (1, 6, 0)),
     (10.8397929099021, (1, 4, 2)),
     (-21.6795858198042, (1, 2, 4)),
     (5.78122288528111, (1, 0, 6))],
    [(-3.31899519333737, (2, 4, 1)),
     (-3.31899519333737, (0, 6, 1)),
     (17.7013076977993, (0, 4, 3)),
     (-10.6207846186796, (0, 2, 5)),
     (3.31899519333737, (6, 0, 1)),
     (3.31899519333737, (4, 2, 1)),
     (-17.7013076977993, (4, 0, 3)),
     (10.6207846186796, (2, 0, 5))],
    [(-0.469376801586882, (5, 2, 0)),
     (-2.34688400793441, (3, 4, 0)),
     (18.7750720634753, (3, 2, 2)),
     (-1.40813040476065, (1, 6, 0)),
     (28.1626080952129, (1, 4, 2)),
     (-37.5501441269506, (1, 2, 4)),
     (0.469376801586882, (7, 0, 0)),
     (-9.38753603173764, (5, 0, 2)),
     (12.5167147089835, (3, 0, 4))],
    [(15.5674673616078, (2, 4, 1)),
     (-3.11349347232156, (0, 6, 1)),
     (10.3783115744052, (0, 4, 3)),
     (15.5674673616078, (4, 2, 1)),
     (-62.2698694464312, (2, 2, 3)),
     (-3.11349347232156, (6, 0, 1)),
     (10.3783115744052, (4, 0, 3))],
    [(2.59457789360130, (3, 4, 0)),
     (-2.59457789360130, (1, 6, 0)),
     (31.1349347232156, (1, 4, 2)),
     (4.67024020848234, (5, 2, 0)),
     (-62.2698694464312, (3, 2, 2)),
     (-0.518915578720260, (7, 0, 0)),
     (6.22698694464312, (5, 0, 2))],
    [(-2.64596066180190, (0, 6, 1)),
     (39.6894099270285, (2, 4, 1)),
     (-39.6894099270285, (4, 2, 1)),
     (2.64596066180190, (6, 0, 1))],
    [(-4.95013912767217, (1, 6, 0)),
     (24.7506956383609, (3, 4, 0)),
     (-14.8504173830165, (5, 2, 0)),
     (0.707162732524596, (7, 0, 0))]]
