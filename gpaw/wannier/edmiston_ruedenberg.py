from math import pi
from typing import Any, Tuple

import numpy as np

import _gpaw
from .overlap import WannierOverlaps

Array2D = Any


class LocalizationNotConvergedError(Exception):
    """Error raised if maxiter exceeded."""


class WannierFunctions:
    def __init__(self, U_nn, centers, values):
        self.U_nn = U_nn
        self.centers = centers
        self.value = value


def localize(overlaps: WannierOverlaps,
             maxiter: int = 100,
             tolerance: float = 1e-5,
             verbose: bool = False) -> Tuple[Array2D, Array2D, float]:
    """"""
    if (np.diag(overlaps.cell.diag()) != overlaps.cell).any():
        raise NotImplementedError('An orthogonal cell is required')
    assert overlaps.monkhors_pack_size == (1, 1, 1)

    Z_nnc = np.empty((overlaps.nbands, overlaps.nbands, 3), complex)
    for c, direction in enumerate([(1, 0, 0), (0, 1, 0), (0, 0, 1)]):
        Z_nnc[:, :, c] = overlaps.overlap(bz_index=0, direction=direction)

    U_nn = np.identity(overlaps.nbands)

    old = -1.0
    for iter in range(maxiter):
        value = _gpaw.localize(Z_nnc, U_nn)
        if verbose:
            print(iter, value)
        if value - old < tolerance:
            break
    else:
        raise LocalizationNotConvergedError(
            f'Did not converge in {maxiter} iterations')

    # Find centers:
    scaled_nc = -np.angle(Z_nnc.diagonal()).T / (2 * pi)
    centers_nv = (scaled_nc % 1.0).dot(overlaps.cell)

    return WannierFunctions(U_nn, centers_nv, value)


if __name__ == '__main__':
    import sys
    from gpaw import GPAW
    from .overlaps import calculate_overlaps

    calc = GPAW(sys.argv[1])
    overlaps = calculate_overlaps(calc)
    wan = localize(overlaps, verbose=True)
    print(wan.centers)
