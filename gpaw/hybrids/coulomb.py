from math import pi

import numpy as np

from io import StringIO
from gpaw.response.wstc import WignerSeitzTruncatedCoulomb as WSTC


def coulomb_inteaction(omega, gd, kd):
    if omega:
        return ShortRangeCoulomb(omega)
    # Wigner-Seitz truncated Coulomb:
    output = StringIO()
    coulomb = WSTC(gd.cell_cv, kd.N_c, txt=output)
    coulomb.description = output.getvalue()
    return coulomb


class ShortRangeCoulomb:
    def __init__(self, omega):
        self.omega = omega
        self.description = f'Short-range Coulomb (omega={omega} bohr^-1)'

    def get_potential(self, pd):
        G2_G = pd.G2_qG[0]
        x_G = 1 - np.exp(-G2_G / (4 * self.omega**2))
        with np.errstate(invalid='ignore'):
            v_G = 4 * pi * x_G / G2_G
        if pd.kd.gamma:
            v_G[0] = pi / self.omega**2
        return v_G
