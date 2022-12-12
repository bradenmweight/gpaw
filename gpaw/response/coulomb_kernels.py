"""This module defines Coulomb and XC kernels for the response model.
"""

import numpy as np
from ase.dft import monkhorst_pack


class CoulombKernel:
    def __init__(self, truncation, gs):
        from gpaw.response.wstc import WignerSeitzTruncatedCoulomb
        self.truncation = truncation
        assert self.truncation in {None, 'wigner-seitz', '2D'}
        self._gs = gs

        if self.truncation == 'wigner-seitz':
            self._wstc = WignerSeitzTruncatedCoulomb(self._gs.gd.cell_cv,
                                                     self._gs.kd.N_c)
        else:
            self._wstc = None

    def description(self):
        if self.truncation is None:
            return 'No Coulomb truncation'
        elif self._wstc is not None:
            return '\n'.join(
                ['Wigner–Seitz truncated Coulomb truncation',
                 self._wstc.get_description()])
        else:
            return f'Using {self.truncation} Coulomb truncation'

    def sqrtV(self, pd, q_v):
        return self.V(pd, q_v)**0.5

    def V(self, pd, q_v):
        return get_coulomb_kernel(
            pd, self._gs.kd.N_c, q_v=q_v,
            truncation=self.truncation,
            wstc=self._wstc)


def get_coulomb_kernel(pd, N_c, truncation=None, q_v=None, wstc=None):
    """Factory function that calls the specified flavour
    of the Coulomb interaction"""

    qG_Gv = pd.get_reciprocal_vectors(add_q=True)
    if q_v is not None:
        assert pd.kd.gamma
        qG_Gv += q_v

    if truncation is None:
        if pd.kd.gamma and q_v is None:
            v_G = np.zeros(len(pd.G2_qG[0]))
            v_G[0] = 4 * np.pi
            v_G[1:] = 4 * np.pi / (qG_Gv[1:]**2).sum(axis=1)
        else:
            v_G = 4 * np.pi / (qG_Gv**2).sum(axis=1)

    elif truncation == '2D':
        v_G = calculate_2D_truncated_coulomb(pd, q_v=q_v, N_c=N_c)
        if pd.kd.gamma and q_v is None:
            v_G[0] = 0.0

    elif truncation in {'1D', '0D'}:
        raise feature_removed()

    elif truncation == 'wigner-seitz':
        v_G = wstc.get_potential(pd, q_v=q_v)
        if pd.kd.gamma and q_v is None:
            v_G[0] = 0.0
    else:
        raise ValueError('Truncation scheme %s not implemented' % truncation)

    return v_G.astype(complex)


def calculate_2D_truncated_coulomb(pd, q_v=None, N_c=None):
    """ Simple 2D truncation of Coulomb kernel PRB 73, 205119.
    The non-periodic direction is determined from k-point grid.
    """

    qG_Gv = pd.get_reciprocal_vectors(add_q=True)
    if pd.kd.gamma:
        if q_v is not None:
            qG_Gv += q_v
        else:  # only to avoid warning. Later set to zero in factory function
            qG_Gv[0] = [1., 1., 1.]

    # The non-periodic direction is determined from k-point grid
    Nn_c = np.where(N_c == 1)[0]
    Np_c = np.where(N_c != 1)[0]
    if len(Nn_c) != 1:
        # The k-point grid does not fit with boundary conditions
        Nn_c = [2]  # choose reduced cell vectors 0, 1
        Np_c = [0, 1]  # choose reduced cell vector 2
    # Truncation length is half of cell vector in non-periodic direction
    R = pd.gd.cell_cv[Nn_c[0], Nn_c[0]] / 2.

    qGp_G = ((qG_Gv[:, Np_c[0]])**2 + (qG_Gv[:, Np_c[1]]**2))**0.5
    qGn_G = qG_Gv[:, Nn_c[0]]

    v_G = 4 * np.pi / (qG_Gv**2).sum(axis=1)
    if np.allclose(qGn_G[0], 0) or pd.kd.gamma:
        """sin(qGn_G * R) = 0 when R = L/2 and q_n = 0.0"""
        v_G *= 1.0 - np.exp(-qGp_G * R) * np.cos(qGn_G * R)
    else:
        """Normal component of q is not zero"""
        a_G = qGn_G / qGp_G * np.sin(qGn_G * R) - np.cos(qGn_G * R)
        v_G *= 1. + np.exp(-qGp_G * R) * a_G

    return v_G.astype(complex)


def get_integrated_kernel(pd, N_c, truncation=None, N=100, reduced=False):
    from scipy.special import j1, k0, j0, k1

    B_cv = 2 * np.pi * pd.gd.icell_cv
    Nf_c = np.array([N, N, N])
    if reduced:
        # Only integrate periodic directions if truncation is used
        Nf_c[np.where(N_c == 1)[0]] = 1
    q_qc = monkhorst_pack(Nf_c) / N_c
    q_qc += pd.kd.ibzk_kc[0]
    q_qv = np.dot(q_qc, B_cv)

    if truncation is None:
        V_q = 4 * np.pi / np.sum(q_qv**2, axis=1)
    elif truncation == '2D':
        # The non-periodic direction is determined from k-point grid
        Nn_c = np.where(N_c == 1)[0]
        Np_c = np.where(N_c != 1)[0]
        if len(Nn_c) != 1:
            # The k-point grid does not fit with boundary conditions
            Nn_c = [2]  # choose reduced cell vectors 0, 1
            Np_c = [0, 1]  # choose reduced cell vector 2
        # Truncation length is half of cell vector in non-periodic direction
        R = pd.gd.cell_cv[Nn_c[0], Nn_c[0]] / 2.

        qp_q = ((q_qv[:, Np_c[0]])**2 + (q_qv[:, Np_c[1]]**2))**0.5
        qn_q = q_qv[:, Nn_c[0]]

        V_q = 4 * np.pi / (q_qv**2).sum(axis=1)
        a_q = qn_q / qp_q * np.sin(qn_q * R) - np.cos(qn_q * R)
        V_q *= 1. + np.exp(-qp_q * R) * a_q
    elif truncation == '1D':
        # The non-periodic direction is determined from k-point grid
        Nn_c = np.where(N_c == 1)[0]
        Np_c = np.where(N_c != 1)[0]

        if len(Nn_c) != 2:
            # The k-point grid does not fit with boundary conditions
            Nn_c = [0, 1]    # Choose reduced cell vectors 0, 1
            Np_c = [2]       # Choose reduced cell vector 2
        # The radius is determined from area of non-periodic part of cell
        Acell_cv = pd.gd.cell_cv[Nn_c, :][:, Nn_c]
        R = abs(np.linalg.det(Acell_cv) / np.pi)**0.5

        qnR_q = (q_qv[:, Nn_c[0]]**2 + q_qv[:, Nn_c[1]]**2)**0.5 * R
        qpR_q = abs(q_qv[:, Np_c[0]]) * R
        V_q = 4 * np.pi / (q_qv**2).sum(axis=1)
        V_q *= (1.0 + qnR_q * j1(qnR_q) * k0(qpR_q)
                - qpR_q * j0(qnR_q) * k1(qpR_q))
    elif truncation == '0D' or 'wigner-seitz':
        R = (3 * pd.gd.volume / (4 * np.pi))**(1. / 3.)
        q2_q = (q_qv**2).sum(axis=1)
        V_q = 4 * np.pi / q2_q
        V_q *= 1.0 - np.cos(q2_q**0.5 * R)

    return np.sum(V_q) / len(V_q), np.sum(V_q**0.5) / len(V_q)


def feature_removed():
    return RuntimeError(
        '0D and 1D truncation have been removed due to not being tested.  '
        'If you need them, please find them in '
        'ec9e49e25613bb99cd69eec9d2613e38b9f6e6e1 '
        'and make sure to add tests in order to have them re-added.')
