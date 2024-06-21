from dataclasses import dataclass
from collections.abc import Sequence
import numpy as np

from gpaw.cgpaw import GG_shuffle

from gpaw.response.symmetry import QSymmetries
from gpaw.response.pair_functions import SingleQPWDescriptor


@dataclass
class QSymmetryOperators(Sequence):
    symmetries: QSymmetries

    def __len__(self):
        return len(self.symmetries)


@dataclass
class HeadSymmetryOperators(QSymmetryOperators):
    cell_cv: np.ndarray
    icell_cv: np.ndarray

    @classmethod
    def from_gd(cls, symmetries, gd):
        return cls(symmetries, gd.cell_cv, gd.icell_cv)

    def __getitem__(self, s):
        U_cc, sign, _ = self.symmetries[s]
        M_vv = self.cell_cv.T @ U_cc.T @ self.icell_cv
        return M_vv, sign

    def symmetrize_wvv(self, A_wvv):
        """Symmetrize chi0_wvv"""
        tmp_wvv = np.zeros_like(A_wvv)
        for M_vv, sign in self:
            tmp = np.dot(np.dot(M_vv.T, A_wvv), M_vv)
            if sign == 1:
                tmp_wvv += np.transpose(tmp, (1, 0, 2))
            elif sign == -1:  # transpose head
                tmp_wvv += np.transpose(tmp, (1, 2, 0))
        # Overwrite the input
        A_wvv[:] = tmp_wvv / len(self)


@dataclass
class BodySymmetryOperators(QSymmetryOperators):
    qpd: SingleQPWDescriptor

    def __post_init__(self):
        self.G_sG = initialize_G_maps(self.symmetries, self.qpd)

    def __getitem__(self, s):
        return self.G_sG[s], self.symmetries.sign_s[s]

    def symmetrize_wGG(self, A_wGG):
        """Symmetrize an array in GG'."""
        for A_GG in A_wGG:
            tmp_GG = np.zeros_like(A_GG, order='C')
            for G_G, sign in self:
                # Numpy:
                # if sign == 1:
                #     tmp_GG += A_GG[G_G, :][:, G_G]
                # if sign == -1:
                #     tmp_GG += A_GG[G_G, :][:, G_G].T
                # C:
                GG_shuffle(G_G, sign, A_GG, tmp_GG)
            A_GG[:] = tmp_GG / len(self)

    # Set up complex frequency alias
    symmetrize_zGG = symmetrize_wGG


@dataclass
class WingSymmetryOperators(QSymmetryOperators):
    qpd: SingleQPWDescriptor

    def __post_init__(self):
        self.head_operators = HeadSymmetryOperators.from_gd(
            self.symmetries, self.qpd.gd)
        self.G_sG = initialize_G_maps(self.symmetries, self.qpd)

    def __getitem__(self, s):
        M_vv, sign = self.head_operators[s]
        return M_vv, sign, self.G_sG[s]

    def symmetrize_wvv(self, *args):
        self.head_operators.symmetrize_wvv(*args)

    def symmetrize_wxvG(self, A_wxvG):
        """Symmetrize chi0_wxvG"""
        tmp_wxvG = np.zeros_like(A_wxvG)
        for M_vv, sign, G_G in self:
            if sign == 1:
                tmp = sign * np.dot(M_vv.T, A_wxvG[..., G_G])
            elif sign == -1:  # transpose wings
                tmp = sign * np.dot(M_vv.T, A_wxvG[:, ::-1, :, G_G])
            tmp_wxvG += np.transpose(tmp, (1, 2, 0, 3))
        # Overwrite the input
        A_wxvG[:] = tmp_wxvG / len(self)


def initialize_G_maps(symmetries: QSymmetries, qpd: SingleQPWDescriptor):
    """Calculate the Gvector mappings."""
    assert np.allclose(symmetries.q_c, qpd.q_c)
    B_cv = 2.0 * np.pi * qpd.gd.icell_cv
    G_Gv = qpd.get_reciprocal_vectors(add_q=False)
    G_Gc = np.dot(G_Gv, np.linalg.inv(B_cv))
    Q_G = qpd.Q_qG[0]

    G_sG = []
    for U_cc, sign, shift_c in symmetries:
        iU_cc = np.linalg.inv(U_cc).T
        UG_Gc = np.dot(G_Gc - shift_c, sign * iU_cc)

        assert np.allclose(UG_Gc.round(), UG_Gc)
        UQ_G = np.ravel_multi_index(UG_Gc.round().astype(int).T,
                                    qpd.gd.N_c, 'wrap')

        G_G = len(Q_G) * [None]
        for G, UQ in enumerate(UQ_G):
            try:
                G_G[G] = np.argwhere(Q_G == UQ)[0][0]
            except IndexError:
                print('This should not be possible but' +
                      'a G-vector was mapped outside the sphere')
                raise IndexError
        G_sG.append(np.array(G_G, dtype=np.int32))
    return np.array(G_sG)
