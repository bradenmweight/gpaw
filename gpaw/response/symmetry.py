from typing import Union
from dataclasses import dataclass
from collections.abc import Sequence
from functools import cached_property

import numpy as np
from scipy.spatial import Delaunay, cKDTree

from gpaw.bztools import get_reduced_bz, unique_rows
from gpaw.cgpaw import GG_shuffle


class KPointFinder:
    def __init__(self, bzk_kc):
        self.kdtree = cKDTree(self._round(bzk_kc))

    @staticmethod
    def _round(bzk_kc):
        return np.mod(np.mod(bzk_kc, 1).round(6), 1)

    def find(self, kpt_c):
        distance, k = self.kdtree.query(self._round(kpt_c))
        if distance > 1.e-6:
            raise ValueError('Requested k-point is not on the grid. '
                             'Please check that your q-points of interest '
                             'are commensurate with the k-point grid.')

        return k


@dataclass
class QSymmetries(Sequence):
    """Symmetry operations for a given q-point.

    We operate with several different symmetry indices:
      * u: indices the unitary symmetries of the system. Length is nU.
      * S: extended symmetry index. In addition to the unitary symmetries
           (first nU indices) it includes also symmetries generated by a
           unitary symmetry transformation *followed* by a time-reversal.
           Length is 2 * nU.
      * s: reduced symmetry index. Includes all the "S-symmetries" which map
           the q-point in question onto itself (up to a reciprocal lattice
           vector). May be reduced further, if some of the symmetries have been
           disabled. Length is q-dependent and depends on user input.
    """
    q_c: np.ndarray
    U_ucc: np.ndarray  # unitary symmetry transformations
    S_s: np.ndarray  # extended symmetry index for each q-symmetry
    shift_sc: np.ndarray  # reciprocal lattice shifts, G = (T)Uq - q

    def __post_init__(self):
        self.nU = len(self.U_ucc)

    def __len__(self):
        return len(self.S_s)

    def __getitem__(self, s):
        S = self.S_s[s]
        return self.unioperator(S), self.sign(S), self.shift_sc[s]

    def unioperator(self, S):
        return self.U_ucc[S % self.nU]

    def timereversal(self, S):
        """Does the extended index S involve a time-reversal symmetry?"""
        return bool(S // self.nU)

    def sign(self, S):
        """Flip the sign under time-reversal."""
        if self.timereversal(S):
            return -1
        return 1

    @cached_property
    def ndirect(self):
        """Number of direct symmetries."""
        return sum(np.array(self.S_s) < self.nU)

    @property
    def nindirect(self):
        """Number of indirect symmetries."""
        return len(self) - self.ndirect

    def description(self) -> str:
        """Return string description of symmetry operations."""
        isl = ['\n']
        nx = 6  # You are not allowed to use non-symmorphic syms (value 3)
        y = 0
        for y in range((len(self) + nx - 1) // nx):
            for c in range(3):
                tisl = []
                for x in range(nx):
                    s = x + y * nx
                    if s == len(self):
                        break
                    U_cc, sign, _ = self[s]
                    op_c = sign * U_cc[c]
                    tisl.append(f'  ({op_c[0]:2d} {op_c[1]:2d} {op_c[2]:2d})')
                tisl.append('\n')
                isl.append(''.join(tisl))
            isl.append('\n')
        return ''.join(isl[:-1])


@dataclass
class QSymmetryAnalyzer:
    """Identifies symmetries of the k-grid, under which q is invariant.

    Parameters
    ----------
    point_group : bool
        Use point group symmetry.
    time_reversal : bool
        Use time-reversal symmetry (if applicable).
    """
    point_group: bool = True
    time_reversal: bool = True

    @property
    def disabled(self):
        return not (self.point_group or self.time_reversal)

    @property
    def disabled_symmetry_info(self):
        if self.disabled:
            txt = ''
        elif not self.point_group:
            txt = 'point-group '
        elif not self.time_reversal:
            txt = 'time-reversal '
        else:
            return ''
        txt += 'symmetry has been manually disabled'
        return txt

    def analysis_info(self, symmetries):
        dsinfo = self.disabled_symmetry_info
        return '\n'.join([
            '',
            f'Symmetries of q_c{f" ({dsinfo})" if len(dsinfo) else ""}:',
            f'    Direct symmetries (Uq -> q): {symmetries.ndirect}',
            f'    Indirect symmetries (TUq -> q): {symmetries.nindirect}',
            f'In total {len(symmetries)} allowed symmetries.',
            symmetries.description()])

    def analyze(self, q_c, kpoints, context):
        """Analyze symmetries and set up KPointDomainGenerator."""
        symmetries = self.analyze_symmetries(q_c, kpoints.kd)
        generator = KPointDomainGenerator(symmetries, kpoints)
        context.print(self.analysis_info(symmetries))
        context.print(generator.get_infostring())
        return symmetries, generator

    def analyze_symmetries(self, q_c, kd):
        r"""Determine allowed symmetries.

        An direct symmetry U must fulfill::

          U \mathbf{q} = q + \Delta

        Under time-reversal (indirect) it must fulfill::

          -U \mathbf{q} = q + \Delta

        where :math:`\Delta` is a reciprocal lattice vector.
        """
        # Map q-point for each unitary symmetry
        U_ucc = kd.symmetry.op_scc  # here s is the unitary symmetry index
        Uq_uc = np.dot(U_ucc, q_c)

        # Direct and indirect -> global symmetries
        nU = len(U_ucc)
        nS = 2 * nU
        shift_Sc = np.zeros((nS, 3), int)
        is_qsymmetry_S = np.zeros(nS, bool)

        # Identify direct symmetries
        # Check whether U q - q is integer (reciprocal lattice vector)
        dshift_uc = Uq_uc - q_c[np.newaxis]
        is_direct_symm_u = (dshift_uc == dshift_uc.round()).all(axis=1)
        is_qsymmetry_S[:nU][is_direct_symm_u] = True
        shift_Sc[:nU] = dshift_uc

        # Identify indirect symmetries
        # Check whether -U q - q is integer (reciprocal lattice vector)
        idshift_uc = -Uq_uc - q_c
        is_indirect_symm_u = (idshift_uc == idshift_uc.round()).all(axis=1)
        is_qsymmetry_S[nU:][is_indirect_symm_u] = True
        shift_Sc[nU:] = idshift_uc

        # The indices of the allowed symmetries
        S_s = is_qsymmetry_S.nonzero()[0]

        # Set up symmetry filters
        def is_not_point_group(S):
            return (U_ucc[S % nU] == np.eye(3)).all()

        def is_not_time_reversal(S):
            return not bool(S // nU)

        def is_not_non_symmorphic(S):
            return not bool(kd.symmetry.ft_sc[S % nU].any())

        # Filter out point-group symmetry, if disabled
        if not self.point_group:
            S_s = list(filter(is_not_point_group, S_s))

        # Filter out time-reversal, if inapplicable or disabled
        if not kd.symmetry.time_reversal or \
           kd.symmetry.has_inversion or \
           not self.time_reversal:
            S_s = list(filter(is_not_time_reversal, S_s))

        # We always filter out non-symmorphic symmetries
        S_s = list(filter(is_not_non_symmorphic, S_s))

        return QSymmetries(q_c, U_ucc, S_s, shift_Sc[S_s])


QSymmetryInput = Union[QSymmetryAnalyzer, dict, bool]


def ensure_qsymmetry(qsymmetry: QSymmetryInput) -> QSymmetryAnalyzer:
    if not isinstance(qsymmetry, QSymmetryAnalyzer):
        if isinstance(qsymmetry, dict):
            qsymmetry = QSymmetryAnalyzer(**qsymmetry)
        else:
            qsymmetry = QSymmetryAnalyzer(
                point_group=qsymmetry, time_reversal=qsymmetry)
    return qsymmetry


class KPointDomainGenerator:
    def __init__(self, symmetries, kpoints):
        self.symmetries = symmetries

        self.kd = kpoints.kd
        self.kptfinder = kpoints.kptfinder

    def how_many_symmetries(self):
        # temporary backwards compatibility for external calls
        return len(self.symmetries)

    def get_infostring(self):
        # Maybe we can avoid calling this somehow, we're only using
        # it to print:
        K_gK = self.group_kpoints()
        ng = len(K_gK)
        txt = f'{ng} groups of equivalent kpoints. '
        percent = (1. - (ng + 0.) / self.kd.nbzkpts) * 100
        txt += f'{percent}% reduction.\n'
        return txt

    def group_kpoints(self, K_k=None):
        """Group kpoints according to the reduced symmetries"""
        if K_k is None:
            K_k = np.arange(self.kd.nbzkpts)
        bz2bz_kS = self.kd.bz2bz_ks  # on kd, s is the global symmetry index
        nk = len(bz2bz_kS)
        sbz2sbz_ks = bz2bz_kS[K_k][:, self.symmetries.S_s]  # s: q-symmetries
        # Avoid -1 (see documentation in gpaw.symmetry)
        sbz2sbz_ks[sbz2sbz_ks == -1] = nk

        smallestk_k = np.sort(sbz2sbz_ks)[:, 0]
        k2g_g = np.unique(smallestk_k, return_index=True)[1]

        K_gs = sbz2sbz_ks[k2g_g]
        K_gK = [np.unique(K_s[K_s != nk]) for K_s in K_gs]

        return K_gK

    def get_kpt_domain(self):
        k_kc = np.array([self.kd.bzk_kc[K_K[0]] for
                         K_K in self.group_kpoints()])
        return k_kc

    def get_tetrahedron_ikpts(self, *, pbc_c, cell_cv):
        """Find irreducible k-points for tetrahedron integration."""
        U_scc = np.array([  # little group of q
            sign * U_cc for U_cc, sign, _ in self.symmetries])

        # Determine the irreducible BZ
        bzk_kc, ibzk_kc, _ = get_reduced_bz(cell_cv,
                                            U_scc,
                                            False,
                                            pbc_c=pbc_c)

        n = 3
        N_xc = np.indices((n, n, n)).reshape((3, n**3)).T - n // 2

        # Find the irreducible kpoints
        tess = Delaunay(ibzk_kc)
        ik_kc = []
        for N_c in N_xc:
            k_kc = self.kd.bzk_kc + N_c
            k_kc = k_kc[tess.find_simplex(k_kc) >= 0]
            if not len(ik_kc) and len(k_kc):
                ik_kc = unique_rows(k_kc)
            elif len(k_kc):
                ik_kc = unique_rows(np.append(k_kc, ik_kc, axis=0))

        return ik_kc

    def get_tetrahedron_kpt_domain(self, *, pbc_c, cell_cv):
        ik_kc = self.get_tetrahedron_ikpts(pbc_c=pbc_c, cell_cv=cell_cv)
        if pbc_c.all():
            k_kc = ik_kc
        else:
            k_kc = np.append(ik_kc,
                             ik_kc + (~pbc_c).astype(int),
                             axis=0)
        return k_kc

    def get_kpoint_weight(self, k_c):
        K = self.kptfinder.find(k_c)
        iK = self.kd.bz2ibz_k[K]
        K_k = self.unfold_ibz_kpoint(iK)
        K_gK = self.group_kpoints(K_k)

        for K_k in K_gK:
            if K in K_k:
                return len(K_k)

    def unfold_ibz_kpoint(self, ik):
        """Return kpoints related to irreducible kpoint."""
        kd = self.kd
        K_k = np.unique(kd.bz2bz_ks[kd.ibz2bz_k[ik]])
        K_k = K_k[K_k != -1]
        return K_k


@dataclass
class HeadSymmetryOperators(Sequence):
    symmetries: QSymmetries
    cell_cv: np.ndarray
    icell_cv: np.ndarray

    @classmethod
    def from_gd(cls, symmetries, gd):
        return cls(symmetries, gd.cell_cv, gd.icell_cv)

    def __len__(self):
        return len(self.symmetries)

    def __getitem__(self, s):
        U_cc, sign, _ = self.symmetries[s]
        M_vv = self.cell_cv.T @ U_cc.T @ self.icell_cv
        return M_vv, sign

    def symmetrize_wvv(self, A_wvv):
        tmp_wvv = np.zeros_like(A_wvv)
        for M_vv, sign in self:
            tmp = np.dot(np.dot(M_vv.T, A_wvv), M_vv)
            if sign == 1:
                tmp_wvv += np.transpose(tmp, (1, 0, 2))
            elif sign == -1:  # transpose head
                tmp_wvv += np.transpose(tmp, (1, 2, 0))
        # Overwrite the input
        A_wvv[:] = tmp_wvv / len(self)


class PWSymmetrizer:
    def __init__(self, symmetries: QSymmetries, qpd):
        assert np.allclose(symmetries.q_c, qpd.q_c)
        self.symmetries = symmetries
        self.qpd = qpd
        self.head_operators = HeadSymmetryOperators.from_gd(symmetries, qpd.gd)
        self.G_sG = self.initialize_G_maps()

    def symmetrize_wvv(self, A_wvv):
        self.head_operators.symmetrize_wvv(A_wvv)

    def symmetrize_wGG(self, A_wGG):
        """Symmetrize an array in GG'."""

        for A_GG in A_wGG:
            tmp_GG = np.zeros_like(A_GG, order='C')
            # tmp2_GG = np.zeros_like(A_GG)

            for (_, sign, _), G_G in zip(self.symmetries, self.G_sG):
                GG_shuffle(G_G, sign, A_GG, tmp_GG)

                # This is the exact operation that GG_shuffle does.
                # Uncomment lines involving tmp2_GG to test the
                # implementation in action:
                #
                # if sign == 1:
                #     tmp2_GG += A_GG[G_G, :][:, G_G]
                # if sign == -1:
                #     tmp2_GG += A_GG[G_G, :][:, G_G].T

            # assert np.allclose(tmp_GG, tmp2_GG)
            A_GG[:] = tmp_GG / len(self.symmetries)

    # Set up complex frequency alias
    symmetrize_zGG = symmetrize_wGG

    def symmetrize_wxvG(self, A_wxvG):
        """Symmetrize chi0_wxvG"""
        A_cv = self.qpd.gd.cell_cv
        iA_cv = self.qpd.gd.icell_cv

        tmp_wxvG = np.zeros_like(A_wxvG)
        for (U_cc, sign, _), G_G in zip(self.symmetries, self.G_sG):
            M_vv = np.dot(np.dot(A_cv.T, U_cc.T), iA_cv)
            if sign == 1:
                tmp = sign * np.dot(M_vv.T, A_wxvG[..., G_G])
            elif sign == -1:  # transpose wings
                tmp = sign * np.dot(M_vv.T, A_wxvG[:, ::-1, :, G_G])
            tmp_wxvG += np.transpose(tmp, (1, 2, 0, 3))

        # Overwrite the input
        A_wxvG[:] = tmp_wxvG / len(self.symmetries)

    def initialize_G_maps(self):
        """Calculate the Gvector mappings."""
        qpd = self.qpd
        B_cv = 2.0 * np.pi * qpd.gd.icell_cv
        G_Gv = qpd.get_reciprocal_vectors(add_q=False)
        G_Gc = np.dot(G_Gv, np.linalg.inv(B_cv))
        Q_G = qpd.Q_qG[0]

        G_sG = []
        for U_cc, sign, shift_c in self.symmetries:
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
