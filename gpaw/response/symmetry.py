from typing import Union
from dataclasses import dataclass

import numpy as np
from scipy.spatial import Delaunay, cKDTree

from gpaw.bztools import get_reduced_bz, unique_rows
from gpaw.cgpaw import GG_shuffle

from gpaw.response import timer


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
class QSymmetries:
    U_ucc: np.ndarray  # u: unitary symmetry index
    S_s: np.ndarray  # S: global symmetry index, s: q-symmetry index
    shift_Sc: np.ndarray

    def __post_init__(self):
        self.nU = len(self.U_ucc)

    def __len__(self):
        return len(self.S_s)

    def timereversal(self, S):
        """Is the global index S a time-reversal symmetry?"""
        return bool(S // self.nU)

    def sign(self, S):
        """Flip the sign under time-reversal."""
        if self.timereversal(S):
            return -1
        return 1

    def get_symmetry_operator(self, S):
        """Return symmetry operator s."""
        return self.U_ucc[S % self.nU], self.sign(S)


@dataclass
class QSymmetryAnalyzer:
    """K-point symmetry analyzer for transitions k -> k + q.

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

    def analyze(self, kpoints, qpd, context):
        symmetries = self.analyze_symmetries(qpd.q_c, kpoints.kd)
        return PWSymmetryAnalyzer(
            symmetries,
            kpoints, qpd, context, not self.point_group,
            not self.time_reversal)

    def analyze_symmetries(self, q_c, kd):
        r"""Determine allowed symmetries.

        An direct symmetry U must fulfill::

          U \mathbf{q} = q + \Delta

        Under time-reversal (indirect) it must fulfill::

          -U \mathbf{q} = q + \Delta

        where :math:`\Delta` is a reciprocal lattice vector.
        """
        # Map q-point for each unitary symmetry
        U_ucc = kd.symmetry.op_scc
        newq_uc = np.dot(U_ucc, q_c)

        # Direct and indirect -> global symmetries
        nU = len(U_ucc)
        nS = 2 * nU
        shift_Sc = np.zeros((nS, 3), int)
        is_qsymmetry_S = np.zeros(nS, bool)

        # Identify direct symmetries
        # Check whether U q - q is integer (reciprocal lattice vector)
        dshift_uc = newq_uc - q_c[np.newaxis]
        is_direct_symm_u = (dshift_uc == dshift_uc.round()).all(axis=1)
        is_qsymmetry_S[:nU][is_direct_symm_u] = True
        shift_Sc[:nU] = dshift_uc

        # Identify indirect symmetries
        # Check whether -U q - q is integer (reciprocal lattice vector)
        idshift_uc = -newq_uc - q_c
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

        return QSymmetries(U_ucc, S_s, shift_Sc)


QSymmetryInput = Union[QSymmetryAnalyzer, dict, bool]


def ensure_qsymmetry(qsymmetry: QSymmetryInput) -> QSymmetryAnalyzer:
    if not isinstance(qsymmetry, QSymmetryAnalyzer):
        if isinstance(qsymmetry, dict):
            qsymmetry = QSymmetryAnalyzer(**qsymmetry)
        else:
            qsymmetry = QSymmetryAnalyzer(
                point_group=qsymmetry, time_reversal=qsymmetry)
    return qsymmetry


class PWSymmetryAnalyzer:
    """Class for handling planewave symmetries."""

    def __init__(self, symmetries, kpoints, qpd, context,
                 disable_point_group=False,
                 disable_time_reversal=False):
        """Creates a PWSymmetryAnalyzer object.

        Determines which of the symmetries of the atomic structure
        that is compatible with the reciprocal lattice. Contains the
        necessary functions for mapping quantities between kpoints,
        and or symmetrizing arrays.

        kd: KPointDescriptor
            The kpoint descriptor containing the
            information about symmetries and kpoints.
        qpd: SingleQPWDescriptor
            Plane wave descriptor that contains the reciprocal
            lattice .
        context: ResponseContext
        disable_point_group: bool
            Switch for disabling point group symmetries.
        disable_time_reversal:
            Switch for disabling time reversal.
        """
        self.symmetries = symmetries
        self.qpd = qpd
        self.kd = kd = kpoints.kd
        self.context = context

        # Settings
        self.disable_point_group = disable_point_group
        self.disable_time_reversal = disable_time_reversal
        if (kd.symmetry.has_inversion or not kd.symmetry.time_reversal) and \
           not self.disable_time_reversal:
            self.context.print('\nThe ground calculation does not support time'
                               '-reversal symmetry possibly because it has an '
                               'inversion center or that it has been manually '
                               'deactivated.\n')
            self.disable_time_reversal = True

        self.disable_symmetries = (self.disable_point_group and
                                   self.disable_time_reversal)

        # Number of symmetries
        self.nsym = 2 * self.symmetries.nU
        self.use_time_reversal = not self.disable_time_reversal

        self.kptfinder = kpoints.kptfinder

        self.G_sG = self.initialize_G_maps()

        self.context.print(self.get_infostring())
        self.context.print(self.symmetry_description())

    def how_many_symmetries(self):
        # temporary backwards compatibility for external calls
        return len(self.symmetries)

    def get_infostring(self):
        txt = ''

        if self.disable_point_group:
            txt += 'Point group not included. '
        else:
            txt += 'Point group included. '

        if self.disable_time_reversal:
            txt += 'Time reversal not included. '
        else:
            txt += 'Time reversal included. '

        txt += 'Disabled non-symmorphic symmetries. '

        if self.disable_symmetries:
            txt += 'All symmetries have been disabled. '

        txt += f'Found {len(self.symmetries)} allowed symmetries. '

        # Maybe we can avoid calling this somehow, we're only using
        # it to print:
        K_gK = self.group_kpoints()
        ng = len(K_gK)
        txt += f'{ng} groups of equivalent kpoints. '
        percent = (1. - (ng + 0.) / self.kd.nbzkpts) * 100
        txt += f'{percent}% reduction. '
        return txt

    def symmetry_description(self) -> str:
        """Return string description of symmetry operations."""
        isl = ['\n']
        nx = 6  # You are not allowed to use non-symmorphic syms (value 3)
        ns = len(self.symmetries)
        y = 0
        for y in range((ns + nx - 1) // nx):
            for c in range(3):
                tisl = []
                for x in range(nx):
                    s = x + y * nx
                    if s == ns:
                        break
                    op_cc, sign = self.symmetries.get_symmetry_operator(
                        # little ugly this, symmetries can do the indexing XXX
                        self.symmetries.S_s[s])
                    op_c = sign * op_cc[c]
                    tisl.append(f'  ({op_c[0]:2d} {op_c[1]:2d} {op_c[2]:2d})')
                tisl.append('\n')
                isl.append(''.join(tisl))
            isl.append('\n')
        return ''.join(isl)

    @timer('Group kpoints')
    def group_kpoints(self, K_k=None):
        """Group kpoints according to the reduced symmetries"""
        if K_k is None:
            K_k = np.arange(self.kd.nbzkpts)
        S_s = self.symmetries.S_s
        bz2bz_ks = self.kd.bz2bz_ks
        nk = len(bz2bz_ks)
        sbz2sbz_ks = bz2bz_ks[K_k][:, S_s]  # Reduced number of symmetries
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

    def get_tetrahedron_ikpts(self, *, pbc_c):
        """Find irreducible k-points for tetrahedron integration."""
        # Get the little group of q
        U_scc = []
        for s in self.symmetries.S_s:
            U_cc, sign = self.symmetries.get_symmetry_operator(s)
            U_scc.append(sign * U_cc)
        U_scc = np.array(U_scc)

        # Determine the irreducible BZ
        bzk_kc, ibzk_kc, _ = get_reduced_bz(self.qpd.gd.cell_cv,
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

    def get_tetrahedron_kpt_domain(self, *, pbc_c):
        ik_kc = self.get_tetrahedron_ikpts(pbc_c=pbc_c)
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

    @timer('symmetrize_wGG')
    def symmetrize_wGG(self, A_wGG):
        """Symmetrize an array in GG'."""

        for A_GG in A_wGG:
            tmp_GG = np.zeros_like(A_GG, order='C')
            # tmp2_GG = np.zeros_like(A_GG)

            for s in self.symmetries.S_s:
                G_G = self.G_sG[s]
                _, sign = self.symmetries.get_symmetry_operator(s)
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

    @timer('symmetrize_wxvG')
    def symmetrize_wxvG(self, A_wxvG):
        """Symmetrize chi0_wxvG"""
        A_cv = self.qpd.gd.cell_cv
        iA_cv = self.qpd.gd.icell_cv

        if self.use_time_reversal:
            # ::-1 corresponds to transpose in wing indices
            AT_wxvG = A_wxvG[:, ::-1]

        tmp_wxvG = np.zeros_like(A_wxvG)
        for s in self.symmetries.S_s:
            G_G = self.G_sG[s]
            U_cc, sign = self.symmetries.get_symmetry_operator(s)
            M_vv = np.dot(np.dot(A_cv.T, U_cc.T), iA_cv)
            if sign == 1:
                tmp = sign * np.dot(M_vv.T, A_wxvG[..., G_G])
            elif sign == -1:
                tmp = sign * np.dot(M_vv.T, AT_wxvG[..., G_G])
            tmp_wxvG += np.transpose(tmp, (1, 2, 0, 3))

        # Overwrite the input
        A_wxvG[:] = tmp_wxvG / len(self.symmetries)

    @timer('symmetrize_wvv')
    def symmetrize_wvv(self, A_wvv):
        """Symmetrize chi_wvv."""
        A_cv = self.qpd.gd.cell_cv
        iA_cv = self.qpd.gd.icell_cv
        tmp_wvv = np.zeros_like(A_wvv)
        if self.use_time_reversal:
            AT_wvv = np.transpose(A_wvv, (0, 2, 1))

        for s in self.symmetries.S_s:
            U_cc, sign = self.symmetries.get_symmetry_operator(s)
            M_vv = np.dot(np.dot(A_cv.T, U_cc.T), iA_cv)
            if sign == 1:
                tmp = np.dot(np.dot(M_vv.T, A_wvv), M_vv)
            elif sign == -1:
                tmp = np.dot(np.dot(M_vv.T, AT_wvv), M_vv)
            tmp_wvv += np.transpose(tmp, (1, 0, 2))

        # Overwrite the input
        A_wvv[:] = tmp_wvv / len(self.symmetries)

    def initialize_G_maps(self):
        """Calculate the Gvector mappings."""
        qpd = self.qpd
        B_cv = 2.0 * np.pi * qpd.gd.icell_cv
        G_Gv = qpd.get_reciprocal_vectors(add_q=False)
        G_Gc = np.dot(G_Gv, np.linalg.inv(B_cv))
        Q_G = qpd.Q_qG[0]

        G_SG = [None] * self.nsym
        for S in self.symmetries.S_s:
            U_cc, sign = self.symmetries.get_symmetry_operator(S)
            iU_cc = np.linalg.inv(U_cc).T
            UG_Gc = np.dot(G_Gc - self.symmetries.shift_Sc[S], sign * iU_cc)

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
            G_SG[S] = np.array(G_G, dtype=np.int32)
        return G_SG

    def unfold_ibz_kpoint(self, ik):
        """Return kpoints related to irreducible kpoint."""
        kd = self.kd
        K_k = np.unique(kd.bz2bz_ks[kd.ibz2bz_k[ik]])
        K_k = K_k[K_k != -1]
        return K_k
