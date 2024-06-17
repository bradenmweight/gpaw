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
class SymmetryAnalyzer:
    point_group: bool = True
    time_reversal: bool = True

    def analyze(self, kpoints, qpd, context):
        return PWSymmetryAnalyzer(
            kpoints, qpd, context, not self.point_group,
            not self.time_reversal)


class PWSymmetryAnalyzer:
    """Class for handling planewave symmetries."""

    def __init__(self, kpoints, qpd, context,
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
        U_scc = kd.symmetry.op_scc
        self.nU = len(U_scc)

        self.nsym = 2 * self.nU
        self.use_time_reversal = not self.disable_time_reversal

        self.kptfinder = kpoints.kptfinder

        self.s_s, self.shift_sc = self.analyze_symmetries()
        self.G_sG = self.initialize_G_maps()

        self.context.print(self.get_infostring())
        self.context.print(self.symmetry_description())

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

        txt += f'Found {len(self.s_s)} allowed symmetries. '

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
        ns = len(self.s_s)
        y = 0
        for y in range((ns + nx - 1) // nx):
            for c in range(3):
                tisl = []
                for x in range(nx):
                    s = x + y * nx
                    if s == ns:
                        break
                    op_cc, sign = self.get_symmetry_operator(self.s_s[s])
                    op_c = sign * op_cc[c]
                    tisl.append(f'  ({op_c[0]:2d} {op_c[1]:2d} {op_c[2]:2d})')
                tisl.append('\n')
                isl.append(''.join(tisl))
            isl.append('\n')
        return ''.join(isl)

    @timer('Analyze symmetries.')
    def analyze_symmetries(self):
        r"""Determine allowed symmetries.

        An direct symmetry U must fulfill::

          U \mathbf{q} = q + \Delta

        Under time-reversal (indirect) it must fulfill::

          -U \mathbf{q} = q + \Delta

        where :math:`\Delta` is a reciprocal lattice vector.
        """
        qpd = self.qpd

        # Shortcuts
        q_c = qpd.q_c
        kd = self.kd

        U_scc = kd.symmetry.op_scc
        nU = self.nU
        nsym = self.nsym

        shift_sc = np.zeros((nsym, 3), int)
        conserveq_s = np.zeros(nsym, bool)

        newq_sc = np.dot(U_scc, q_c)

        # Direct symmetries
        dshift_sc = (newq_sc - q_c[np.newaxis]).round().astype(int)
        inds_s = np.argwhere((newq_sc == q_c[np.newaxis] + dshift_sc).all(1))
        conserveq_s[inds_s] = True

        shift_sc[:nU] = dshift_sc

        # Time reversal
        trshift_sc = (-newq_sc - q_c[np.newaxis]).round().astype(int)
        trinds_s = np.argwhere((-newq_sc == q_c[np.newaxis] +
                                trshift_sc).all(1)) + nU
        conserveq_s[trinds_s] = True
        shift_sc[nU:nsym] = trshift_sc

        # The indices of the allowed symmetries
        s_s = conserveq_s.nonzero()[0]

        # Filter out disabled symmetries
        if self.disable_point_group:
            s_s = list(filter(self.is_not_point_group, s_s))

        if self.disable_time_reversal:
            s_s = list(filter(self.is_not_time_reversal, s_s))

        # You are not allowed to use non-symmorphic syms, sorry. So we remove
        # the option and always filter those symmetries out.
        s_s = list(filter(self.is_not_non_symmorphic, s_s))

#        stmp_s = []
#        for s in s_s:
#            if self.kd.bz2bz_ks[0, s] == -1:
#                assert (self.kd.bz2bz_ks[:, s] == -1).all()
#            else:
#                stmp_s.append(s)

#        s_s = stmp_s

        return s_s, shift_sc

    def is_not_point_group(self, s):
        U_scc = self.kd.symmetry.op_scc
        nU = self.nU
        return (U_scc[s % nU] == np.eye(3)).all()

    def is_not_time_reversal(self, s):
        nU = self.nU
        return not bool(s // nU)

    def is_not_non_symmorphic(self, s):
        ft_sc = self.kd.symmetry.ft_sc
        nU = self.nU
        return not bool(ft_sc[s % nU].any())

    def how_many_symmetries(self):
        """Return number of symmetries."""
        return len(self.s_s)

    @timer('Group kpoints')
    def group_kpoints(self, K_k=None):
        """Group kpoints according to the reduced symmetries"""
        if K_k is None:
            K_k = np.arange(self.kd.nbzkpts)
        s_s = self.s_s
        bz2bz_ks = self.kd.bz2bz_ks
        nk = len(bz2bz_ks)
        sbz2sbz_ks = bz2bz_ks[K_k][:, s_s]  # Reduced number of symmetries
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
        for s in self.s_s:
            U_cc, sign = self.get_symmetry_operator(s)
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

            for s in self.s_s:
                G_G = self.G_sG[s]
                _, sign = self.get_symmetry_operator(s)
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
            A_GG[:] = tmp_GG / self.how_many_symmetries()

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
        for s in self.s_s:
            G_G = self.G_sG[s]
            U_cc, sign = self.get_symmetry_operator(s)
            M_vv = np.dot(np.dot(A_cv.T, U_cc.T), iA_cv)
            if sign == 1:
                tmp = sign * np.dot(M_vv.T, A_wxvG[..., G_G])
            elif sign == -1:
                tmp = sign * np.dot(M_vv.T, AT_wxvG[..., G_G])
            tmp_wxvG += np.transpose(tmp, (1, 2, 0, 3))

        # Overwrite the input
        A_wxvG[:] = tmp_wxvG / self.how_many_symmetries()

    @timer('symmetrize_wvv')
    def symmetrize_wvv(self, A_wvv):
        """Symmetrize chi_wvv."""
        A_cv = self.qpd.gd.cell_cv
        iA_cv = self.qpd.gd.icell_cv
        tmp_wvv = np.zeros_like(A_wvv)
        if self.use_time_reversal:
            AT_wvv = np.transpose(A_wvv, (0, 2, 1))

        for s in self.s_s:
            U_cc, sign = self.get_symmetry_operator(s)
            M_vv = np.dot(np.dot(A_cv.T, U_cc.T), iA_cv)
            if sign == 1:
                tmp = np.dot(np.dot(M_vv.T, A_wvv), M_vv)
            elif sign == -1:
                tmp = np.dot(np.dot(M_vv.T, AT_wvv), M_vv)
            tmp_wvv += np.transpose(tmp, (1, 0, 2))

        # Overwrite the input
        A_wvv[:] = tmp_wvv / self.how_many_symmetries()

    def timereversal(self, s):
        """Is this a time-reversal symmetry?"""
        tr = bool(s // self.nU)
        return tr

    def get_symmetry_operator(self, s):
        """Return symmetry operator s."""
        U_scc = self.kd.symmetry.op_scc

        reds = s % self.nU
        if self.timereversal(s):
            sign = -1
        else:
            sign = 1

        return U_scc[reds], sign

    def initialize_G_maps(self):
        """Calculate the Gvector mappings."""
        qpd = self.qpd
        B_cv = 2.0 * np.pi * qpd.gd.icell_cv
        G_Gv = qpd.get_reciprocal_vectors(add_q=False)
        G_Gc = np.dot(G_Gv, np.linalg.inv(B_cv))
        Q_G = qpd.Q_qG[0]

        G_sG = [None] * self.nsym
        for s in self.s_s:
            U_cc, sign = self.get_symmetry_operator(s)
            iU_cc = np.linalg.inv(U_cc).T
            UG_Gc = np.dot(G_Gc - self.shift_sc[s], sign * iU_cc)

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
            G_sG[s] = np.array(G_G, dtype=np.int32)
        return G_sG

    def unfold_ibz_kpoint(self, ik):
        """Return kpoints related to irreducible kpoint."""
        kd = self.kd
        K_k = np.unique(kd.bz2bz_ks[kd.ibz2bz_k[ik]])
        K_k = K_k[K_k != -1]
        return K_k
