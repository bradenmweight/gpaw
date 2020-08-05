"""Tetrahedron method for Brillouin-zone integrations.

See::

    Improved tetrahedron method for Brillouin-zone integrations.

    Peter E. Blöchl, O. Jepsen, and O. K. Andersen
    Phys. Rev. B 49, 16223 – Published 15 June 1994

    DOI:https://doi.org/10.1103/PhysRevB.49.16223
"""

from math import nan
from typing import List, Tuple, Any
import numpy as np
from scipy.spatial import Delaunay

from gpaw.occupations import (ZeroWidth, findroot, collect_eigelvalues,
                              distribute_occupation_numbers,
                              OccupationNumberCalculator, ParallelLayout)
from gpaw.mpi import broadcast_float

Array = Any


def bja1(e1: Array, e2: Array, e3: Array, e4: Array) -> Tuple[float, float]:
    """Eq. (A2) and (C2) from Blöchl, Jepsen and Andersen."""
    x = 1.0 / ((e2 - e1) * (e3 - e1) * (e4 - e1))
    return (-(e1**3).dot(x),
            3 * (e1**2).dot(x))


def bja2(e1: Array, e2: Array, e3: Array, e4: Array) -> Tuple[float, float]:
    """Eq. (A3) and (C3) from Blöchl, Jepsen and Andersen."""
    x = 1.0 / ((e3 - e1) * (e4 - e1))
    y = (e3 - e1 + e4 - e2) / ((e3 - e2) * (e4 - e2))
    return (x.dot((e2 - e1)**2
                  - 3 * (e2 - e1) * e2
                  + 3 * e2**2
                  + y * e2**3),
            x.dot(3 * (e2 - e1)
                  - 6 * e2
                  - 3 * y * e2**2))


def bja3(e1: Array, e2: Array, e3: Array, e4: Array) -> Tuple[float, float]:
    """Eq. (A4) and (C4) from Blöchl, Jepsen and Andersen."""
    x = 1.0 / ((e4 - e1) * (e4 - e2) * (e4 - e3))
    return (len(e1) - x.dot(e4**3),
            3 * x.dot(e4**2))


def bja1b(e1: Array, e2: Array, e3: Array, e4: Array) -> Array:
    """Eq. (B2)-(B5) from Blöchl, Jepsen and Andersen."""
    C = -0.25 * e1**3 / ((e2 - e1) * (e3 - e1) * (e4 - e1))
    w2 = -C * e1 / (e2 - e1)
    w3 = -C * e1 / (e3 - e1)
    w4 = -C * e1 / (e4 - e1)
    w1 = 4 * C - w2 - w3 - w4
    return np.array([w1, w2, w3, w4])


def bja2b(e1: Array, e2: Array, e3: Array, e4: Array) -> Array:
    """Eq. (B7)-(B10) from Blöchl, Jepsen and Andersen."""
    C1 = 0.25 * e1**2 / ((e4 - e1) * (e3 - e1))
    C2 = 0.25 * e1 * e2 * e3 / ((e4 - e1) * (e3 - e2) * (e3 - e1))
    C3 = 0.25 * e2**2 * e4 / ((e4 - e2) * (e3 - e2) * (e4 - e1))
    w1 = C1 + (C1 + C2) * e3 / (e3 - e1) + (C1 + C2 + C3) * e4 / (e4 - e1)
    w2 = C1 + C2 + C3 + (C2 + C3) * e3 / (e3 - e2) + C3 * e4 / (e4 - e2)
    w3 = (C1 + C2) * e1 / (e1 - e3) - (C2 + C3) * e2 / (e3 - e2)
    w4 = (C1 + C2 + C3) * e1 / (e1 - e4) + C3 * e2 / (e2 - e4)
    return np.array([w1, w2, w3, w4])


def bja3b(e1: Array, e2: Array, e3: Array, e4: Array) -> Array:
    """Eq. (B14)-(B17) from Blöchl, Jepsen and Andersen."""
    C = 0.25 * e4**3 / ((e4 - e1) * (e4 - e2) * (e4 - e3))
    w1 = 0.25 - C * e4 / (e4 - e1)
    w2 = 0.25 - C * e4 / (e4 - e2)
    w3 = 0.25 - C * e4 / (e4 - e3)
    w4 = 1.0 - 4 * C - w1 - w2 - w3
    return np.array([w1, w2, w3, w4])


def triangulate_submesh(rcell_cv: Array) -> Array:
    """Find the 6 tetrahedra."""
    ABC_sc = np.array([[A, B, C]
                       for A in [0, 1] for B in [0, 1] for C in [0, 1]])
    dt = Delaunay(ABC_sc.dot(rcell_cv))
    s_tq = dt.simplices
    assert s_tq.shape == (6, 4)
    ABC_tqc = ABC_sc[s_tq]
    return ABC_tqc


def triangulate_everything(size_c: Array,
                           ABC_tqc: Array,
                           i_k: Array) -> Array:
    nbzk = size_c.prod()
    ABC_ck = np.unravel_index(np.arange(nbzk), size_c)
    ABC_tqck = ABC_tqc[..., np.newaxis] + ABC_ck
    ABC_cktq = np.transpose(ABC_tqck, (2, 3, 0, 1))
    k_ktq = np.ravel_multi_index(ABC_cktq.reshape((3, nbzk * 6 * 4)),
                                 size_c,
                                 mode='wrap').reshape((nbzk, 6, 4))
    i_ktq = i_k[k_ktq]
    return i_ktq


class TetrahedronMethod(OccupationNumberCalculator):
    name = 'tetrahedron-method'

    def __init__(self,
                 rcell: List[List[float]],
                 size: Tuple[int, int, int],
                 bz2ibzmap: List[int] = None,
                 parallel_layout: ParallelLayout = None):
        """Tetrahedron method for calculating occupation numbers.

        The reciprocal cell (rcell) can be given in arbitrary units
        (only the shape matters) and size is the size of the
        Monkhorst-Pack grid.  If k-points have been symmetry-reduced
        the bz2ibzmap mapping BZ k-point indizes to IBZ -k-point indices.
        """

        OccupationNumberCalculator.__init__(self, parallel_layout)

        self.rcell_cv = np.asarray(rcell)
        self.size_c = np.asarray(size)

        nbzk = self.size_c.prod()

        if bz2ibzmap is None:
            bz2ibzmap = np.arange(nbzk)

        self.i_k = np.asarray(bz2ibzmap)

        assert self.size_c.shape == (3,)
        assert self.rcell_cv.shape == (3, 3)
        assert self.i_k.shape == (nbzk,)

        ABC_tqc = triangulate_submesh(
            self.rcell_cv / self.size_c[:, np.newaxis])

        self.i_ktq = triangulate_everything(self.size_c, ABC_tqc, self.i_k)

    def _calculate(self,
                   nelectrons,
                   eig_qn,
                   weight_q,
                   f_qn,
                   fermi_level_guess=nan):

        if np.isnan(fermi_level_guess):
            zero = ZeroWidth()
            fermi_level_guess, _ = zero._calculate(
                nelectrons, eig_qn, weight_q, f_qn)
            if np.isinf(fermi_level_guess):
                return fermi_level_guess, 0.0

        x = fermi_level_guess

        eig_in, weight_i, nkpts_r = collect_eigelvalues(eig_qn, weight_q,
                                                        self.bd, self.kpt_comm)

        assert self.i_k.max() == len(eig_in) - 1

        if eig_in is not None:
            def func(x, eig_in=eig_in):
                n, dn = count(x, eig_in, self.i_ktq)
                return n - nelectrons, dn

            fermi_level, niter = findroot(func, x)
            f_in = weights(eig_in - fermi_level, self.i_ktq)
            f_in *= 1 / (weight_i[:, np.newaxis] * len(self.i_k))
        else:
            f_in = None
            fermi_level = nan

        distribute_occupation_numbers(f_in, f_qn, nkpts_r,
                                      self.bd, self.kpt_comm)

        if self.kpt_comm.rank == 0:
            fermi_level = broadcast_float(fermi_level, self.bd.comm)
        fermi_level = broadcast_float(fermi_level, self.kpt_comm)

        return fermi_level, 0.0


def count(fermi_level: float,
          eig_in: Array,
          i_ktq: Array) -> Tuple[float, float]:
    """Count electrons.

    Return number of electrons and derivative with respect to fermi level.
    """
    eig_in = eig_in - fermi_level
    nocc_i = (eig_in < 0.0).sum(axis=1)
    n1 = nocc_i.min()
    n2 = nocc_i.max()

    ne = n1
    dnedef = 0.0

    if n1 == n2:
        return ne, dnedef

    ntetra = 6 * i_ktq.shape[0]
    eig_Tq = eig_in[i_ktq, n1:n2].transpose((0, 1, 3, 2)).reshape(
        (ntetra * (n2 - n1), 4))
    eig_Tq.sort(axis=1)

    eig_Tq = eig_Tq[eig_Tq[:, 0] < 0.0]

    mask1_T = eig_Tq[:, 1] > 0.0
    mask2_T = ~mask1_T & (eig_Tq[:, 2] > 0.0)
    mask3_T = ~mask1_T & ~mask2_T & (eig_Tq[:, 3] > 0.0)

    for mask_T, bjaa in [(mask1_T, bja1), (mask2_T, bja2), (mask3_T, bja3)]:
        n, dn = bjaa(*eig_Tq[mask_T].T)
        ne += n / ntetra
        dnedef += dn / ntetra

    mask4_T = ~mask1_T & ~mask2_T & ~mask3_T
    ne += mask4_T.sum() / ntetra

    return ne, dnedef


def weights(eig_in: Array, i_ktq: Array) -> Array:
    """Calculate occupation numbers."""
    nocc_i = (eig_in < 0.0).sum(axis=1)
    n1 = nocc_i.min()
    n2 = nocc_i.max()

    f_in = np.zeros_like(eig_in)

    for i in i_ktq[:, 0, 0]:
        f_in[i, :n1] += 6.0

    if n1 == n2:
        return

    ntetra = 6 * i_ktq.shape[0]
    eig_Tq = eig_in[i_ktq, n1:n2].transpose((0, 1, 3, 2)).reshape(
        (ntetra * (n2 - n1), 4))
    q_Tq = eig_Tq.argsort(axis=1)
    eig_Tq = np.take_along_axis(eig_Tq, q_Tq, 1)
    f_Tq = np.zeros_like(eig_Tq)

    mask0_T = eig_Tq[:, 0] > 0.0
    mask1_T = ~mask0_T & (eig_Tq[:, 1] > 0.0)
    mask2_T = ~mask0_T & ~mask1_T & (eig_Tq[:, 2] > 0.0)
    mask3_T = ~mask0_T & ~mask1_T & ~mask2_T & (eig_Tq[:, 3] > 0.0)

    for mask_T, bja in [(mask1_T, bja1b), (mask2_T, bja2b), (mask3_T, bja3b)]:
        w_qT = bja(*eig_Tq[mask_T].T)
        f_Tq[mask_T] += w_qT.T

    mask4_T = ~mask0_T & ~mask1_T & ~mask2_T & ~mask3_T
    f_Tq[mask4_T] += 0.25

    ktn_T = np.array(np.unravel_index(np.arange(len(eig_Tq)),
                                      (len(i_ktq), 6, n2 - n1))).T
    for f_q, q_q, (k, t, n) in zip(f_Tq, q_Tq, ktn_T):
        for q, f in zip(q_q, f_q):
            f_in[i_ktq[k, t, q], n1 + n] += f

    f_in *= 1 / 6

    return f_in
