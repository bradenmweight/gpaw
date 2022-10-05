import stat
from collections import defaultdict
from curses import wrapper
from functools import wraps
from typing import List, Tuple, Union

import numpy as np
import numpy.typing as npt
from ase.units import Hartree
from gpaw import GPAW
from gpaw.lcao.pwf2 import LCAOwrap
from gpaw.lcao.tightbinding import TightBinding  # as LCAOTightBinding
from gpaw.lcao.tools import get_bfi, get_lcao_hamiltonian
from gpaw.utilities.blas import r2k
from gpaw.utilities.tools import lowdin, tri2full
from scipy.linalg import eigh

# Numeric type
Numeric = Union[int, float, complex]


def get_subspace(A_MM, index):
    """Get the subspace spanned by the basis function listed in index."""
    assert A_MM.ndim == 2 and A_MM.shape[0] == A_MM.shape[1]
    return A_MM.take(index, 0).take(index, 1)


def get_orthonormal_subspace(H_MM, S_MM, index=None):
    """Get orthonormal eigen-values and -vectors of subspace listed in index."""
    if index is not None:
        h_ww = get_subspace(H_MM, index)
        s_ww = get_subspace(S_MM, index)
    else:
        h_ww = H_MM
        s_ww = S_MM
    eps, v = eigh(h_ww, s_ww)
    return eps, v


def subdiagonalize(H_MM, S_MM, block_lists):
    """Subdiagonalize blocks."""
    nM = len(H_MM)
    v_MM = np.eye(nM)
    eps_M = np.zeros(nM)
    mask_M = np.ones(nM, dtype=int)
    for block in block_lists:
        eps, v = get_orthonormal_subspace(H_MM, S_MM, block)
        v_MM[np.ix_(block, block)] = v
        eps_M[block] = eps
        mask_M[block] = 0
    epsx_M = np.ma.masked_array(eps_M, mask=mask_M)
    return epsx_M, v_MM


def subdiagonalize_atoms(calc, H_MM, S_MM, atom_list=None):
    """Subdiagonalize atomic sub-spaces."""
    if atom_list is None:
        atom_list = range(len(calc.atoms))
    if isinstance(atom_list, int):
        atom_list = [atom_list]
    block_lists = []
    for a in atom_list:
        M = calc.wfs.basis_functions.M_a[a]
        block = range(M, M + calc.wfs.setups[a].nao)
        block_lists.append(block)
    return subdiagonalize(H_MM, S_MM, block_lists)


def get_orbitals(calc, U_Mw: npt.NDArray, q=0):
    """Get orbitals from AOs coefficients.

    Parameters
    ----------
    calc : GPAW
        LCAO calculator
    U_Mw : array_like
        LCAO expansion coefficients.
    """
    Nw = U_Mw.shape[1]
    C_wM = np.ascontiguousarray(U_Mw.T).astype(calc.wfs.dtype)
    w_wG = calc.wfs.gd.zeros(Nw, dtype=calc.wfs.dtype)
    calc.wfs.basis_functions.lcao_to_grid(C_wM, w_wG, q=q)
    return w_wG


def get_xc(calc, v_wG, P_awi=None):
    """Get exchange-correlation part of the Hamiltonian."""
    if calc.density.nt_sg is None:
        calc.density.interpolate_pseudo_density()
    nt_sg = calc.density.nt_sg
    vxct_sg = calc.density.finegd.zeros(calc.wfs.nspins)
    calc.hamiltonian.xc.calculate(calc.density.finegd, nt_sg, vxct_sg)
    vxct_G = calc.wfs.gd.empty()
    calc.hamiltonian.restrict_and_collect(vxct_sg[0], vxct_G)

    # Integrate pseudo part
    Nw = len(v_wG)
    xc_ww = np.empty((Nw, Nw))
    r2k(0.5 * calc.wfs.gd.dv, v_wG, vxct_G * v_wG, 0.0, xc_ww)
    tri2full(xc_ww, "L")

    # Atomic PAW corrections required? XXX
    if P_awi is not None:
        raise NotImplementedError("""Atomic PAW corrections not included.
                                  Have a look at pwf2::get_xc2 for inspiration.""")

    return xc_ww * Hartree


def get_Fcore():
    pass


def dots(*args):
    """Multi dot."""
    x = args[0]
    for M in args[1:]:
        x = np.dot(x, M)
    return x


def rotate_matrix(M, U):
    """Rotate a matrix."""
    return dots(U.T.conj(), M, U)


class BasisTransform:
    """Class to perform a basis transformation.

    Attributes
    ----------
    U_MM : np.ndarray
        Rotation matrix between 2 basis.
    indices : array_like, optional
        Subset of indices in the new basis.
    U_Mw : np.ndarray
        Same as `U_MM` but the include only `indices` of the new basis.

    Methods
    -------
    rotate_matrx(A_MM, keep_rest=False)
        Rotate a matrix.
    rotate_projections(P_aMi, keep_rest=False)
        Rotate PAW atomic projects.
    rotate_function(P_aMi, keep_rest=False)
        Rotate PAW atomic projects.

    """

    def __init__(self, U_MM: npt.NDArray, indices: npt.ArrayLike = None) -> None:
        self.U_MM = U_MM
        self.indices = indices
        if indices is not None:
            self.U_Mw = np.ascontiguousarray(U_MM[:, indices])
        else:
            self.U_Mw = None

    def get_rotation(self, keep_rest=False):
        if keep_rest or self.U_Mw is None:
            return self.U_MM
        return self.U_Mw

    def rotate_matrix(self, A_MM, keep_rest=False):
        U_Mx = self.get_rotation(keep_rest)
        return dots(U_Mx.T.conj(), A_MM, U_Mx)

    def rotate_projections(self, P_aMi, keep_rest=False):
        U_Mx = self.get_rotation(keep_rest)
        P_awi = {}
        for a, P_Mi in P_aMi.items():
            P_awi[a] = np.tensordot(U_Mx, P_Mi, axes=[[0], [0]])
        return P_awi

    def rotate_function(self, Psi_MG, keep_rest=False):
        U_Mx = self.get_rotation(keep_rest)
        return np.tensordot(U_Mx, Psi_MG, axes=[[0], [0]])


class EffectiveModel(BasisTransform):
    """Class for an effective model."""

    def __init__(self, U_MM: npt.NDArray, indices: npt.ArrayLike, S_MM: npt.NDArray = None) -> None:
        """

        Parameters:
        ----------
        S_MM : (M,M) array_like, optional
            LCAO overlap matrix. If provided, the resulting basis
            is orthogonalized.
        """
        if S_MM is not None:
            lowdin(self.U_Mw, self.rotate_matrix(S_MM))
            np.testing.assert_allclose(self.rotate_matrix(
                S_MM), np.eye(self.U_Mw.shape[1]))
            U_MM = U_MM[:]
            U_MM = U_MM[:, indices] = self.U_Mw

        super().__init__(U_MM, indices)

    def get_static_correction(self, H_MM: npt.NDArray, S_MM: npt.NDArray, z: Numeric = 0. + 1e-5j):
        """Get static correction to model Hamiltonian."""
        w = self.indices  # Alias

        Hp_MM = self.rotate_matrix(H_MM, keep_rest=True)
        Sp_MM = self.rotate_matrix(S_MM, keep_rest=True)
        Up_Mw = Sp_MM[:, w].dot(np.linalg.inv(Sp_MM[np.ix_(w, w)]))

        H_ww = self.rotate_matrix(H_MM)
        S_ww = self.rotate_matrix(S_MM)

        # Coupled
        G = np.linalg.inv(z * Sp_MM - Hp_MM)
        G_inv = np.linalg.inv(rotate_matrix(G, Up_Mw))
        # Uncoupled
        G0_inv = z * S_ww - H_ww
        # Hybridization
        D0 = G0_inv - G_inv
        return D0.real


class Subdiagonalization(BasisTransform):
    """Class to perform a subdiagonalization of the Hamiltonian.

    Attributes
    ----------
    H_MM, S_MMM : np.ndarray
        LCAO Hamiltonian and overlap matrices.
    U_MM : np.ndarray
        Rotation matrix that subdiagonalizes the LCAO Hamiltonian.

    Methods:
    --------
    get_effective_model(index_list, ortho=None)
        Builds and effective model.

    """

    def __init__(self, H_MM: npt.NDArray, S_MM: npt.NDArray, blocks: List[List]) -> None:
        """
        Parameters
        ----------
        H_MM : (M,M) NDArray
            LCAO Hamiltonian.
        S_MM : (M,M) NDArray
            LCAO Overlap.
        blocks : list of lists
            Subdiagonalization blocks given as a list of list indices.
        """
        self.blocks = blocks
        self.H_MM = H_MM
        self.S_MM = S_MM
        self.eps_M, U_MM = subdiagonalize(
            self.H_MM, self.S_MM, blocks)
        super().__init__(U_MM)

    def get_model(self, indices, ortho=False):
        """Get a local orbital model for the subdiagonalized space."""
        return EffectiveModel(self.U_MM, indices, self.S_MM if ortho else None)


class LocalOrbitals(TightBinding):
    """Local Orbitals.

    Attributes
    ----------
    subdiag : Subdiagonalization
    model : EffectiveModel

    H_NMM, S_NMM: np.ndarray
        LCAO Hamiltonian and overlap matrices.
    N0: int
        Home unit cell index.

    Methods:
    --------
    set_effmodel(self, list_indices)
        Set the effective model from the a list of indices.
    h_and_s:
        Returns H_NMM in Hartree and S_NNM.
    get_hamiltonian(self)
        Returns H_NMM[N0]
    get_overlap(self)
        Returns S_NMM[N0]
    get_orbitals(self)
        Returns the orbitals on the grid.
    get_projections(self, q=0)
        Returns the rotated projections.
    get_xc(self)
        Returns the exchenge and correlation.
    get_Fcore(self)
        pass

    """

    def __init__(self, calc: GPAW):

        self.calc = calc
        self.gamma = calc.wfs.kd.gamma  # Gamma point calculation

        if self.gamma:
            self.calc = calc
            h = self.calc.hamiltonian
            wfs = self.calc.wfs
            kpt = wfs.kpt_u[0]

            H_MM = wfs.eigensolver.calculate_hamiltonian_matrix(h, wfs, kpt)
            S_MM = wfs.S_qMM[kpt.q]
            # XXX Converting to full matrices here
            tri2full(H_MM)
            tri2full(S_MM)
            self.H_NMM = H_MM[None, ...] * Hartree  # eV
            self.S_NMM = S_MM[None, ...]
            self.N0 = 0
        else:
            super().__init__(calc.atoms, calc)
            # Bloch to real
            self.H_NMM, self.S_NMM = TightBinding.h_and_s(self)
            self.H_NMM *= Hartree  # eV
            try:
                self.N0 = np.argwhere(
                    self.R_cN.T.dot(self.R_cN) < 1e-13).flat[0]
            except Exception as exc:
                raise RuntimeError(
                    "Must include central unit cell, i.e. R=[0,0,0].") from exc

    def subdiagonalize(self, symbols=None, blocks=None):
        """Subdiagonalize Hamiltonian and overlap matrices."""
        if symbols is not None:
            atoms = self.calc.atoms.symbols.search(symbols)
            blocks = [get_bfi(self.calc, [c]) for c in atoms]
        if blocks is None:
            raise RuntimeError("""User must provide either the element(s)
                               or a list of blocks to subdiagonalize.""")
        self.blocks = blocks
        self.subdiag = Subdiagonalization(
            self.H_NMM[self.N0], self.S_NMM[self.N0], blocks)

        # Group orbitals with same symmetry.
        eps = self.subdiag.eps_M.round(1)
        show = np.where(~eps.mask)[0]
        groups = defaultdict(list)
        for index in show:
            groups[eps[index]].append(index)

        self.groups = groups

    def take_model(self, indices=None, minimal=True, cutoff=1e-3, ortho=False):
        """Build an effective model."""
        if self.subdiag is None:
            raise RuntimeError("""Not yet subdiagonalized.""")

        eps = self.subdiag.eps_M.round(1)

        if indices is None:
            # Find active orbitals with energy closest to Fermi.

            fermi = round(self.calc.get_fermi_level(), 1)
            # diffs = [] # Min distance from fermi for each block
            indices = []  # Min distance index for each block
            for block in self.blocks:
                eb = eps[block]
                ib = np.abs(eb - fermi).argmin()
                indices.append(block[ib])
                # diffs.append(abs(eb[ib]))

        if not minimal:
            # Find orbitals that connect to active with a matrix element larger than cutoff

            # Look at gamma of 1st neighbor
            H_MM = self.H_NMM[(self.N0 + 1) % len(self.H_NMM)]
            H_MM = self.subdiag.rotate_matrix(H_MM)
            # H_MM = dots(self.subdiag.U_MM.T.conj(), H_MM, self.subdiag.U_MM)

            extend = []
            for group in self.groups.values():
                if np.isin(group, indices).any():
                    continue
                if np.abs(H_MM[np.ix_(indices, group)]).max() > cutoff:
                    extend += group

            # Expand model
            indices += extend

        self.indices = indices
        self.model = self.subdiag.get_model(indices)

        if self.gamma:
            H_Nww = self.model.rotate_matrix(self.H_NMM[0])[None, ...]
            S_Nww = self.model.rotate_matrix(self.S_NMM[0])[None, ...]

        else:
            # Bypass parent's LCAO construction.
            shape = (self.R_cN.shape[1],) + 2 * (len(self.indices),)
            dtype = self.H_NMM.dtype
            H_Nww = np.empty(shape, dtype)
            S_Nww = np.empty(shape, dtype)

            for N, (H_MM, S_MM) in enumerate(zip(self.H_NMM, self.S_NMM)):
                H_Nww[N] = self.model.rotate_matrix(H_MM)
                S_Nww[N] = self.model.rotate_matrix(S_MM)
        self.H_Nww = H_Nww
        self.S_Nww = S_Nww

        if minimal:
            # Add static correction of hybridization to minimal model.
            self.H_Nww[self.N0] += self.model.get_static_correction(
                self.H_NMM[self.N0], self.S_NMM[self.N0])

    def h_and_s(self):
        """Returns the Hamiltonian and overlap in Hartree units."""

        # Bypass TightBinding method.
        eV = 1 / Hartree
        return self.H_Nww * eV, self.S_Nww

    def get_hamiltonian(self):
        return self.H_Nww[self.N0]

    def get_overlap(self):
        return self.S_Nww[self.N0]

    def get_orbitals(self):
        return get_orbitals(self.calc, self.model.U_Mw)

    def get_projections(self, q=0):
        P_aMi = {a: P_aqMi[q] for a, P_aqMi in self.calc.wfs.P_aqMi.items()}
        return self.model.rotate_projections(P_aMi)

    def get_xc(self):
        return get_xc(self.calc, self.get_orbitals(), self.get_projections())

    def get_Fcore(self):
        pass
