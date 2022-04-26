"""
A class for finding optimal
orbitals of the KS-DFT or PZ-SIC
functionals using exponential transformation
direct minimization in lcao mode

arXiv:2101.12597 [physics.comp-ph]
Comput. Phys. Commun. 267, 108047 (2021).
https://doi.org/10.1016/j.cpc.2021.108047
"""

import numpy as np
from copy import deepcopy
from gpaw.lcao.eigensolver import DirectLCAO
from gpaw.utilities.tools import tri2full


class DirectMinLCAO(DirectLCAO):

    def __init__(self, wfs, ham, nkpts, diagonalizer=None,
                 orthonormalization='gramschmidt',
                 need_init_orbs=True,
                 constraints=[]):

        super(DirectMinLCAO, self).__init__(diagonalizer)
        super(DirectMinLCAO, self).initialize(wfs.gd, wfs.dtype,
                                              wfs.setups.nao, wfs.ksl)
        self.orthonormalization = orthonormalization
        self.need_init_orbs = need_init_orbs
        self.nkpts = nkpts
        self.reference_orbitals = {}
        self.initialize_orbitals(wfs, ham)
        self.constraints = constraints

    def __repr__(self):
        pass

    def set_reference_orbitals(self, wfs, n_dim):
        for kpt in wfs.kpt_u:
            u = self.kpointval(kpt)
            self.reference_orbitals[u] = np.copy(kpt.C_nM[:n_dim[u]])

    def appy_transformation_kpt(self, wfs, u_mat, kpt, c_ref=None,
                                broadcast=True,
                                update_proj=True):
        """
        If c_ref are not provided then
        kpt.C_nM <- u_mat kpt.C_nM
        otherwise kpt.C_nM <- u_mat c_ref
        """

        dimens1 = u_mat.shape[1]
        dimens2 = u_mat.shape[0]

        if c_ref is None:
            kpt.C_nM[:dimens2] = u_mat @ kpt.C_nM[:dimens1]
        else:
            kpt.C_nM[:dimens2] = u_mat @ c_ref[:dimens1]

        if broadcast:
            with wfs.timer('Broadcast coefficients'):
                wfs.gd.comm.broadcast(kpt.C_nM, 0)
        if update_proj:
            with wfs.timer('Calculate projections'):
                wfs.atomic_correction.calculate_projections(wfs, kpt)

    def initialize_orbitals(self, wfs, ham):

        """
        If it is the first use of the scf then initialize
        coefficient matrix using eigensolver
        and then localise orbitals
        """

        # if it is the first use of the scf then initialize
        # coefficient matrix using eigensolver
        orthname = self.orthonormalization
        need_canon_coef = \
            (not wfs.coefficients_read_from_file and self.need_init_orbs)
        if need_canon_coef or orthname == 'diag':
            super(DirectMinLCAO, self).iterate(ham, wfs)
        else:
            wfs.orthonormalize(type=orthname)
        wfs.coefficients_read_from_file = False
        self.need_init_orbs = False

    def calc_grad(self, wfs, ham, kpt, func, evecs, evals, matrix_exp,
                  representation, ind_up, constraints):

        """
        Gradient w.r.t. skew-Hermitian matrices
        """

        h_mm = self.calculate_hamiltonian_matrix(ham, wfs, kpt)
        # make matrix hermitian
        tri2full(h_mm)
        # calc gradient and eigenstate error
        g_mat, error = func.get_gradients(
            h_mm, kpt.C_nM, kpt.f_n, evecs, evals,
            kpt, wfs, wfs.timer, matrix_exp,
            representation, ind_up, constraints)

        return g_mat, error

    def update_to_canonical_orbitals(self, wfs, ham, kpt,
                                     update_ref_orbs_canonical, restart):
        """
        Choose canonical orbitals
        """

        h_mm = self.calculate_hamiltonian_matrix(ham, wfs, kpt)
        tri2full(h_mm)
        if update_ref_orbs_canonical or restart:
            # Diagonalize entire Hamiltonian matrix
            with wfs.timer('Diagonalize and rotate'):
                kpt.C_nM, kpt.eps_n = rotate_subspace(h_mm, kpt.C_nM)
        else:
            # Diagonalize equally occupied subspaces separately
            n_init = 0
            while True:
                n_fin = \
                    find_equally_occupied_subspace(kpt.f_n, n_init)
                with wfs.timer('Diagonalize and rotate'):
                    kpt.C_nM[n_init:n_init + n_fin, :], \
                        kpt.eps_n[n_init:n_init + n_fin] = \
                        rotate_subspace(
                            h_mm, kpt.C_nM[n_init:n_init + n_fin, :])
                n_init += n_fin
                if n_init == len(kpt.f_n):
                    break
                elif n_init > len(kpt.f_n):
                    raise RuntimeError('Bug is here!')

        with wfs.timer('Calculate projections'):
            self.update_projections(wfs, kpt)

    def sort_orbitals(self, wfs, kpt, ind):
        """
        sort orbitals according to indices stored in ind
        """
        kpt.C_nM[np.arange(len(ind)), :] = kpt.C_nM[ind, :]
        self.update_projections(wfs, kpt)

    def update_projections(self, wfs, kpt):
        """
        calculate projections kpt.P_ani
        """

        wfs.atomic_correction.calculate_projections(wfs, kpt)

    def orbital_energies(self, wfs, ham, kpt):
        """
        diagonal elements of hamiltonian matrix in orbital representation
        """

        h_mm = self.calculate_hamiltonian_matrix(ham, wfs, kpt)
        tri2full(h_mm)
        # you probably need only diagonal terms?
        # wouldn't "for" be faster?
        h_mm = kpt.C_nM.conj() @ h_mm.conj() @ kpt.C_nM.T

        return h_mm.diagonal().real.copy()

    def kpointval(self, kpt):
        return self.nkpts * kpt.s + kpt.q


def rotate_subspace(h_mm, c_nm):
    """
    choose canonical orbitals
    """
    l_nn = (c_nm @ h_mm @ c_nm.conj().T).conj()
    # check if diagonal then don't rotate? it could save a bit of time
    eps, w = np.linalg.eigh(l_nn)
    return w.T.conj() @ c_nm, eps


def find_equally_occupied_subspace(f_n, index=0):
    return np.searchsorted(f_n[index] - f_n[index:], 1.0e-10)
