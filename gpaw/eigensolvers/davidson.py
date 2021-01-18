from functools import partial

from ase.utils.timing import timer
from ase.parallel import world, parprint
import numpy as np

from gpaw import debug
from gpaw.eigensolvers.eigensolver import Eigensolver
from gpaw.matrix import matrix_matrix_multiply as mmm
from gpaw.blacs import BlacsGrid# , general_diagonalize_dc


class DummyArray:
    def __getitem__(self, x):
        return np.empty((0, 0))


class Davidson(Eigensolver):
    """Simple Davidson eigensolver

    It is expected that the trial wave functions are orthonormal
    and the integrals of projector functions and wave functions
    ``nucleus.P_uni`` are already calculated.

    Solution steps are:

    * Subspace diagonalization
    * Calculate all residuals
    * Add preconditioned residuals to the subspace and diagonalize
    """

    def __init__(self, niter=2, smin=None, normalize=True):
        Eigensolver.__init__(self)
        self.niter = niter
        self.smin = smin
        self.normalize = normalize
        self.i = 0


        if smin is not None:
            raise NotImplementedError(
                'See https://trac.fysik.dtu.dk/projects/gpaw/ticket/248')

        self.orthonormalization_required = False
        self.H_NN = DummyArray()
        self.S_NN = DummyArray()
        self.eps_N = DummyArray()

    def __repr__(self):
        return 'Davidson(niter=%d, smin=%r, normalize=%r)' % (
            self.niter, self.smin, self.normalize)

    def todict(self):
        return {'name': 'dav', 'niter': self.niter}

    def initialize(self, wfs):
        Eigensolver.initialize(self, wfs)

        if wfs.gd.comm.rank == 0 and wfs.bd.comm.rank == 0:
            # Allocate arrays
            B = self.nbands
            self.H_NN = np.zeros((2 * B, 2 * B), self.dtype)
            self.S_NN = np.zeros((2 * B, 2 * B), self.dtype)
            self.eps_N = np.zeros(2 * B)

    def estimate_memory(self, mem, wfs):
        Eigensolver.estimate_memory(self, mem, wfs)
        nbands = wfs.bd.nbands
        mem.subnode('H_nn', nbands * nbands * mem.itemsize[wfs.dtype])
        mem.subnode('S_nn', nbands * nbands * mem.itemsize[wfs.dtype])
        mem.subnode('H_2n2n', 4 * nbands * nbands * mem.itemsize[wfs.dtype])
        mem.subnode('S_2n2n', 4 * nbands * nbands * mem.itemsize[wfs.dtype])
        mem.subnode('eps_2n', 2 * nbands * mem.floatsize)

    @timer('Davidson')
    def iterate_one_k_point(self, ham, wfs, kpt, weights):
        """Do Davidson iterations for the kpoint"""
        bd = wfs.bd
        B = bd.nbands

        H_NN = self.H_NN
        S_NN = self.S_NN
        eps_N = self.eps_N

        slcomm, nrows, ncols, slsize = wfs.scalapack_parameters
        grid = BlacsGrid(slcomm, nrows, ncols)
        block_desc = grid.new_descriptor(2 * B, 2 * B, 8, 8)
        local_desc = grid.new_descriptor(2 * B, 2 * B, 2 * B, 2 * B)

        eps_M = np.zeros([2 * B])
        Hsc_MM = local_desc.zeros(dtype=self.dtype)
        Ssc_MM = local_desc.zeros(dtype=self.dtype)
        Csc_MM = local_desc.zeros(dtype=self.dtype)

        Hsc_mm = block_desc.zeros(dtype=self.dtype)
        Ssc_mm = block_desc.zeros(dtype=self.dtype)
        Csc_mm = block_desc.zeros(dtype=self.dtype)

        def integrate(a_G):
            if wfs.collinear:
                return np.real(wfs.integrate(a_G, a_G, global_integral=False))
            return sum(np.real(wfs.integrate(b_G, b_G, global_integral=False))
                       for b_G in a_G)

        self.subspace_diagonalize(ham, wfs, kpt)

        psit = kpt.psit
        psit2 = psit.new(buf=wfs.work_array)
        P = kpt.projections
        P2 = P.new()
        P3 = P.new()
        M = wfs.work_matrix_nn
        dS = wfs.setups.dS
        comm = wfs.gd.comm

        if bd.comm.size > 1:
            M0 = M.new(dist=(bd.comm, 1, 1))
        else:
            M0 = M

        if comm.rank == 0:
            e_N = bd.collect(kpt.eps_n)
            if e_N is not None:
                eps_N[:B] = e_N

        Ht = partial(wfs.apply_pseudo_hamiltonian, kpt, ham)

        if self.keep_htpsit:
            R = psit.new(buf=self.Htpsit_nG)
        else:
            R = psit.apply(Ht)

        self.calculate_residuals(kpt, wfs, ham, psit, P, kpt.eps_n, R, P2)

        precond = self.preconditioner

        for nit in range(self.niter):
            if nit == self.niter - 1:
                error = np.dot(weights, [integrate(R_G) for R_G in R.array])

            for psit_G, R_G, psit2_G in zip(psit.array, R.array, psit2.array):
                ekin = precond.calculate_kinetic_energy(psit_G, kpt)
                precond(R_G, kpt, ekin, out=psit2_G)

            # Calculate projections
            psit2.matrix_elements(wfs.pt, out=P2)

            def copy(M, C_nn):
                comm.sum(M.array, 0)
                if comm.rank == 0:
                    M.complex_conjugate()
                    M.redist(M0)
                    if bd.comm.rank == 0:
                        C_nn[:] = M0.array

            with self.timer('calc. matrices'):
                # <psi2 | H | psi2>
                psit2.matrix_elements(operator=Ht, result=R, out=M,
                                      symmetric=True, cc=True)
                ham.dH(P2, out=P3)
                mmm(1.0, P2, 'N', P3, 'C', 1.0, M)  # , symmetric=True)
                copy(M, H_NN[B:, B:])

                # <psi2 | H | psi>
                R.matrix_elements(psit, out=M, cc=True)
                mmm(1.0, P3, 'N', P, 'C', 1.0, M)
                copy(M, H_NN[B:, :B])

                # <psi2 | S | psi2>
                psit2.matrix_elements(out=M, symmetric=True, cc=True)
                dS.apply(P2, out=P3)
                mmm(1.0, P2, 'N', P3, 'C', 1.0, M)
                copy(M, S_NN[B:, B:])

                # <psi2 | S | psi> 
                psit2.matrix_elements(psit, out=M, cc=True)
                mmm(1.0, P3, 'N', P, 'C', 1.0, M)
                copy(M, S_NN[B:, :B])


            if slcomm.rank == 0:
                H_NN[:B, :B] = np.diag(eps_N[:B])
                S_NN[:B, :B] = np.eye(B)

                assert Hsc_MM.shape == H_NN.shape
                assert Ssc_MM.shape == H_NN.shape
                assert Csc_MM.shape == H_NN.shape
                Hsc_MM[:, :] = H_NN.copy()
                Ssc_MM[:, :] = S_NN.copy()
                # eps_M[:] = eps_N.copy()

            from gpaw.blacs import Redistributor
            redistributor = Redistributor(slcomm, local_desc, block_desc)
            redistributor.redistribute(Hsc_MM, Hsc_mm)
            redistributor.redistribute(Ssc_MM, Ssc_mm)
            redistributor.redistribute(Csc_MM, Csc_mm)

            scalapacked_davidson = False
            if scalapacked_davidson:
                block_desc.general_diagonalize_dc(Hsc_mm.copy(), Ssc_mm.copy(), Csc_mm, eps_M)
            # block_desc.inverse_cholesky(Csc_mm)

            redistributor2 = Redistributor(slcomm, block_desc, local_desc)
            redistributor2.redistribute(Csc_mm, Csc_MM, uplo='G')

            # if slcomm.rank == 0:
            #     import time as t
            #     np.savez(f'scalH{self.i}.npz', Hsc_MM)
            #     np.savez(f'scalCtrans{self.i}.npz', Csc_MM.T)

            with self.timer('diagonalize'):
                if slcomm.rank == 0 and comm.rank == 0 and bd.comm.rank == 0:
                    # H_NN[:B, :B] = np.diag(eps_N[:B])
                    # S_NN[:B, :B] = np.eye(B)
                    # np.savez(f'preteighHmat{self.i}.npz', H_NN)
                    # np.savez(f'preteighSmat{self.i}.npz', S_NN)
                    if debug:
                        H_NN[np.triu_indices(2 * B, 1)] = 42.0
                        S_NN[np.triu_indices(2 * B, 1)] = 42.0
                    from scipy.linalg import eigh, orth
                    
                    eps_N, scipy_eigvecs = eigh(H_NN, S_NN,
                                          lower=True,
                                          check_finite=debug)
                    np.set_printoptions(linewidth=200, precision=4, threshold=np.inf, suppress=True, floatmode='fixed')
                    if scalapacked_davidson
                        H_NN[:, :] = Csc_MM.T.copy()
                        eps_N[:] = eps_M.copy()
                    else:
                        H_NN[:, :] = scipy_eigvecs.copy()
                        # eps_N[:] = eps_M.copy()
                    # print("Scalapack consistency", np.linalg.norm(H_NN[:, :B], axis=0), np.linalg.norm(Csc_MM[:, :B], axis=0), sep='\n')
                    # print("Scalapack consistency", eps_N, eps_M, sep='\n')
                    # np.savez(f'posteighHmat{self.i}.npz', H_NN)
                
                    

            if comm.rank == 0:
                bd.distribute(eps_N[:B], kpt.eps_n)
            comm.broadcast(kpt.eps_n, 0)

            with self.timer('rotate_psi'):
                if comm.rank == 0:
                    if bd.comm.rank == 0:
                        M0.array[:] = H_NN[:B, :B].T
                    M0.redist(M)
                comm.broadcast(M.array, 0)
                mmm(1.0, M, 'N', psit, 'N', 0.0, R)
                mmm(1.0, M, 'N', P, 'N', 0.0, P3)
                if comm.rank == 0:
                    if bd.comm.rank == 0:
                        M0.array[:] = H_NN[B:, :B].T
                    M0.redist(M)
                comm.broadcast(M.array, 0)
                mmm(1.0, M, 'N', psit2, 'N', 1.0, R)
                mmm(1.0, M, 'N', P2, 'N', 1.0, P3)
                psit[:] = R
                P, P3 = P3, P
                kpt.projections = P

            if nit < self.niter - 1:
                psit.apply(Ht, out=R)
                self.calculate_residuals(kpt, wfs, ham, psit,
                                         P, kpt.eps_n, R, P2)

        error = wfs.gd.comm.sum(error)
        self.i += 1
        return error
