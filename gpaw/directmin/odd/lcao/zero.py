import numpy as np
from gpaw.utilities.blas import mmm
from ase.units import Hartree
from gpaw.directmin.lcao.tools import D_matrix

class ZeroCorrectionsLcao:
    """
    Don't apply any corrections

    """
    def __init__(self, wfs, dens, ham, **kwargs):
        self.name = 'Zero'
        self.n_kps = wfs.kd.nibzkpts
        self.dtype = wfs.dtype
        self.nvalence = wfs.nvalence

    def get_gradients(self, h_mm, c_nm, f_n,
                      evec, evals, kpt, wfs, timer, matrix_exp,
                      repr_name, ind_up,
                      occupied_only=False):
        # FIXME: what to do with occupied only here?
        #  also wfs is not used, but wfs is needed in other
        #  odd corrections see PzCorrectionsLcao
        #

        timer.start('Construct Gradient Matrix')
        hc_mn = np.zeros(shape=(c_nm.shape[1], c_nm.shape[0]),
                         dtype=self.dtype)
        mmm(1.0, h_mm.conj(), 'N', c_nm, 'T', 0.0, hc_mn)
        if c_nm.shape[0] != c_nm.shape[1]:
            h_mm = np.zeros(shape=(c_nm.shape[0], c_nm.shape[0]),
                            dtype=self.dtype)
        mmm(1.0, c_nm.conj(), 'N', hc_mn, 'N', 0.0, h_mm)
        timer.stop('Construct Gradient Matrix')

        # let's also calculate residual here.
        # it's extra calculation though, maybe it's better to use
        # norm of grad as convergence criteria..
        timer.start('Residual')
        n_occ = 0
        nbands = len(f_n)
        while n_occ < nbands and f_n[n_occ] > 1e-10:
            n_occ += 1
        # what if there are empty states between occupied?
        rhs = np.zeros(shape=(c_nm.shape[1], n_occ),
                       dtype=self.dtype)
        rhs2 = np.zeros(shape=(c_nm.shape[1], n_occ),
                        dtype=self.dtype)
        mmm(1.0, kpt.S_MM.conj(), 'N', c_nm[:n_occ], 'T', 0.0, rhs)
        mmm(1.0, rhs, 'N', h_mm[:n_occ, :n_occ], 'N', 0.0, rhs2)
        hc_mn = hc_mn[:, :n_occ] - rhs2[:, :n_occ]
        norm = []
        for i in range(n_occ):
            norm.append(np.dot(hc_mn[:,i].conj(),
                               hc_mn[:,i]).real * kpt.f_n[i])
            # needs to be contig. to use this:
            # x = np.ascontiguousarray(hc_mn[:,i])
            # norm.append(dotc(x, x).real * kpt.f_n[i])

        error = sum(norm) * Hartree ** 2 / self.nvalence
        del rhs, rhs2, hc_mn, norm
        timer.stop('Residual')

        # continue with gradients
        timer.start('Construct Gradient Matrix')
        h_mm = f_n * h_mm - f_n[:, np.newaxis] * h_mm
        if matrix_exp in ['pade_approx', 'egdecomp2']:
            # timer.start('Frechet derivative')
            # frechet derivative, unfortunately it calculates unitary
            # matrix which we already calculated before. Could it be used?
            # it also requires a lot of memory so don't use it now
            # u, grad = expm_frechet(a_mat, h_mm,
            #                        compute_expm=True,
            #                        check_finite=False)
            # grad = grad @ u.T.conj()
            # timer.stop('Frechet derivative')
            grad = np.ascontiguousarray(h_mm)
        elif matrix_exp == 'egdecomp':
            timer.start('Use Eigendecomposition')
            grad = np.dot(evec.T.conj(), np.dot(h_mm, evec))
            grad = grad * D_matrix(evals)
            grad = np.dot(evec, np.dot(grad, evec.T.conj()))
            for i in range(grad.shape[0]):
                grad[i][i] *= 0.5
            timer.stop('Use Eigendecomposition')
        else:
            raise NotImplementedError('Check the keyword '
                                      'for matrix_exp. \n'
                                      'Must be '
                                      '\'pade_approx\' or '
                                      '\'egdecomp\'')
        if self.dtype == float:
            grad = grad.real
        if repr_name in ['sparse', 'u_invar']:
            grad = grad[ind_up]
        timer.stop('Construct Gradient Matrix')

        return 2.0 * grad, error
