from ase.units import Hartree
import numpy as np
from gpaw.utilities.blas import mmm  # , dotc, dotu
from gpaw.directmin.lcao.tools import D_matrix, expm_ed, expm_ed_unit_inv
from gpaw.lcao.eigensolver import DirectLCAO
from scipy.linalg import expm  # , expm_frechet
from gpaw.utilities.tools import tri2full
from gpaw.directmin.lcao import search_direction, line_search_algorithm
from gpaw.xc import xc_string_to_dict
from ase.utils import basestring
from gpaw.directmin.odd.lcao import odd_corrections
from gpaw.directmin.lcao.tools import loewdin
from gpaw.pipekmezey.pipek_mezey_wannier import PipekMezey as pm
from gpaw.pipekmezey.wannier_basic import WannierLocalization as wl


class DirectMinLCAO(DirectLCAO):

    def __init__(self, diagonalizer=None, error=np.inf,
                 searchdir_algo='LBFGS_P',
                 linesearch_algo='SwcAwc',
                 initial_orbitals=None,
                 initial_rotation='zero',  # not used right now
                 update_ref_orbs_counter=20,
                 update_precond_counter=1000,
                 use_prec=True, matrix_exp='pade_approx',
                 representation='sparse',
                 odd_parameters='Zero',
                 init_from_ks_eigsolver=False):

        super(DirectMinLCAO, self).__init__(diagonalizer, error)

        self.sda = searchdir_algo
        self.lsa = linesearch_algo
        self.initial_rotation = initial_rotation
        self.initial_orbitals = initial_orbitals
        self.eg_count = 0
        self.update_ref_orbs_counter = update_ref_orbs_counter
        self.update_precond_counter = update_precond_counter
        self.use_prec = use_prec
        self.matrix_exp = matrix_exp
        self.representation = representation
        self.iters = 0
        self.name = 'direct_min'

        self.a_mat_u = None  # skew-hermitian matrix to be exponented
        self.g_mat_u = None  # gradient matrix
        self.c_nm_ref = None  # reference orbitals to be rotated

        self.odd_parameters = odd_parameters
        self.init_from_ks_eigsolver = init_from_ks_eigsolver

        if isinstance(self.odd_parameters, basestring):
            self.odd_parameters = xc_string_to_dict(self.odd_parameters)
        if isinstance(self.sda, basestring):
            self.sda = xc_string_to_dict(self.sda)
        if isinstance(self.lsa, basestring):
            self.lsa = xc_string_to_dict(self.lsa)
            self.lsa['method'] = self.sda['name']

        if isinstance(self.representation, basestring):
            assert self.representation in ['sparse', 'u_invar', 'full'], \
                'Value Error'
            self.representation = \
                xc_string_to_dict(self.representation)

        if self.odd_parameters['name'] == 'PZ_SIC':
            if self.initial_orbitals is None:
                self.initial_orbitals = 'FB'

        if self.sda['name'] == 'LBFGS_P' and not self.use_prec:
            raise ValueError('Use LBFGS_P with use_prec=True')

        if matrix_exp == 'egdecomp2':
            assert self.representation['name'] == 'u_invar', \
                'Use u_invar representation with egdecomp2'

    def __repr__(self):

        sds = {'SD': 'Steepest Descent',
               'FRcg': 'Fletcher-Reeves conj. grad. method',
               'HZcg': 'Hager-Zhang conj. grad. method',
               'QuickMin': 'Molecular-dynamics based algorithm',
               'LBFGS': 'LBFGS algorithm',
               'LBFGS_P': 'LBFGS algorithm with preconditioning',
               'LBFGS_P2': 'LBFGS algorithm with preconditioning',
               'LSR1P': 'Limited-memory SR1P algorithm'}

        lss = {'UnitStep': 'step size equals one',
               'Parabola': 'Parabolic line search',
               'SwcAwc': 'Inexact line search based '
                         'on cubic interpolation,\n'
                         '                    strong'
                         ' and approximate Wolfe conditions'}

        repr_string = 'Direct minimisation using exponential ' \
                      'transformation.\n'
        repr_string += '       ' \
                       'Search ' \
                       'direction: {}\n'.format(sds[self.sda['name']])
        repr_string += '       ' \
                       'Line ' \
                       'search: {}\n'.format(lss[self.lsa['name']])
        repr_string += '       ' \
                       'Preconditioning: {}\n'.format(self.use_prec)
        repr_string += '       ' \
                       'WARNING: do not use it for metals as ' \
                       'occupation numbers are\n' \
                       '                ' \
                       'not found variationally\n'

        return repr_string

    def initialize_2(self, wfs, dens, ham):

        self.dtype = wfs.dtype
        self.n_kps = wfs.kd.nks // wfs.kd.nspins

        # dimensionality of the problem.
        # this implementation rotates among all bands
        self.n_dim = {}
        for kpt in wfs.kpt_u:
            u = kpt.s * self.n_kps + kpt.q
            self.n_dim[u] = wfs.bd.nbands

        # values: matrices, keys: kpt number
        self.a_mat_u = {}  # skew-hermitian matrix to be exponented
        self.g_mat_u = {}  # gradient matrix
        self.c_nm_ref = {}  # reference orbitals to be rotated

        self.evecs = {}   # eigendecomposition for a
        self.evals = {}
        self.ind_up = {}

        if self.representation['name'] in ['sparse', 'u_invar']:
            # Matrices are sparse and Skew-Hermitian.
            # They have this structure:
            #  A_BigMatrix =
            #
            # (  A_1          A_2 )
            # ( -A_2.T.conj() 0   )
            #
            # where 0 is a zero-matrix of size of (M-N) * (M-N)
            #
            # A_1 i skew-hermitian matrix of N * N,
            # N-number of occupied states
            # A_2 is matrix of size of (M-N) * N,
            # M - number of basis functions
            #
            # if the energy functional is unitary invariant
            # then A_1 = 0
            # (see Hutter J et. al, J. Chem. Phys. 101, 3862 (1994))
            #
            # We will keep A_1 as we would like to work with metals,
            # SIC, and molecules with different occupation numbers.
            # this corresponds to 'sparse' representation
            #
            # Thus, for the 'sparse' we need to store upper
            # triangular part of A_1, and matrix A_2, so in total
            # (M-N) * N + N * (N - 1)/2 = N * (M - (N + 1)/2) elements
            #
            # we will store these elements as a vector and
            # also will store indices of the A_BigMatrix
            # which correspond to these elements.
            #
            # 'u_invar' corresponds to the case when we want to
            # store only A_2, that is this representaion is sparser

            M = wfs.bd.nbands  # M - one dimension of the A_BigMatrix
            if self.representation['name'] == 'sparse':
                # let's take all upper triangular indices
                # of A_BigMatrix and then delete indices from ind_up
                # which correspond to 0 matrix in in A_BigMatrix.
                ind_up = np.triu_indices(M, 1)
                for kpt in wfs.kpt_u:
                    n_occ = get_n_occ(kpt)
                    u = self.n_kps * kpt.s + kpt.q
                    zero_ind = ((M - n_occ) * (M - n_occ - 1)) // 2
                    self.ind_up[u] = (ind_up[0][:-zero_ind].copy(),
                                      ind_up[1][:-zero_ind].copy())
                del ind_up
            else:
                # take indices of A_2 only
                for kpt in wfs.kpt_u:
                    n_occ = get_n_occ(kpt)
                    u = self.n_kps * kpt.s + kpt.q
                    i1, i2 = [], []
                    for i in range(n_occ):
                        for j in range(n_occ, M):
                            i1.append(i)
                            i2.append(j)
                    self.ind_up[u] = (np.asarray(i1), np.asarray(i2))

        for kpt in wfs.kpt_u:
            u = self.n_kps * kpt.s + kpt.q
            if self.representation['name'] in ['sparse', 'u_invar']:
                shape_of_arr = len(self.ind_up[u][0])
            else:
                self.ind_up[u] = None
                shape_of_arr = (self.n_dim[u], self.n_dim[u])

            self.a_mat_u[u] = np.zeros(shape=shape_of_arr,
                                       dtype=self.dtype)
            self.g_mat_u[u] = np.zeros(shape=shape_of_arr,
                                       dtype=self.dtype)
            # use initial KS orbitals, but can be others
            self.c_nm_ref[u] = np.copy(kpt.C_nM[:self.n_dim[u]])
            self.evecs[u] = None
            self.evals[u] = None

        self.alpha = 1.0  # step length
        self.phi_2i = [None, None]  # energy at last two iterations
        self.der_phi_2i = [None, None]  # energy gradient w.r.t. alpha
        self.precond = None

        self.iters = 1
        self.nvalence = wfs.nvalence
        self.nbands = wfs.bd.nbands
        self.kd_comm = wfs.kd.comm
        self.hess = {}  # hessian for LBFGS-P
        self.precond = {}  # precondiner for other methods

        # choose search direction and line search algorithm
        if isinstance(self.sda, (basestring, dict)):
            self.search_direction = search_direction(self.sda, wfs)
        else:
            raise Exception('Check Search Direction Parameters')

        if isinstance(self.lsa, (basestring, dict)):
            self.line_search = \
                line_search_algorithm(self.lsa,
                                      self.evaluate_phi_and_der_phi)
        else:
            raise Exception('Check Search Direction Parameters')

        # odd corrections
        if isinstance(self.odd_parameters, (basestring, dict)):
            self.odd = odd_corrections(self.odd_parameters, wfs,
                                       dens, ham)
        elif self.odd is None:
            pass
        else:
            raise Exception('Check ODD Parameters')
        self.e_sic = 0.0

    def iterate(self, ham, wfs, dens, occ, log):

        assert dens.mixer.driver.name == 'dummy', \
            'Please, use: mixer=DummyMixer()'
        assert wfs.bd.nbands == wfs.basis_functions.Mmax, \
            'Please, use: nbands=\'nao\''
        assert wfs.bd.comm.size == 1, \
            'Band parallelization is not supported'
        assert occ.width < 1.0e-5, \
            'Zero Kelvin only.'

        wfs.timer.start('Direct Minimisation step')

        if self.iters == 0:
            # need to initialize c_nm, eps, f_n and so on.
            self.init_wave_functions(wfs, ham, occ, log)
            self.update_ks_energy(ham, wfs, dens, occ)
            # not sure sort wfs is good here,
            # you need to probably run loop over sort_wfs
            # with update of energy. So, don't use it for now.
            # for kpt in wfs.kpt_u:
            #     self.sort_wavefunctions(ham, wfs, kpt)
            self.initialize_2(wfs, dens, ham)

        wfs.timer.start('Preconditioning:')
        precond = self.update_preconditioning(wfs, self.use_prec)
        wfs.timer.stop('Preconditioning:')
        self.update_ref_orbitals(wfs, ham)

        if str(self.search_direction) == 'LBFGS_P2':
            for kpt in wfs.kpt_u:
                u = kpt.s * self.n_kps + kpt.q
                self.c_nm_ref[u] = kpt.C_nM.copy()

        a_mat_u = self.a_mat_u
        n_dim = self.n_dim
        alpha = self.alpha
        phi_2i = self.phi_2i
        der_phi_2i = self.der_phi_2i
        c_ref = self.c_nm_ref

        if self.iters == 1:
            phi_2i[0], g_mat_u = \
                self.get_energy_and_gradients(a_mat_u, n_dim, ham, wfs,
                                              dens, occ, c_ref)
        else:
            g_mat_u = self.g_mat_u

        wfs.timer.start('Get Search Direction')
        p_mat_u = self.get_search_direction(a_mat_u, g_mat_u, precond,
                                            wfs)
        wfs.timer.stop('Get Search Direction')

        # recalculate derivative with new search direction
        der_phi_2i[0] = 0.0
        for k in g_mat_u.keys():
            if self.representation['name'] in ['sparse', 'u_invar']:
                der_phi_2i[0] += np.dot(g_mat_u[k].conj(),
                                        p_mat_u[k]).real
            else:
                il1 = get_indices(g_mat_u[k].shape[0], self.dtype)
                der_phi_2i[0] += np.dot(g_mat_u[k][il1].conj(),
                                        p_mat_u[k][il1]).real
                # der_phi_c += dotc(g[k][il1], p[k][il1]).real
        der_phi_2i[0] = wfs.kd.comm.sum(der_phi_2i[0])

        if str(self.search_direction) == 'LBFGS_P2':
            for kpt in wfs.kpt_u:
                u = kpt.s * self.n_kps + kpt.q
                a_mat_u[u] = np.zeros_like(a_mat_u[u])

        alpha, phi_alpha, der_phi_alpha, g_mat_u = \
            self.line_search.step_length_update(a_mat_u, p_mat_u,
                                                n_dim, ham, wfs, dens,
                                                occ, c_ref,
                                                phi_0=phi_2i[0],
                                                der_phi_0=der_phi_2i[0],
                                                phi_old=phi_2i[1],
                                                der_phi_old=der_phi_2i[1],
                                                alpha_max=5.0,
                                                alpha_old=alpha,
                                                kpdescr=wfs.kd)

        if wfs.gd.comm.size > 1:
            wfs.timer.start('Broadcast gradients')
            alpha_phi_der_phi = np.array([alpha, phi_2i[0],
                                          der_phi_2i[0]])
            wfs.gd.comm.broadcast(alpha_phi_der_phi, 0)
            alpha = alpha_phi_der_phi[0]
            phi_2i[0] = alpha_phi_der_phi[1]
            der_phi_2i[0] = alpha_phi_der_phi[2]
            for kpt in wfs.kpt_u:
                k = self.n_kps * kpt.s + kpt.q
                wfs.gd.comm.broadcast(g_mat_u[k], 0)
            wfs.timer.stop('Broadcast gradients')

        # calculate new matrices for optimal step length
        for k in a_mat_u.keys():
            if str(self.search_direction) == 'LBFGS_P2':
                a_mat_u[k] = alpha * p_mat_u[k]
            else:
                a_mat_u[k] += alpha * p_mat_u[k]
        self.alpha = alpha
        self.g_mat_u = g_mat_u
        self.iters += 1

        # and 'shift' phi, der_phi for the next iteration
        phi_2i[1], der_phi_2i[1] = phi_2i[0], der_phi_2i[0]
        phi_2i[0], der_phi_2i[0] = phi_alpha, der_phi_alpha,

        wfs.timer.stop('Direct Minimisation step')

    def get_energy_and_gradients(self, a_mat_u, n_dim, ham, wfs, dens,
                                 occ, c_nm_ref):

        """
        Energy E = E[C exp(A)]. Gradients G_ij[C, A] = dE/dA_ij

        :param a_mat_u: A
        :param c_nm_ref: C
        :param n_dim:
        :return:
        """

        self.rotate_wavefunctions(wfs, a_mat_u, n_dim, c_nm_ref)

        e_total = self.update_ks_energy(ham, wfs, dens, occ)

        wfs.timer.start('Calculate gradients')
        g_mat_u = {}
        self._error = 0.0
        self.e_sic = 0.0  # this is odd energy
        for kpt in wfs.kpt_u:
            k = self.n_kps * kpt.s + kpt.q
            if n_dim[k] == 0:
                g_mat_u[k] = np.zeros_like(a_mat_u[k])
                continue
            h_mm = self.calculate_hamiltonian_matrix(ham, wfs, kpt)
            # make matrix hermitian
            tri2full(h_mm)
            g_mat_u[k], error = self.odd.get_gradients(h_mm, kpt.C_nM,
                                                       kpt.f_n,
                                                       self.evecs[k],
                                                       self.evals[k],
                                                       kpt, wfs,
                                                       wfs.timer,
                                                       self.matrix_exp,
                                                       self.representation['name'],
                                                       self.ind_up[k])
            if hasattr(self.odd, 'e_sic_by_orbitals'):
                self.e_sic += self.odd.e_sic_by_orbitals[k].sum()

            self._error += error
        self._error = self.kd_comm.sum(self._error)
        self.e_sic = self.kd_comm.sum(self.e_sic)
        wfs.timer.stop('Calculate gradients')

        self.eg_count += 1

        return e_total + self.e_sic, g_mat_u

    def update_ks_energy(self, ham, wfs, dens, occ):

        wfs.timer.start('Update Kohn-Sham energy')
        dens.update(wfs)
        ham.update(dens, wfs, False)
        wfs.timer.stop('Update Kohn-Sham energy')

        return ham.get_energy(occ, False)

    def get_gradients(self, h_mm, c_nm, f_n, evec, evals, kpt, timer):

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
        if self.matrix_exp in ['pade_approx', 'egdecomp2']:
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
        elif self.matrix_exp == 'egdecomp':
            timer.start('Use Eigendecomposition')
            grad = np.dot(evec.T.conj(), np.dot(h_mm, evec))
            grad = grad * D_matrix(evals)
            grad = np.dot(evec, np.dot(grad, evec.T.conj()))
            timer.stop('Use Eigendecomposition')
            for i in range(grad.shape[0]):
                grad[i][i] *= 0.5
        else:
            raise ValueError('Check the keyword '
                             'for matrix_exp. \n'
                             'Must be '
                             '\'pade_approx\' or '
                             '\'egdecomp\'')

        if self.dtype == float:
            grad = grad.real
        if self.representation['name'] in ['sparse', 'u_invar']:
            u = self.n_kps * kpt.s + kpt.q
            grad = grad[self.ind_up[u]]
        timer.stop('Construct Gradient Matrix')

        return 2.0 * grad, error

    def get_search_direction(self, a_mat_u, g_mat_u, precond, wfs):

        if self.representation['name'] in ['sparse', 'u_invar']:
            p_mat_u = self.search_direction.update_data(wfs, a_mat_u,
                                                        g_mat_u,
                                                        precond)
        else:
            g_vec = {}
            a_vec = {}

            for k in a_mat_u.keys():
                il1 = get_indices(a_mat_u[k].shape[0], self.dtype)
                a_vec[k] = a_mat_u[k][il1]
                g_vec[k] = g_mat_u[k][il1]

            p_vec = self.search_direction.update_data(wfs, a_vec,
                                                      g_vec, precond)
            del a_vec, g_vec

            p_mat_u = {}
            for k in p_vec.keys():
                p_mat_u[k] = np.zeros_like(a_mat_u[k])
                il1 = get_indices(p_mat_u[k].shape[0], self.dtype)
                p_mat_u[k][il1] = p_vec[k]
                # make it skew-hermitian
                il1 = np.tril_indices(p_mat_u[k].shape[0], -1)
                p_mat_u[k][(il1[1], il1[0])] = -p_mat_u[k][il1].conj()

            del p_vec

        return p_mat_u

    def evaluate_phi_and_der_phi(self, a_mat_u, p_mat_u, n_dim, alpha,
                                 ham, wfs, dens, occ, c_ref,
                                 phi=None, g_mat_u=None):
        """
        phi = f(x_k + alpha_k*p_k)
        der_phi = \\grad f(x_k + alpha_k*p_k) \\cdot p_k
        :return:  phi, der_phi # floats
        """
        if phi is None or g_mat_u is None:
            x_mat_u = {k: a_mat_u[k] + alpha * p_mat_u[k]
                       for k in a_mat_u.keys()}
            phi, g_mat_u = \
                self.get_energy_and_gradients(x_mat_u, n_dim,
                                              ham, wfs, dens, occ,
                                              c_ref
                                              )
            del x_mat_u
        else:
            pass

        der_phi = 0.0
        if self.representation['name'] in ['sparse', 'u_invar']:
            for k in p_mat_u.keys():
                der_phi += np.dot(g_mat_u[k].conj(),
                                  p_mat_u[k]).real
        else:
            for k in p_mat_u.keys():

                il1 = get_indices(p_mat_u[k].shape[0], self.dtype)

                der_phi += np.dot(g_mat_u[k][il1].conj(),
                                  p_mat_u[k][il1]).real
                # der_phi += dotc(g_mat_u[k][il1],
                #                 p_mat_u[k][il1]).real

        der_phi = wfs.kd.comm.sum(der_phi)

        return phi, der_phi, g_mat_u

    def update_ref_orbitals(self, wfs, ham):
        if str(self.search_direction) == 'LBFGS_P2':
            return 0

        counter = self.update_ref_orbs_counter
        if self.iters % counter == 0 and self.iters > 1:
            # TODO: you need to recompute the gradients now
            for kpt in wfs.kpt_u:
                u = kpt.s * self.n_kps + kpt.q
                # self.sort_wavefunctions(ham, wfs, kpt)
                self.c_nm_ref[u] = kpt.C_nM.copy()
                self.a_mat_u[u] = np.zeros_like(self.a_mat_u[u])
                # self.sort_wavefunctions(ham, wfs, kpt)

            # choose search direction and line search algorithm
            if isinstance(self.sda, (basestring, dict)):
                self.search_direction = search_direction(self.sda, wfs)
            else:
                raise Exception('Check Search Direction Parameters')

            if isinstance(self.lsa, (basestring, dict)):
                self.line_search = \
                    line_search_algorithm(self.lsa,
                                          self.evaluate_phi_and_der_phi)
            else:
                raise Exception('Check Search Direction Parameters')

    def update_preconditioning(self, wfs, use_prec):
        counter = self.update_precond_counter
        if use_prec:
            if self.sda['name'] != 'LBFGS_P':
                if self.iters % counter == 0 or self.iters == 1:
                    for kpt in wfs.kpt_u:
                        k = self.n_kps * kpt.s + kpt.q
                        hess = self.get_hessian(kpt)
                        if self.dtype is float:
                            self.precond[k] = np.zeros_like(hess)
                            for i in range(hess.shape[0]):
                                if abs(hess[i]) < 1.0e-4:
                                    self.precond[k][i] = 1.0
                                else:
                                    self.precond[k][i] = \
                                        1.0 / (hess[i].real)
                        else:
                            self.precond[k] = np.zeros_like(hess)
                            for i in range(hess.shape[0]):
                                if abs(hess[i]) < 1.0e-4:
                                    self.precond[k][i] = 1.0 + 1.0j
                                else:
                                    self.precond[k][i] = \
                                        1.0 / hess[i].real + \
                                        1.0j / hess[i].imag
                    return self.precond
                else:
                    return self.precond
            else:
                # it's a bit messy, here you store self.heis,
                # but in 'if' above self.precond
                precond = {}
                for kpt in wfs.kpt_u:
                    k = self.n_kps * kpt.s + kpt.q
                    w = kpt.weight / (3.0 - wfs.nspins)
                    if self.iters % counter == 0 or self.iters == 1:
                        self.hess[k] = self.get_hessian(kpt)
                    hess = self.hess[k]
                    if self.dtype is float:
                        precond[k] = \
                            1.0 / (0.75 * hess +
                                   w * 0.25 * self.search_direction.beta_0 ** (-1))
                    else:
                        precond[k] = \
                            1.0 / (0.75 * hess.real +
                                   w * 0.25 * self.search_direction.beta_0 ** (-1)) + \
                            1.0j / (0.75 * hess.imag +
                                   w * 0.25 * self.search_direction.beta_0 ** (-1))
                return precond
        else:
            return None

    def get_hessian(self, kpt):

        f_n = kpt.f_n
        eps_n = kpt.eps_n
        if self.representation['name'] in ['sparse', 'u_invar']:
            u = self.n_kps * kpt.s + kpt.q
            il1 = list(self.ind_up[u])
        else:
            il1 = get_indices(eps_n.shape[0], self.dtype)
            il1 = list(il1)

        hess = np.zeros(len(il1[0]), dtype=self.dtype)
        x = 0
        for l, m in zip(*il1):
            df = f_n[l] - f_n[m]
            hess[x] = -2.0 * (eps_n[l] - eps_n[m]) * df
            if self.dtype is complex:
                hess[x] += 1.0j * hess[x]
                if abs(hess[x]) < 1.0e-10:
                    hess[x] = 0.0 + 0.0j
            else:
                if abs(hess[x]) < 1.0e-10:
                    hess[x] = 0.0
            x += 1

        return hess

    def calculate_residual(self, kpt, H_MM, S_MM, wfs):
        return np.inf

    def get_canonical_representation(self, ham, wfs):

        # choose canonical orbitals which diagonalise
        # lagrange matrix. it's probably necessary
        # to do subspace rotation with equally
        # occupied states.
        # In this case, the total energy remains the same,
        # as it's unitary invariant within equally occupied subspaces.
        wfs.timer.start('Get canonical representation')

        for kpt in wfs.kpt_u:
            h_mm = self.calculate_hamiltonian_matrix(ham, wfs, kpt)
            tri2full(h_mm)
            if self.odd.name == 'Zero':
                n_init = 0
                while True:
                    n_fin = find_equally_occupied_subspace(kpt.f_n, n_init)
                    kpt.C_nM[n_init:n_init + n_fin, :], kpt.eps_n[n_init:n_init + n_fin] = \
                        rotate_subspace(h_mm, kpt.C_nM[n_init:n_init + n_fin, :])
                    n_init += n_fin
                    if n_init == len(kpt.f_n):
                        break
                    elif n_init > len(kpt.f_n):
                        raise SystemExit('Bug is here!')
                    # this function rotates occpuied subspace
                    # n_occ = 0
                    # nbands = len(kpt.f_n)
                    # while n_occ < nbands and kpt.f_n[n_occ] > 1e-10:
                    #     n_occ += 1
                    # n_unocc = len(kpt.f_n) - n_occ
                    #
                    # for x in [0, n_occ]:
                    #     y = (x // n_occ) * (n_unocc - n_occ) + n_occ
                    #     kpt.C_nM[x:x + y, :], kpt.eps_n[x:x + y] = \
                    #         rotate_subspace(h_mm, kpt.C_nM[x:x + y, :])
            elif self.odd.name == 'PZ_SIC':
                self.odd.get_lagrange_matrices(h_mm, kpt.C_nM,
                                               kpt.f_n, kpt, wfs,
                                               update_eigenvalues=True)
            u = kpt.s * self.n_kps + kpt.q
            self.c_nm_ref[u] = kpt.C_nM.copy()
            self.a_mat_u[u] = np.zeros_like(self.a_mat_u[u])

        wfs.timer.stop('Get canonical representation')

    def reset(self):
        super(DirectMinLCAO, self).reset()
        self._error = np.inf
        self.iters = 0

    def sort_wavefunctions(self, ham, wfs, kpt):
        """
        this function sorts wavefunctions according
        it's orbitals energies, not eigenvalues.
        :return:
        """
        wfs.timer.start('Sort WFS')
        h_mm = self.calculate_hamiltonian_matrix(ham, wfs, kpt)
        tri2full(h_mm)
        hc_mn = np.zeros(shape=(kpt.C_nM.shape[1], kpt.C_nM.shape[0]),
                         dtype=kpt.C_nM.dtype)
        mmm(1.0, h_mm.conj(), 'N', kpt.C_nM, 'T', 0.0, hc_mn)
        mmm(1.0, kpt.C_nM.conj(), 'N', hc_mn, 'N', 0.0, h_mm)
        orbital_energies = h_mm.diagonal().real.copy()
        # label each orbital energy
        # add some noise to get rid off degeneracy
        orbital_energies += \
            np.random.rand(len(orbital_energies)) * 1.0e-8
        oe_labeled = {}
        for i, lamb in enumerate(orbital_energies):
            oe_labeled[str(round(lamb, 12))] = i
        # now sort orb energies
        oe_sorted = np.sort(orbital_energies)
        # run over sorted orbital energies and get their label
        ind = []
        for x in oe_sorted:
            i = oe_labeled[str(round(x, 12))]
            ind.append(i)
        # check if it is necessary to sort
        x = np.max(abs(np.array(ind) - np.arange(len(ind))))
        if x > 0:
            # now sort wfs according to orbital energies
            kpt.C_nM[np.arange(len(ind)),:] = kpt.C_nM[ind,:]
            kpt.eps_n[:] = np.sort(h_mm.diagonal().real.copy())
        wfs.timer.stop('Sort WFS')

        return

    def todict(self):
        return {'name': 'direct_min_lcao',
                'searchdir_algo': self.sda,
                'linesearch_algo': self.lsa,
                'initial_orbitals': self.initial_orbitals,
                'initial_rotation': 'zero',
                'update_ref_orbs_counter': self.update_ref_orbs_counter,
                'update_precond_counter': self.update_precond_counter,
                'use_prec': self.use_prec,
                'matrix_exp': self.matrix_exp,
                'representation': self.representation,
                'odd_parameters': self.odd_parameters}

    def get_numerical_gradients(self, n_dim, ham, wfs, dens, occ,
                                c_nm_ref, eps=1.0e-7):

        assert not self.representation['name'] in ['sparse', 'u_invar']
        a_m = {}
        g_n = {}
        if self.matrix_exp == 'pade_approx':
            c_nm_ref = {}
        for kpt in wfs.kpt_u:
            u = self.n_kps * kpt.s + kpt.q
            a = np.random.random_sample(self.a_mat_u[u].shape) + \
                1.0j * np.random.random_sample(self.a_mat_u[u].shape)
            a = a - a.T.conj()
            u_nn = expm(a)
            g_n[u] = np.zeros_like(self.a_mat_u[u])

            if self.matrix_exp == 'pade_approx':
                a_m[u] = np.zeros_like(self.a_mat_u[u])
                c_nm_ref[u] = np.dot(u_nn.T, kpt.C_nM[:u_nn.shape[0]])
            elif self.matrix_exp == 'egdecomp':
                a_m[u] = a

        g_a = self.get_energy_and_gradients(a_m, n_dim, ham, wfs,
                                            dens, occ, c_nm_ref)[1]

        h = [eps, -eps]
        coeif = [1.0, -1.0]

        if self.dtype == complex:
            range_z = 2
            complex_gr = [1.0, 1.0j]
        else:
            range_z = 1
            complex_gr = [1.0]

        for kpt in wfs.kpt_u:
            u = self.n_kps * kpt.s + kpt.q
            dim = a_m[u].shape[0]
            for z in range(range_z):
                for i in range(dim):
                    for j in range(dim):
                        print(u, z, i, j)
                        a = a_m[u][i][j]
                        g = 0.0
                        for l in range(2):
                            if z == 0:
                                if i != j:
                                    a_m[u][i][j] = a + h[l]
                                    a_m[u][j][i] = -np.conjugate(a + h[l])
                            else:
                                a_m[u][i][j] = a + 1.0j * h[l]
                                if i != j:
                                    a_m[u][j][i] = -np.conjugate(a + 1.0j * h[l])

                            E = self.get_energy_and_gradients(a_m, n_dim, ham, wfs, dens, occ, c_nm_ref)[0]

                            g += E * coeif[l]

                        g *= 1.0 / (2.0 * eps)

                        g_n[u][i][j] += g * complex_gr[z]
                        a_m[u][i][j] = a
                        if i != j:
                            a_m[u][j][i] = -np.conjugate(a)

        return g_a, g_n

    def rotate_wavefunctions(self, wfs, a_mat_u, n_dim, c_nm_ref):

        wfs.timer.start('Unitary rotation')
        for kpt in wfs.kpt_u:
            k = self.n_kps * kpt.s + kpt.q
            if n_dim[k] == 0:
                continue

            if self.gd.comm.rank == 0:
                if self.representation['name'] in ['sparse', 'u_invar']:
                    if self.matrix_exp == 'egdecomp2' and \
                            self.representation['name'] == 'u_invar':
                        n_occ = get_n_occ(kpt)
                        n_v = self.nbands - n_occ
                        a = a_mat_u[k].reshape(n_occ, n_v)
                    else:
                        a = np.zeros(shape=(n_dim[k], n_dim[k]),
                                     dtype=self.dtype)
                        a[self.ind_up[k]] = a_mat_u[k]
                        a += -a.T.conj()
                else:
                    a = a_mat_u[k]

                if self.matrix_exp == 'pade_approx':
                    # this function takes a lot of memory
                    # for large matrices... what can we do?
                    wfs.timer.start('Pade Approximants')
                    u_nn = expm(a)
                    wfs.timer.stop('Pade Approximants')
                elif self.matrix_exp == 'egdecomp':
                    # this method is based on diagonalisation
                    wfs.timer.start('Eigendecomposition')
                    u_nn, evecs, evals = \
                        expm_ed(a, evalevec=True)
                    wfs.timer.stop('Eigendecomposition')
                elif self.matrix_exp == 'egdecomp2':
                    assert self.representation['name'] == 'u_invar'
                    wfs.timer.start('Eigendecomposition')
                    u_nn = expm_ed_unit_inv(a, oo_vo_blockonly=False)
                    wfs.timer.stop('Eigendecomposition')

                else:
                    raise ValueError('Check the keyword '
                                     'for matrix_exp. \n'
                                     'Must be '
                                     '\'pade_approx\' or '
                                     '\'egdecomp\'')

                dimens1 = u_nn.shape[0]
                dimens2 = u_nn.shape[1]
                kpt.C_nM[:dimens2] = np.dot(u_nn.T,
                                             c_nm_ref[k][:dimens1])

                del u_nn
                del a

            wfs.timer.start('Broadcast coefficients')
            self.gd.comm.broadcast(kpt.C_nM, 0)
            wfs.timer.stop('Broadcast coefficients')

            if self.matrix_exp == 'egdecomp':
                wfs.timer.start('Broadcast evecs and evals')
                if self.gd.comm.rank != 0:
                    evecs = np.zeros(shape=(n_dim[k], n_dim[k]),
                                     dtype=complex)
                    evals = np.zeros(shape=n_dim[k],
                                     dtype=float)

                self.gd.comm.broadcast(evecs, 0)
                self.gd.comm.broadcast(evals, 0)
                self.evecs[k], self.evals[k] = evecs, evals
                wfs.timer.stop('Broadcast evecs and evals')

            wfs.atomic_correction.calculate_projections(wfs, kpt)

        wfs.timer.stop('Unitary rotation')

    def init_wave_functions(self, wfs, ham, occ, log):

        # if it is the first use of the scf then initialize
        # coefficient matrix using eigensolver
        # and then localise orbitals
        if (not wfs.coefficients_read_from_file and
                self.c_nm_ref is None) or self.init_from_ks_eigsolver:
            super(DirectMinLCAO, self).iterate(ham, wfs)
            occ.calculate(wfs)
            self.localize_wfs(wfs, log)

        # if one want to use coefficients saved in gpw file
        # or to use coefficients from the previous scf circle
        else:
            for kpt in wfs.kpt_u:
                u = kpt.s * wfs.kd.nks // wfs.kd.nspins + kpt.q
                if self.c_nm_ref is not None:
                    C = self.c_nm_ref[u]
                else:
                    C = kpt.C_nM
                kpt.C_nM[:] = loewdin(C, kpt.S_MM.conj())
            wfs.coefficients_read_from_file = False
            occ.calculate(wfs)

    def localize_wfs(self, wfs, log):

        log("Initial Localization: ...", flush=True)
        wfs.timer.start('Initial Localization')
        for kpt in wfs.kpt_u:
            if sum(kpt.f_n) < 1.0e-3:
                continue
            if self.initial_orbitals == 'KS' or \
                    self.initial_orbitals is None:
                continue
            elif self.initial_orbitals == 'PM':
                lf_obj = pm(wfs=wfs, spin=kpt.s,
                            dtype=wfs.dtype)
            elif self.initial_orbitals == 'FB':
                lf_obj = wl(wfs=wfs, spin=kpt.s)
            else:
                raise ValueError('Check initial orbitals.')
            lf_obj.localize(tolerance=1.0e-5)
            if self.initial_orbitals == 'PM':
                U = np.ascontiguousarray(
                    lf_obj.W_k[kpt.q].T)
            else:
                U = np.ascontiguousarray(
                    lf_obj.U_kww[kpt.q].T)
                if kpt.C_nM.dtype == float:
                    U = U.real
            wfs.gd.comm.broadcast(U, 0)
            dim = U.shape[0]
            kpt.C_nM[:dim] = np.dot(U, kpt.C_nM[:dim])
            del lf_obj
        wfs.timer.stop('Initial Localization')
        log("Done", flush=True)


def get_indices(dimens, dtype):

    if dtype == complex:
        il1 = np.tril_indices(dimens)
    else:
        il1 = np.tril_indices(dimens, -1)

    return il1


def get_n_occ(kpt):
    nbands = len(kpt.f_n)
    n_occ = 0
    while n_occ < nbands and kpt.f_n[n_occ] > 1e-10:
        n_occ += 1
    return n_occ


def find_equally_occupied_subspace(f_n, index=0):
    n_occ = 0
    f1 = f_n[index]
    for f2 in f_n[index:]:
        if abs(f1 - f2) < 1.0e-8:
            n_occ += 1
        else:
            return n_occ
    return n_occ


def rotate_subspace(h_mm, c_nm):
    l_nn = np.dot(np.dot(c_nm, h_mm), c_nm.conj().T).conj()
    # check if diagonal then don't rotate? it could save a bit of time
    eps, w = np.linalg.eigh(l_nn)
    return w.T.conj() @ c_nm, eps
