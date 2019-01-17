from ase.units import Hartree
import numpy as np
from gpaw.utilities.blas import mmm  # , dotc, dotu
from gpaw.directmin.tools import D_matrix, expm_ed
from gpaw.lcao.eigensolver import DirectLCAO
from scipy.linalg import expm  # , expm_frechet
from gpaw.utilities.tools import tri2full
from ase.utils import basestring
from gpaw.directmin.odd import odd_corrections
from gpaw.directmin import search_direction, line_search_algorithm
from gpaw.directmin.tools import loewdin
from gpaw.xc.__init__ import xc_string_to_dict


class DirectMinLCAO(DirectLCAO):

    def __init__(self, diagonalizer=None, error=np.inf,
                 search_direction_algorithm='LBFGS_P',
                 line_search_algorithm='SwcAwc',
                 initial_orbitals='KS',
                 initial_rotation='zero',
                 update_ref_orbs_counter=15,
                 update_precond_counter=1000,
                 use_prec=True, matrix_exp='pade_approx',
                 sparse=True, odd_parameters='Zero'):

        super(DirectMinLCAO, self).__init__(diagonalizer, error)

        self.sda = search_direction_algorithm
        self.lsa = line_search_algorithm
        self.initial_rotation = initial_rotation
        self.initial_orbitals = initial_orbitals
        self.get_en_and_grad_iters = 0
        self.update_ref_orbs_counter = update_ref_orbs_counter
        self.update_precond_counter = update_precond_counter
        self.use_prec = use_prec
        self.matrix_exp = matrix_exp
        self.sparse = sparse
        self.iters = 0
        self.name = 'direct_min'

        self.a_mat_u = None  # skew-hermitian matrix to be exponented
        self.g_mat_u = None  # gradient matrix
        self.c_nm_ref = None  # reference orbitals to be rotated

        self.odd_parameters = odd_parameters

        if isinstance(self.sda, basestring):
            self.sda = xc_string_to_dict(self.sda)
        if isinstance(self.lsa, basestring):
            self.lsa = xc_string_to_dict(self.lsa)
            self.lsa['method'] = self.sda['name']

    def __repr__(self):

        sds = {'SD': 'Steepest Descent',
               'FRcg': 'Fletcher-Reeves conj. grad. method',
               'HZcg': 'Hager-Zhang conj. grad. method',
               'QuickMin': 'Molecular-dynamics based algorithm',
               'LBFGS': 'LBFGS algorithm',
               'LBFGS_P': 'LBFGS algorithm with preconditioning'}

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

        self.n_dim = {}  # dimensionality of the problem.
                         # this implementation rotates among all bands
        for kpt in wfs.kpt_u:
            u = kpt.s * self.n_kps + kpt.q
            self.n_dim[u] = wfs.bd.nbands

        self.a_mat_u = {}  # skew-hermitian matrix to be exponented
        self.g_mat_u = {}  # gradient matrix
        self.c_nm_ref = {}  # reference orbitals to be rotated

        self.evecs = {}   # eigendecomposition for a
        self.evals = {}
        self.ind_up = None

        if self.sparse:
            # we may want to use different shapes for different
            # kpts, for example metals or sic, but later..
            max = wfs.bd.nbands
            nax = self.nvalence = wfs.nvalence
            nax = nax // (3 - wfs.nspins) + nax % (3 - wfs.nspins)
            ind_up = np.triu_indices(max, 1)
            x = (max - nax) * (max - nax - 1) // 2
            self.ind_up = (ind_up[0][:-x], ind_up[1][:-x])
            del ind_up

        for kpt in wfs.kpt_u:
            u = self.n_kps * kpt.s + kpt.q
            if self.sparse:
                shape_of_arr = len(self.ind_up[0])
            else:
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
        self.phi = [None, None]  # energy at alpha and alpha old
        self.der_phi = [None, None]  # gradients at alpha and al. old
        self.precond = None

        self.iters = 1
        self.nvalence = wfs.nvalence
        self.kd_comm = wfs.kd.comm
        self.heiss = {}  # heissian for LBFGS-P
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

    def iterate(self, ham, wfs, dens, occ):

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
            # first iteration is diagonilisation using super class
            # super().iterate(ham, wfs)
            if not wfs.coefficients_read_from_file and \
                    self.c_nm_ref is None:
                super().iterate(ham, wfs)
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
            self.update_ks_energy(ham, wfs, dens, occ)
            for kpt in wfs.kpt_u:
                self.sort_wavefunctions(ham, wfs, kpt)
            self.initialize_2(wfs, dens, ham)

        wfs.timer.start('Preconditioning:')
        precond = self.update_preconditioning(wfs, self.use_prec)
        self.update_ref_orbitals(wfs)  #, ham)
        wfs.timer.stop('Preconditioning:')

        a = self.a_mat_u
        n_dim = self.n_dim
        alpha = self.alpha
        phi = self.phi
        c_ref = self.c_nm_ref
        der_phi = self.der_phi

        if self.iters == 1:
            phi[0], g = self.get_energy_and_gradients(a, n_dim, ham,
                                                      wfs, dens, occ,
                                                      c_ref)
        else:
            g = self.g_mat_u

        wfs.timer.start('Get Search Direction')
        p = self.get_search_direction(a, g, precond, wfs)
        wfs.timer.stop('Get Search Direction')
        der_phi_c = 0.0
        for k in g.keys():
            if self.sparse:
                der_phi_c += np.dot(g[k].conj(),
                                    p[k]).real
            else:
                if self.dtype is complex:
                    il1 = np.tril_indices(g[k].shape[0])
                else:
                    il1 = np.tril_indices(g[k].shape[0], -1)
                der_phi_c += np.dot(g[k][il1].conj(), p[k][il1]).real
                # der_phi_c += dotc(g[k][il1], p[k][il1]).real
        der_phi_c = wfs.kd.comm.sum(der_phi_c)
        der_phi[0] = der_phi_c

        phi_c, der_phi_c = phi[0], der_phi[0]

        alpha, phi[0], der_phi[0], g = \
            self.line_search.step_length_update(a, p, n_dim,
                                                ham, wfs, dens, occ,
                                                c_ref,
                                                phi_0=phi[0],
                                                der_phi_0=der_phi[0],
                                                phi_old=phi[1],
                                                der_phi_old=der_phi[1],
                                                alpha_max=3.0,
                                                alpha_old=alpha)

        if wfs.gd.comm.size > 1:
            wfs.timer.start('Broadcast gradients')
            alpha_phi_der_phi = np.array([alpha, phi[0], der_phi[0]])
            wfs.gd.comm.broadcast(alpha_phi_der_phi, 0)
            alpha = alpha_phi_der_phi[0]
            phi[0] = alpha_phi_der_phi[1]
            der_phi[0] = alpha_phi_der_phi[2]
            for kpt in wfs.kpt_u:
                k = self.n_kps * kpt.s + kpt.q
                wfs.gd.comm.broadcast(g[k], 0)
            wfs.timer.stop('Broadcast gradients')

        phi[1], der_phi[1] = phi_c, der_phi_c

        # calculate new matrices for optimal step length
        for k in a.keys():
            a[k] += alpha * p[k]
        self.g_mat_u = g
        self.iters += 1

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

        wfs.timer.start('Unitary rotation')
        for kpt in wfs.kpt_u:
            k = self.n_kps * kpt.s + kpt.q
            if n_dim[k] == 0:
                continue

            if self.gd.comm.rank == 0:
                if self.sparse:
                    a = np.zeros(shape=(n_dim[k], n_dim[k]),
                                 dtype=self.dtype)
                    a[self.ind_up] = a_mat_u[k]
                    a += -a.T.conj()
                else:
                    a = a_mat_u[k]

                if self.matrix_exp == 'pade_approx':
                    # this function takes a lot of memory
                    # for large matrices... what can we do?
                    wfs.timer.start('Pade Approximants')
                    u_nn = expm(a)
                    del a
                    wfs.timer.stop('Pade Approximants')
                elif self.matrix_exp == 'eigendecomposition':
                    # this method is based on diagonalisation
                    wfs.timer.start('Eigendecomposition')
                    u_nn, self.evecs[k], self.evals[k] =\
                        expm_ed(a, evalevec=True)
                    wfs.timer.stop('Eigendecomposition')
                else:
                    raise NotImplementedError('Check the keyword '
                                              'for matrix_exp. \n'
                                              'Must be '
                                              '\'pade_approx\' or '
                                              '\'eigendecomposition\'')
                kpt.C_nM[:n_dim[k]] = np.dot(u_nn.T,
                                             c_nm_ref[k][:n_dim[k]])
                del u_nn

            self.gd.comm.broadcast(kpt.C_nM, 0)
            wfs.atomic_correction.calculate_projections(wfs, kpt)
        wfs.timer.stop('Unitary rotation')

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
                                                       self.sparse,
                                                       self.ind_up)
            if hasattr(self.odd, 'e_sic_by_orbitals'):
                self.e_sic += self.odd.e_sic_by_orbitals[k].sum()

            self._error += error
        self._error = self.kd_comm.sum(self._error)
        self.e_sic = self.kd_comm.sum(self.e_sic)
        wfs.timer.stop('Calculate gradients')

        self.get_en_and_grad_iters += 1

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
        if self.matrix_exp == 'pade_approx':
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
        elif self.matrix_exp == 'eigendecomposition':
            timer.start('Use Eigendecomposition')
            grad = evec.T.conj() @ h_mm @ evec
            grad = grad * D_matrix(evals)
            grad = evec @ grad @ evec.T.conj()
            timer.stop('Use Eigendecomposition')
            for i in range(grad.shape[0]):
                grad[i][i] *= 0.5
        else:
            raise NotImplementedError('Check the keyword '
                                      'for matrix_exp. \n'
                                      'Must be '
                                      '\'pade_approx\' or '
                                      '\'eigendecomposition\'')
        if self.dtype == float:
            grad = grad.real
        if self.sparse:
            grad = grad[self.ind_up]
        timer.stop('Construct Gradient Matrix')

        return 2.0 * grad, error

    def get_search_direction(self, a_mat_u, g_mat_u, precond, wfs):

        if self.sparse:
            p_mat_u = self.search_direction.update_data(wfs, a_mat_u,
                                                        g_mat_u,
                                                        precond)
        else:
            g_vec = {}
            a_vec = {}

            for k in a_mat_u.keys():
                if self.dtype is complex:
                    il1 = np.tril_indices(a_mat_u[k].shape[0])
                else:
                    il1 = np.tril_indices(a_mat_u[k].shape[0], -1)

                a_vec[k] = a_mat_u[k][il1]
                g_vec[k] = g_mat_u[k][il1]

            p_vec = self.search_direction.update_data(wfs, a_vec,
                                                      g_vec, precond)
            del a_vec, g_vec

            p_mat_u = {}
            for k in p_vec.keys():
                p_mat_u[k] = np.zeros_like(a_mat_u[k])
                if self.dtype is complex:
                    il1 = np.tril_indices(a_mat_u[k].shape[0])
                else:
                    il1 = np.tril_indices(a_mat_u[k].shape[0], -1)
                p_mat_u[k][il1] = p_vec[k]
                # make it skew-hermitian
                ind_l = np.tril_indices(p_mat_u[k].shape[0], -1)
                p_mat_u[k][(ind_l[1], ind_l[0])] = -p_mat_u[k][ind_l].conj()
            del p_vec

        return p_mat_u

    def evaluate_phi_and_der_phi(self, a_mat_u, p_mat_u, n_dim, alpha,
                                 ham, wfs, dens, occ, c_ref,
                                 phi=None, g_mat_u=None):
        """
        phi = f(x_k + alpha_k*p_k)
        der_phi = \grad f(x_k + alpha_k*p_k) \cdot p_k
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
        if self.sparse:
            for k in p_mat_u.keys():
                der_phi += np.dot(g_mat_u[k].conj(),
                                  p_mat_u[k]).real
        else:
            for k in p_mat_u.keys():
                if self.dtype is complex:
                    il1 = np.tril_indices(p_mat_u[k].shape[0])
                else:
                    il1 = np.tril_indices(p_mat_u[k].shape[0], -1)

                der_phi += np.dot(g_mat_u[k][il1].conj(),
                                  p_mat_u[k][il1]).real
                # der_phi += dotc(g_mat_u[k][il1],
                #                 p_mat_u[k][il1]).real

        der_phi = wfs.kd.comm.sum(der_phi)

        return phi, der_phi, g_mat_u

    def update_ref_orbitals(self, wfs):  # , ham):
        counter = self.update_ref_orbs_counter
        if self.iters % counter == 0 and self.iters > 1:
            for kpt in wfs.kpt_u:
                u = kpt.s * self.n_kps + kpt.q
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
                        heiss = self.get_hessian(kpt)
                        if self.dtype is float:
                            self.precond[k] = np.zeros_like(heiss)
                            for i in range(heiss.shape[0]):
                                if abs(heiss[i]) < 1.0e-4:
                                    self.precond[k][i] = 1.0
                                else:
                                    self.precond[k][i] = \
                                        1.0 / (heiss[i].real)
                        else:
                            self.precond[k] = np.zeros_like(heiss)
                            for i in range(heiss.shape[0]):
                                if abs(heiss[i]) < 1.0e-4:
                                    self.precond[k][i] = 1.0 + 1.0j
                                else:
                                    self.precond[k][i] = 1.0 / \
                                                         heiss[i].real + \
                                                         1.0j / \
                                                         heiss[i].imag
                    return self.precond
                else:
                    return self.precond
            else:
                # it's a bit messy, here you store self.heis,
                # but in 'if' above self.precond
                precond = {}
                for kpt in wfs.kpt_u:
                    k = self.n_kps * kpt.s + kpt.q
                    if self.iters % counter == 0 or self.iters == 1:
                        self.heiss[k] = self.get_hessian(kpt)
                    heiss = self.heiss[k]
                    if self.dtype is float:
                        precond[k] = 1.0 / (
                                0.75 * heiss +
                                0.25 * self.search_direction.beta_0 ** (-1))
                    else:
                        precond[k] = \
                            1.0 / (0.75 * heiss.real +
                                   0.25 * self.search_direction.beta_0 ** (-1)) + \
                            1.0j / (0.75 * heiss.imag +
                                    0.25 * self.search_direction.beta_0 ** (-1))
                return precond
        else:
            return None

    def get_hessian(self, kpt):

        f_n = kpt.f_n
        eps_n = kpt.eps_n
        if self.sparse:
            il1 = list(self.ind_up)
        else:
            if self.dtype is complex:
                il1 = np.tril_indices(eps_n.shape[0])
            else:
                il1 = np.tril_indices(eps_n.shape[0], -1)
            il1 = list(il1)

        heiss = np.zeros(len(il1[0]), dtype=self.dtype)
        x = 0
        for l, m in zip(*il1):
            df = f_n[l] - f_n[m]
            heiss[x] = -2.0 * (eps_n[l] - eps_n[m]) * df
            if self.dtype is complex:
                heiss[x] += 1.0j * heiss[x]
                if abs(heiss[x]) < 1.0e-10:
                    heiss[x] = 0.0 + 0.0j
            else:
                if abs(heiss[x]) < 1.0e-10:
                    heiss[x] = 0.0
            x += 1

        return heiss

    def calculate_residual(self, kpt, H_MM, S_MM, wfs):
        return np.inf

    def get_canonical_representation(self, ham, wfs, dens, occ):

        # choose canonical orbitals which diagonalise
        # lagrange matrix. need to do subspace rotation with equally
        # occupied states?
        # I tried to call function below, but due to instability
        # of eigensolvers
        # for some systems, it can 'mess' the solution.
        # this usually happens in metals,
        # the so-called charge-sloshing problem..
        wfs.timer.start('Get Canonical Representation')
        if self.odd.name == 'PZ_SIC':
            for kpt in wfs.kpt_u:
                h_mm = self.calculate_hamiltonian_matrix(ham, wfs, kpt)
                tri2full(h_mm)
                self.odd.get_lagrange_matrices(h_mm, kpt.C_nM,
                                               kpt.f_n, kpt, wfs,
                                               update_eigenvalues=True)
        elif self.odd.name == 'Zero':
            super().iterate(ham, wfs)
            occ.calculate(wfs)
            self.initialize_2(wfs, dens, ham)
            self.update_ks_energy(ham, wfs, dens, occ)

        for kpt in wfs.kpt_u:
            u = kpt.s * self.n_kps + kpt.q
            self.c_nm_ref[u] = kpt.C_nM.copy()
            self.a_mat_u[u] = np.zeros_like(self.a_mat_u[u])

        wfs.timer.stop('Get Canonical Representation')

    def reset(self):
        super().reset()
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
                'search_direction_algorithm': self.sda,
                'line_search_algorithm': self.lsa,
                'initial_orbitals': 'KS',
                'initial_rotation':'zero',
                'update_ref_orbs_counter': self.update_ref_orbs_counter,
                'update_precond_counter': self.update_precond_counter,
                'use_prec': self.use_prec,
                'matrix_exp': self.matrix_exp,
                'sparse': self.sparse,
                'odd_parameters': self.odd_parameters}
