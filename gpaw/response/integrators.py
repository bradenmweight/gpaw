from __future__ import print_function, division

import sys

import numpy as np
from scipy.spatial import Delaunay

from ase.utils import devnull
from ase.utils.timing import timer, Timer

from _gpaw import tetrahedron_weight

import gpaw.mpi as mpi
from gpaw.occupations import FermiDirac
from gpaw.utilities.blas import gemm, rk, czher, mmm
from gpaw.utilities.progressbar import ProgressBar
from gpaw.utilities.memory import maxrss
from functools import partial


class Integrator():

    def __init__(self, cell_cv, comm=mpi.world,
                 txt=sys.stdout, timer=None,  nblocks=1):
        """Baseclass for Brillouin zone integration and band summation.

        Simple class to calculate integrals over Brilloun zones
        and summation of bands.

        comm: mpi.communicator
        nblocks: block parallelization
        """

        self.comm = comm
        self.nblocks = nblocks
        self.vol = abs(np.linalg.det(cell_cv))
        if nblocks == 1:
            self.blockcomm = self.comm.new_communicator([comm.rank])
            self.kncomm = comm
        else:
            assert comm.size % nblocks == 0, comm.size
            rank1 = comm.rank // nblocks * nblocks
            rank2 = rank1 + nblocks
            self.blockcomm = self.comm.new_communicator(range(rank1, rank2))
            ranks = range(comm.rank % nblocks, comm.size, nblocks)
            self.kncomm = self.comm.new_communicator(ranks)

        if comm.rank != 0:
            txt = devnull
        elif isinstance(txt, str):
            txt = open(txt, 'w')
        self.fd = txt

        self.timer = timer or Timer()

    def distribute_domain(self, domain_dl):
        """Distribute integration domain. """
        domainsize = [len(domain_l) for domain_l in domain_dl]
        nterms = np.prod(domainsize)
        size = self.kncomm.size
        rank = self.kncomm.rank

        n = (nterms + size - 1) // size
        i1 = rank * n
        i2 = min(i1 + n, nterms)
        mydomain = []
        for i in range(i1, i2):
            unravelled_d = np.unravel_index(i, domainsize)
            arguments = []
            for domain_l, index in zip(domain_dl, unravelled_d):
                arguments.append(domain_l[index])
            mydomain.append(tuple(arguments))

        print('Distributing domain %s' % (domainsize, ),
              'over %d process%s' %
              (self.kncomm.size, ['es', ''][self.kncomm.size == 1]),
              file=self.fd)
        print('Number of blocks:', self.blockcomm.size, file=self.fd)

        return mydomain

    def integrate(self, *args, **kwargs):
        raise NotImplementedError


class PointIntegrator(Integrator):

    def __init__(self, *args, **kwargs):

        """Integrate brillouin zone using a broadening technique.

        The broadening technique consists of smearing out the
        delta functions appearing in many integrals by some factor
        eta. In this code we use Lorentzians."""
        Integrator.__init__(self, *args, **kwargs)

    def integrate(self, kind=None, *args, **kwargs):
        if kind is None:
            return self.pointwise_integration(*args, **kwargs)
        elif kind == 'spectral function':
            return self.response_function_integration(hilbert=True,
                                                      *args, **kwargs)
        elif kind == 'hermitian response function':
            return self.response_function_integration(hermitian=True,
                                                      *args, **kwargs)
        elif kind == 'response function':
            return self.response_function_integration(*args, **kwargs)
        else:
            raise NotImplementedError

    def response_function_integration(self, domain=None, integrand=None,
                                      x=None, kwargs=None, out_wxx=None,
                                      timeordered=False, hermitian=False,
                                      intraband=False, hilbert=True):
        """Integrate a response function over bands and kpoints.

        func: method
        omega_w: ndarray
        out: np.ndarray
        timeordered: Bool
        """
        if out_wxx is None:
            raise NotImplementedError

        mydomain_t = self.distribute_domain(domain)
        nbz = len(domain[0])
        get_matrix_element, get_eigenvalues = integrand
        # The kwargs contain any constant
        # arguments provided by the user
        if kwargs is not None:
            get_matrix_element = partial(get_matrix_element,
                                         **kwargs[0])
            get_eigenvalues = partial(get_eigenvalues,
                                      **kwargs[1])

        # Sum kpoints
        # Calculate integrations weight
        pb = ProgressBar(self.fd)
        for _, arguments in pb.enumerate(mydomain_t):
            n_MG = get_matrix_element(*arguments)
            deps_M = get_eigenvalues(*arguments)

            if intraband:
                self.update_intraband(n_MG, out_wxx)
            elif hermitian:
                self.update_hermitian(n_MG, deps_M, x, out_wxx)
            elif hilbert:
                self.update_hilbert(n_MG, deps_M, x, out_wxx)
            else:
                self.update(n_MG, deps_M, x, out_wxx,
                            timeordered=timeordered)
        # Sum over
        for out_xx in out_wxx:
            self.kncomm.sum(out_xx)

        if (hermitian or hilbert) and self.blockcomm.size == 1:
            # Fill in upper/lower triangle also:
            nx = out_wxx.shape[1]
            il = np.tril_indices(nx, -1)
            iu = il[::-1]
            if hilbert:
                for out_xx in out_wxx:
                    out_xx[il] = out_xx[iu].conj()
            else:
                for out_xx in out_wxx:
                    out_xx[iu] = out_xx[il].conj()

        out_wxx *= (2 * np.pi)**3 / self.vol / nbz

    @timer('CHI_0 update')
    def update(self, n_mG, deps_m, wd, chi0_wGG, timeordered=False):
        """Update chi."""

        omega_w = wd.omega_w
        if timeordered:
            deps1_m = deps_m + 1j * self.eta * np.sign(deps_m)
            deps2_m = deps1_m
        else:
            deps1_m = deps_m + 1j * self.eta
            deps2_m = deps_m - 1j * self.eta

        for omega, chi0_GG in zip(omega_w, chi0_wGG):
            x_m = (1 / (omega + deps1_m) - 1 / (omega - deps2_m))
            if self.blockcomm.size > 1:
                nx_mG = n_mG[:, self.Ga:self.Gb] * x_m[:, np.newaxis]
            else:
                nx_mG = n_mG * x_m[:, np.newaxis]
            gemm(1.0, n_mG.conj(), np.ascontiguousarray(nx_mG.T),
                 1.0, chi0_GG)

    @timer('CHI_0 hermetian update')
    def update_hermitian(self, n_mG, deps_m, wd, chi0_wGG):
        """If eta=0 use hermitian update."""
        omega_w = wd.omega_w
        for w, omega in enumerate(omega_w):
            if self.blockcomm.size == 1:
                x_m = (-2 * deps_m / (omega.imag**2 + deps_m**2))**0.5
                nx_mG = n_mG.conj() * x_m[:, np.newaxis]
                rk(-1.0, nx_mG, 1.0, chi0_wGG[w], 'n')
            else:
                x_m = 2 * deps_m / (omega.imag**2 + deps_m**2)
                mynx_mG = n_mG[:, self.Ga:self.Gb] * x_m[:, np.newaxis]
                mmm(1.0, mynx_mG, 'c', n_mG, 'n', 1.0, chi0_wGG[w])

    @timer('CHI_0 spectral function update')
    def update_hilbert(self, n_mG, deps_m, wd, chi0_wGG):
        """Update spectral function.

        Updates spectral function A_wGG and saves it to chi0_wGG for
        later hilbert-transform."""

        self.timer.start('prep')
        o_m = abs(deps_m)
        w_m = wd.get_closest_index(o_m)
        o1_m = wd.omega_w[w_m]
        o2_m = wd.omega_w[w_m + 1]
        p_m = 1 / (o2_m - o1_m)**2  # XXX There was abs around df_m here
        p1_m = p_m * (o2_m - o_m)
        p2_m = p_m * (o_m - o1_m)
        self.timer.stop('prep')

        if self.blockcomm.size > 1:
            for p1, p2, n_G, w in zip(p1_m, p2_m, n_mG, w_m):
                myn_G = n_G[self.Ga:self.Gb].reshape((-1, 1))
                gemm(p1, n_G.reshape((-1, 1)), myn_G, 1.0, chi0_wGG[w], 'c')
                gemm(p2, n_G.reshape((-1, 1)), myn_G, 1.0, chi0_wGG[w + 1],
                     'c')
            return

        for p1, p2, n_G, w in zip(p1_m, p2_m, n_mG, w_m):
            czher(p1, n_G.conj(), chi0_wGG[w])
            czher(p2, n_G.conj(), chi0_wGG[w + 1])

    @timer('CHI_0 optical limit update')
    def update_optical_limit(self, n0_mv, deps_m, df_m, n_mG,
                             chi0_wxvG, chi0_wvv):
        """Optical limit update of chi."""

        if self.hilbert:  # Do something special when hilbert transforming
            self.update_optical_limit_hilbert(n0_mv, deps_m, df_m, n_mG,
                                              chi0_wxvG, chi0_wvv)
            return

        if self.timeordered:
            # avoid getting a zero from np.sign():
            deps1_m = deps_m + 1j * self.eta * np.sign(deps_m + 1e-20)
            deps2_m = deps1_m
        else:
            deps1_m = deps_m + 1j * self.eta
            deps2_m = deps_m - 1j * self.eta

        for w, omega in enumerate(omega_w):
            x_m = df_m * (1 / (omega + deps1_m) -
                          1 / (omega - deps2_m))

            chi0_wvv[w] += np.dot(x_m * n0_mv.T, n0_mv.conj())
            chi0_wxvG[w, 0, :, 1:] += np.dot(x_m * n0_mv.T, n_mG[:, 1:].conj())
            chi0_wxvG[w, 1, :, 1:] += np.dot(x_m * n0_mv.T.conj(), n_mG[:, 1:])

    @timer('CHI_0 optical limit hilbert-update')
    def update_optical_limit_hilbert(self, n0_mv, deps_m, df_m, n_mG,
                                     chi0_wxvG, chi0_wvv):
        """Optical limit update of chi-head and -wings."""

        beta = (2**0.5 - 1) * self.domega0 / self.omega2
        for deps, df, n0_v, n_G in zip(deps_m, df_m, n0_mv, n_mG):
            o = abs(deps)
            w = int(o / (self.domega0 + beta * o))
            if w + 2 > len(omega_w):
                break
            o1, o2 = omega_w[w:w + 2]
            assert o1 <= o <= o2, (o1, o, o2)

            p = abs(df) / (o2 - o1)**2  # XXX abs()?
            p1 = p * (o2 - o)
            p2 = p * (o - o1)
            x_vv = np.outer(n0_v, n0_v.conj())
            chi0_wvv[w] += p1 * x_vv
            chi0_wvv[w + 1] += p2 * x_vv
            x_vG = np.outer(n0_v, n_G[1:].conj())
            chi0_wxvG[w, 0, :, 1:] += p1 * x_vG
            chi0_wxvG[w + 1, 0, :, 1:] += p2 * x_vG
            chi0_wxvG[w, 1, :, 1:] += p1 * x_vG.conj()
            chi0_wxvG[w + 1, 1, :, 1:] += p2 * x_vG.conj()

    @timer('CHI_0 intraband update')
    def update_intraband(self, vel_mv, chi0_wvv):
        """Add intraband contributions"""

        for vel_v in vel_mv:
            x_vv = np.outer(vel_v, vel_v)
            chi0_wvv[0] += x_vv


class TetrahedronIntegrator(Integrator):
    """Integrate brillouin zone using tetrahedron integration.

    Tetrahedron integration uses linear interpolation of
    the eigenenergies and of the matrix elements
    between the vertices of the tetrahedron."""

    def __init__(self, *args, **kwargs):
        Integrator.__init__(self, *args, **kwargs)

    @timer('Tesselate')
    def tesselate(self, vertices):
        """Get tesselation descriptor."""
        td = Delaunay(vertices)

        td.volumes_s = None
        return td

    def get_simplex_volume(self, td, S):
        """Get volume of simplex S"""

        if td.volumes_s is not None:
            return td.volumes_s[S]

        td.volumes_s = np.zeros(td.nsimplex, float)
        for s in range(td.nsimplex):
            K_k = td.simplices[s]
            k_kc = td.points[K_k]
            volume = np.abs(np.linalg.det(k_kc[1:] - k_kc[0])) / 6.
            td.volumes_s[s] = volume

        return self.get_simplex_volume(td, S)

    def integrate(self, kind, *args, **kwargs):
        if kind == 'spectral function':
            return self.spectral_function_integration(*args, **kwargs)
        else:
            raise NotImplementedError

    @timer('Spectral function integration')
    def spectral_function_integration(self, domain=None, integrand=None,
                                      x=None, kwargs=None, out_wxx=None):
        """Integrate response function.

        Assume that the integral has the
        form of a response function. For the linear tetrahedron
        method it is possible calculate frequency dependent weights
        and do a point summation using these weights."""

        # Input domain
        td = self.tesselate(domain[0])
        args = domain[1:]
        get_matrix_element, get_eigenvalues = integrand

        # The kwargs contain any constant
        # arguments provided by the user
        if kwargs is not None:
            get_matrix_element = partial(get_matrix_element,
                                         **kwargs[0])
            get_eigenvalues = partial(get_eigenvalues,
                                      **kwargs[1])

        # Relevant quantities
        bzk_kc = td.points
        nk = len(bzk_kc)

        with self.timer('pts'):
            # Point to simplex
            pts_k = [[] for n in xrange(nk)]
            for s, K_k in enumerate(td.simplices):
                A_kv = np.append(td.points[K_k],
                                 np.ones(4)[:, np.newaxis], axis=1)

                D_kv = np.append((A_kv[:, :-1]**2).sum(1)[:, np.newaxis],
                                 A_kv, axis=1)
                a = np.linalg.det(D_kv[:, np.arange(5) != 0])

                if np.abs(a) < 1e-10:
                    continue

                for K in K_k:
                    pts_k[K].append(s)

            # Change to numpy arrays:
            for k in xrange(nk):
                pts_k[k] = np.array(pts_k[k], int)

        with self.timer('neighbours'):
            # Nearest neighbours
            neighbours_k = [None for n in xrange(nk)]

            for k in xrange(nk):
                neighbours_k[k] = np.unique(td.simplices[pts_k[k]])

        # Distribute everything
        myterms_t = self.distribute_domain(list(args) +
                                           [list(range(nk))])

        with self.timer('eigenvalues'):
            # Store eigenvalues
            deps_tMk = None  # t for term
            shape = [len(domain_l) for domain_l in args]
            nterms = int(np.prod(shape))

            for t in range(nterms):
                if nterms == 1:
                    arguments = ()
                else:
                    arguments = np.unravel_index(t, shape)
                for K in range(nk):
                    k_c = bzk_kc[K]
                    deps_M = get_eigenvalues(k_c, *arguments)
                    if deps_tMk is None:
                        deps_tMk = np.zeros([nterms] +
                                            list(deps_M.shape) +
                                            [nk], float)
                    deps_tMk[t, :, K] = deps_M

        with self.timer('Get indices'):
            # Store indices for frequencies
            indices_tMki = np.zeros(list(deps_tMk.shape) + [2], int)
            for t, deps_Mk in enumerate(deps_tMk):
                for K in xrange(nk):
                    teteps_Mk = deps_Mk[:, neighbours_k[K]]
                    try:
                        emin_M, emax_M = teteps_Mk.min(1), teteps_Mk.max(1)
                    except:
                        print(neighbours_k[K])
                        print(teteps_Mk)
                        raise
                    i0_M, i1_M = x.get_index_range(emin_M, emax_M)
                    indices_tMki[t, :, K, 0] = i0_M
                    indices_tMki[t, :, K, 1] = i1_M

        omega_w = x.get_data()

        # Calculate integrations weight
        pb = ProgressBar(self.fd)
        for _, arguments in pb.enumerate(myterms_t):
            K = arguments[-1]
            if len(shape) == 0:
                t = 0
            else:
                t = np.ravel_multi_index(arguments[:-1], shape)
            deps_Mk = deps_tMk[t]
            indices_Mi = indices_tMki[t, :, K]
            n_MG = get_matrix_element(bzk_kc[K],
                                      *arguments[:-1])
            for n_G, deps_k, I_i in zip(n_MG, deps_Mk,
                                        indices_Mi):
                i0, i1 = I_i[0], I_i[1]
                if i0 == i1:
                    continue

                W_w = self.get_kpoint_weight(K, deps_k,
                                             pts_k, omega_w[i0:i1],
                                             td)

                for iw, weight in enumerate(W_w):
                    czher(weight, n_G.conj(), out_wxx[i0 + iw])

        self.kncomm.sum(out_wxx)

        # Fill in upper/lower triangle also:
        nx = out_wxx.shape[1]
        il = np.tril_indices(nx, -1)
        iu = il[::-1]
        for out_xx in out_wxx:
            out_xx[il] = out_xx[iu].conj()

    @timer('Get kpoint weight')
    def get_kpoint_weight(self, K, deps_k, pts_k,
                          omega_w, td):
        # Find appropriate index range
        simplices_s = pts_k[K]
        W_w = np.zeros(len(omega_w), float)
        vol_s = self.get_simplex_volume(td, simplices_s)
        with self.timer('Tetrahedron weight'):
            tetrahedron_weight(deps_k, td.simplices, K,
                               simplices_s,
                               W_w, omega_w, vol_s)

        return W_w
