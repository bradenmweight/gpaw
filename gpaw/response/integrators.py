from abc import ABC, abstractmethod
import numpy as np
from gpaw.response import timer
from scipy.spatial import Delaunay
from scipy.linalg.blas import zher

import _gpaw
from gpaw.utilities.blas import rk, mmm
from gpaw.utilities.progressbar import ProgressBar
from gpaw.response.pw_parallelization import Blocks1D, block_partition


class Integrand(ABC):
    @abstractmethod
    def matrix_element(self, k_v, s):
        ...

    @abstractmethod
    def eigenvalues(self, k_v, s):
        ...


def czher(alpha: float, x, A) -> None:
    """Hermetian rank-1 update of upper half of A.

    A += alpha * np.outer(x.conj(), x)

    """
    AT = A.T
    out = zher(alpha, x, 1, 1, 0, len(x), AT, 1)
    assert out is AT


class Integrator:
    def __init__(self, cell_cv, context, nblocks=1, eshift=0.0):
        """Baseclass for Brillouin zone integration and band summation.

        Simple class to calculate integrals over Brilloun zones
        and summation of bands.

        context: ResponseContext
        nblocks: block parallelization
        """

        self.context = context
        self.eshift = eshift
        self.nblocks = nblocks
        self.vol = abs(np.linalg.det(cell_cv))

        self.blockcomm, self.kncomm = block_partition(self.context.comm,
                                                      nblocks)

    def distribute_domain(self, domain_dl):
        """Distribute integration domain. """
        domainsize = [len(domain_l) for domain_l in domain_dl]
        nterms = np.prod(domainsize)
        size = self.kncomm.size
        rank = self.kncomm.rank

        n = (nterms + size - 1) // size
        i1 = min(rank * n, nterms)
        i2 = min(i1 + n, nterms)
        assert i1 <= i2
        mydomain = []
        for i in range(i1, i2):
            unravelled_d = np.unravel_index(i, domainsize)
            arguments = []
            for domain_l, index in zip(domain_dl, unravelled_d):
                arguments.append(domain_l[index])
            mydomain.append(tuple(arguments))

        self.context.print('Distributing domain %s' % (domainsize,),
                           'over %d process%s' %
                           (self.kncomm.size,
                            ['es', ''][self.kncomm.size == 1]),
                           flush=False)
        self.context.print('Number of blocks:', self.blockcomm.size)

        return mydomain

    def integrate(self, kind, **kwargs):
        raise NotImplementedError

    def _blocks1d(self, nG):
        return Blocks1D(self.blockcomm, nG)


class PointIntegrator(Integrator):
    """Integrate brillouin zone using a broadening technique.

    The broadening technique consists of smearing out the
    delta functions appearing in many integrals by some factor
    eta. In this code we use Lorentzians."""

    def integrate(self, *, kind, task, wd, **kwargs):
        self.context.print('Integral kind:', kind)
        assert kind == task.kind or task.kind == 'intraband', (kind, task.kind)
        if kind == 'hermitian response function':
            assert task is not None
            dct = dict(hermitian=True,
                       hilbert=False,
                       wings=False)
        elif kind == 'hermitian response function wings':
            assert task is not None
            dct = dict(hermitian=True,
                       hilbert=False,
                       wings=True)
        elif kind == 'spectral function':
            assert task is not None
            dct = dict(hilbert=True)
        elif kind == 'spectral function wings':
            assert task is not None
            dct = dict(hilbert=True, wings=True)
        elif kind == 'response function':
            assert task is not None
            dct = dict(hilbert=False)
        elif kind == 'response function wings':
            assert task is not None
            dct = dict(hilbert=False,
                       wings=True)
        else:
            raise ValueError(kind)

        dct.update(kwargs)
        return self.response_function_integration(task=task, wd=wd, **dct)

    def response_function_integration(self, *, domain, integrand,
                                      wd, out_wxx,
                                      hermitian=False,
                                      hilbert=False,
                                      wings=False,
                                      eta=None,
                                      task):
        """Integrate a response function over bands and kpoints.

        func: method
        omega_w: ndarray
        out: np.ndarray
        """
        mydomain_t = self.distribute_domain(domain)
        nbz = len(domain[0])

        prefactor = (2 * np.pi)**3 / self.vol / nbz
        out_wxx /= prefactor

        # Sum kpoints
        # Calculate integrations weight
        pb = ProgressBar(self.context.fd)
        for _, arguments in pb.enumerate(mydomain_t):
            n_MG = integrand.matrix_element(*arguments)
            if n_MG is None:
                continue
            deps_M = integrand.eigenvalues(*arguments)

            task.run(wd, n_MG, deps_M, out_wxx)

        # Sum over
        # Can this really be valid, if the original input out_wxx is nonzero?
        # This smells and should be investigated XXX
        # There could also be similar errors elsewhere... XXX
        self.kncomm.sum(out_wxx)

        symmetrizable = (hermitian or hilbert) and not wings
        assert symmetrizable == task.symmetrizable_unless_blocked or task.kind == 'intraband', task

        if self.blockcomm.size == 1 and symmetrizable:
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

        out_wxx *= prefactor


from abc import ABC, abstractmethod

class IntegralKind(ABC):
    # Unique string for each kind of integral:
    kind = '(unset)'

    # Some integrals kinds like to calculate upper or lower half of the output
    # when nblocks==1.  In that case, this boolean signifies to the
    # integrator that the output array should be symmetrized.
    symmetrizable_unless_blocked = False

    @abstractmethod
    def run(self, wd, n_mG, deps_m, out_wxx):
        """Add contribution from one point to out_wxx."""


class GenericUpdate(IntegralKind):
    kind = 'response function'
    symmetrizable_unless_blocked = False

    def __init__(self, eta, integrator):
        self.eta = eta
        self.integrator = integrator

    # @timer('CHI_0 update')
    def run(self, wd, n_mG, deps_m, chi0_wGG):
        """Update chi."""

        deps_m += self.integrator.eshift * np.sign(deps_m)
        deps1_m = deps_m + 1j * self.eta
        deps2_m = deps_m - 1j * self.eta

        blocks1d = self.integrator._blocks1d(chi0_wGG.shape[2])

        for omega, chi0_GG in zip(wd.omega_w, chi0_wGG):
            x_m = (1 / (omega + deps1_m) - 1 / (omega - deps2_m))
            if blocks1d.blockcomm.size > 1:
                nx_mG = n_mG[:, blocks1d.myslice] * x_m[:, np.newaxis]
            else:
                nx_mG = n_mG * x_m[:, np.newaxis]

            mmm(1.0, np.ascontiguousarray(nx_mG.T), 'N', n_mG.conj(), 'N',
                1.0, chi0_GG)


class Hermitian(IntegralKind):
    kind = 'hermitian response function'
    symmetrizable_unless_blocked = True

    def __init__(self, integrator):
        self.integrator = integrator

    # @timer('CHI_0 hermetian update')
    def run(self, wd, n_mG, deps_m, chi0_wGG):
        """If eta=0 use hermitian update."""
        deps_m += self.integrator.eshift * np.sign(deps_m)

        blocks1d = self.integrator._blocks1d(chi0_wGG.shape[2])

        for w, omega in enumerate(wd.omega_w):
            if blocks1d.blockcomm.size == 1:
                x_m = np.abs(2 * deps_m / (omega.imag**2 + deps_m**2))**0.5
                nx_mG = n_mG.conj() * x_m[:, np.newaxis]
                rk(-1.0, nx_mG, 1.0, chi0_wGG[w], 'n')
            else:
                x_m = np.abs(2 * deps_m / (omega.imag**2 + deps_m**2))
                mynx_mG = n_mG[:, blocks1d.myslice] * x_m[:, np.newaxis]
                mmm(-1.0, mynx_mG, 'T', n_mG.conj(), 'N', 1.0, chi0_wGG[w])


class Hilbert(IntegralKind):
    kind = 'spectral function'
    symmetrizable_unless_blocked = True

    def __init__(self, integrator):
        self.integrator = integrator
        self.eshift = integrator.eshift

    # @timer('CHI_0 spectral function update (new)')
    def run(self, wd, n_mG, deps_m, chi0_wGG):
        """Update spectral function.

        Updates spectral function A_wGG and saves it to chi0_wGG for
        later hilbert-transform."""

        deps_m += self.eshift * np.sign(deps_m)
        o_m = abs(deps_m)
        w_m = wd.get_floor_index(o_m)

        # blocks1d = self._blocks1d(chi0_wGG.shape[2])
        blocks1d = self.integrator._blocks1d(chi0_wGG.shape[2])

        # Sort frequencies
        argsw_m = np.argsort(w_m)
        sortedo_m = o_m[argsw_m]
        sortedw_m = w_m[argsw_m]
        sortedn_mG = n_mG[argsw_m]

        index = 0
        while 1:
            if index == len(sortedw_m):
                break

            w = sortedw_m[index]
            startindex = index
            while 1:
                index += 1
                if index == len(sortedw_m):
                    break
                if w != sortedw_m[index]:
                    break

            endindex = index

            # Here, we have same frequency range w, for set of
            # electron-hole excitations from startindex to endindex.
            o1 = wd.omega_w[w]
            o2 = wd.omega_w[w + 1]
            p = np.abs(1 / (o2 - o1)**2)
            p1_m = np.array(p * (o2 - sortedo_m[startindex:endindex]))
            p2_m = np.array(p * (sortedo_m[startindex:endindex] - o1))

            if blocks1d.blockcomm.size > 1 and w + 1 < wd.wmax:
                x_mG = sortedn_mG[startindex:endindex, blocks1d.myslice]
                mmm(1.0,
                    np.concatenate((p1_m[:, None] * x_mG,
                                    p2_m[:, None] * x_mG),
                                   axis=1).T.copy(),
                    'N',
                    sortedn_mG[startindex:endindex].T.copy(),
                    'C',
                    1.0,
                    chi0_wGG[w:w + 2].reshape((2 * blocks1d.nlocal,
                                               blocks1d.N)))

            if blocks1d.blockcomm.size <= 1 and w + 1 < wd.wmax:
                x_mG = sortedn_mG[startindex:endindex]
                l_Gm = (p1_m[:, None] * x_mG).T.copy()
                r_Gm = x_mG.T.copy()
                mmm(1.0, r_Gm, 'N', l_Gm, 'C', 1.0, chi0_wGG[w])
                l_Gm = (p2_m[:, None] * x_mG).T.copy()
                mmm(1.0, r_Gm, 'N', l_Gm, 'C', 1.0, chi0_wGG[w + 1])


class Intraband(IntegralKind):
    kind = 'intraband'
    symmetrizable_unless_blocked = False

    # @timer('CHI_0 intraband update')
    def run(self, wd, vel_mv, deps_M, chi0_wvv):
        """Add intraband contributions"""
        # Intraband is a little bit special, we use neither wd nor deps_M

        for vel_v in vel_mv:
            x_vv = np.outer(vel_v, vel_v)
            chi0_wvv[0] += x_vv


class OpticalLimit(IntegralKind):
    kind = 'response function wings'
    symmetrizable_unless_blocked = False

    def __init__(self, eta):
        self.eta = eta

    # @timer('CHI_0 optical limit update')
    def run(self, wd, n_mG, deps_m, chi0_wxvG):
        """Optical limit update of chi."""
        deps1_m = deps_m + 1j * self.eta
        deps2_m = deps_m - 1j * self.eta

        for w, omega in enumerate(wd.omega_w):
            x_m = (1 / (omega + deps1_m) - 1 / (omega - deps2_m))
            chi0_wxvG[w, 0] += np.dot(x_m * n_mG[:, :3].T, n_mG.conj())
            chi0_wxvG[w, 1] += np.dot(x_m * n_mG[:, :3].T.conj(), n_mG)


class HermitianOpticalLimit(IntegralKind):
    kind = 'hermitian response function wings'
    symmetrizable_unless_blocked = False

    # @timer('CHI_0 hermitian optical limit update')
    def run(self, wd, n_mG, deps_m, chi0_wxvG):
        """Optical limit update of hermitian chi."""
        for w, omega in enumerate(wd.omega_w):
            x_m = - np.abs(2 * deps_m / (omega.imag**2 + deps_m**2))
            chi0_wxvG[w, 0] += np.dot(x_m * n_mG[:, :3].T, n_mG.conj())
            chi0_wxvG[w, 1] += np.dot(x_m * n_mG[:, :3].T.conj(), n_mG)


class HilbertOpticalLimit(IntegralKind):
    kind = 'spectral function wings'
    symmetrizable_unless_blocked = False

    # @timer('CHI_0 optical limit hilbert-update')
    def run(self, wd, n_mG, deps_m, chi0_wxvG):
        """Optical limit update of chi-head and -wings."""

        for deps, n_G in zip(deps_m, n_mG):
            o = abs(deps)
            w = wd.get_floor_index(o)
            if w + 1 >= wd.wmax:
                continue
            o1, o2 = wd.omega_w[w:w + 2]
            if o > o2:
                continue
            else:
                assert o1 <= o <= o2, (o1, o, o2)

            p = 1 / (o2 - o1)**2
            p1 = p * (o2 - o)
            p2 = p * (o - o1)
            x_vG = np.outer(n_G[:3], n_G.conj())
            chi0_wxvG[w, 0, :, :] += p1 * x_vG
            chi0_wxvG[w + 1, 0, :, :] += p2 * x_vG
            chi0_wxvG[w, 1, :, :] += p1 * x_vG.conj()
            chi0_wxvG[w + 1, 1, :, :] += p2 * x_vG.conj()


class TetrahedronIntegrator(Integrator):
    """Integrate brillouin zone using tetrahedron integration.

    Tetrahedron integration uses linear interpolation of
    the eigenenergies and of the matrix elements
    between the vertices of the tetrahedron."""

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

    @timer('Spectral function integration')
    def integrate(self, *, kind, domain, integrand, wd, out_wxx, task):
        """Integrate response function.

        Assume that the integral has the
        form of a response function. For the linear tetrahedron
        method it is possible calculate frequency dependent weights
        and do a point summation using these weights."""

        assert kind == task.kind, (kind, task.kind)
        assert type(task) == type(choose_tetrahedron_integral_kind(kind))

        blocks1d = self._blocks1d(out_wxx.shape[2])

        # Input domain
        td = self.tesselate(domain[0])
        args = domain[1:]

        # Relevant quantities
        bzk_kc = td.points
        nk = len(bzk_kc)

        with self.context.timer('pts'):
            # Point to simplex
            pts_k = [[] for n in range(nk)]
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
            for k in range(nk):
                pts_k[k] = np.array(pts_k[k], int)

        with self.context.timer('neighbours'):
            # Nearest neighbours
            neighbours_k = [None for n in range(nk)]

            for k in range(nk):
                neighbours_k[k] = np.unique(td.simplices[pts_k[k]])

        # Distribute everything
        myterms_t = self.distribute_domain(list(args) +
                                           [list(range(nk))])

        with self.context.timer('eigenvalues'):
            # Store eigenvalues
            deps_tMk = None  # t for term
            shape = [len(domain_l) for domain_l in args]
            nterms = int(np.prod(shape))

            for t in range(nterms):
                if len(shape) == 0:
                    arguments = ()
                else:
                    arguments = np.unravel_index(t, shape)
                for K in range(nk):
                    k_c = bzk_kc[K]
                    deps_M = -integrand.eigenvalues(k_c, *arguments)
                    if deps_tMk is None:
                        deps_tMk = np.zeros([nterms] +
                                            list(deps_M.shape) +
                                            [nk], float)
                    deps_tMk[t, :, K] = deps_M

        # Calculate integrations weight
        pb = ProgressBar(self.context.fd)
        for _, arguments in pb.enumerate(myterms_t):
            K = arguments[-1]
            if len(shape) == 0:
                t = 0
            else:
                t = np.ravel_multi_index(arguments[:-1], shape)
            deps_Mk = deps_tMk[t]
            teteps_Mk = deps_Mk[:, neighbours_k[K]]
            n_MG = integrand.matrix_element(bzk_kc[K],
                                            *arguments[:-1])

            # Generate frequency weights
            i0_M, i1_M = wd.get_index_range(
                teteps_Mk.min(1), teteps_Mk.max(1))
            W_Mw = []
            for deps_k, i0, i1 in zip(deps_Mk, i0_M, i1_M):
                W_w = self.get_kpoint_weight(K, deps_k,
                                             pts_k, wd.omega_w[i0:i1],
                                             td)
                W_Mw.append(W_w)

            task.run(n_MG, deps_Mk, W_Mw, i0_M, i1_M, out_wxx, blocks1d)

        self.kncomm.sum(out_wxx)

        if self.blockcomm.size == 1 and task.symmetrizable_unless_blocked:
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
        with self.context.timer('Tetrahedron weight'):
            _gpaw.tetrahedron_weight(deps_k, td.simplices, K,
                                     simplices_s,
                                     W_w, omega_w, vol_s)

        return W_w


class HilbertTetrahedron:
    kind = 'spectral function'
    symmetrizable_unless_blocked = True

    def run(self, n_MG, deps_Mk, W_Mw, i0_M, i1_M, out_wxx, blocks1d):
        """Update output array with dissipative part."""
        for n_G, deps_k, W_w, i0, i1 in zip(n_MG, deps_Mk, W_Mw,
                                            i0_M, i1_M):
            if i0 == i1:
                continue

            for iw, weight in enumerate(W_w):
                if blocks1d.blockcomm.size > 1:
                    myn_G = n_G[blocks1d.myslice].reshape((-1, 1))
                    # gemm(weight, n_G.reshape((-1, 1)), myn_G,
                    #      1.0, out_wxx[i0 + iw], 'c')
                    mmm(weight, myn_G, 'N', n_G.reshape((-1, 1)), 'C',
                        1.0, out_wxx[i0 + iw])
                else:
                    czher(weight, n_G.conj(), out_wxx[i0 + iw])


class HilbertOpticalLimitTetrahedron:
    kind = 'spectral function wings'
    symmetrizable_unless_blocked = False

    def run(self, n_MG, deps_Mk, W_Mw, i0_M, i1_M, out_wxvG, blocks1d):
        """Update optical limit output array with dissipative part of the head
        and wings."""
        for n_G, deps_k, W_w, i0, i1 in zip(n_MG, deps_Mk, W_Mw,
                                            i0_M, i1_M):
            assert blocks1d.blockcomm.size == 1
            if i0 == i1:
                continue

            for iw, weight in enumerate(W_w):
                x_vG = np.outer(n_G[:3], n_G.conj())
                out_wxvG[i0 + iw, 0, :, :] += weight * x_vG
                out_wxvG[i0 + iw, 1, :, :] += weight * x_vG.conj()


def choose_tetrahedron_integral_kind(kind):
    if kind == 'spectral function':
        return HilbertTetrahedron()
    if kind == 'spectral function wings':
        return HilbertOpticalLimitTetrahedron()

    raise ValueError("Expected kind='spectral function'",
                     f"or 'spectral function wings', got: {kind}")
