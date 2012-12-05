from gpaw.poisson import PoissonSolver
from gpaw.transformers import Transformer
from gpaw.fd_operators import Laplace, Gradient
from gpaw.wfd_operators import WeightedFDOperator
from gpaw.utilities.gauss import Gaussian
from gpaw.utilities import erf
import numpy


class SolvationPoissonSolver(PoissonSolver):
    def __init__(self, nn=3, relax='J', eps=2e-10, dielectric=None):
        """
        initialize the Poisson solver

        dielectric -- list [epsr dx_epsr, dy_epsr, dz_epsr]
        """
        self.dielectric = dielectric
        PoissonSolver.__init__(self, nn, relax, eps)

    def load_gauss(self):
        # XXX Check if update is needed (dielectric changed)?
        epsr, dx_epsr, dy_epsr, dz_epsr = self.dielectric
        gauss = Gaussian(self.gd)
        rho_g = gauss.get_gauss(0)
        phi_g = gauss.get_gauss_pot(0)
        x, y, z = gauss.xyz
        fac = 2. * numpy.sqrt(gauss.a) * numpy.exp(-gauss.a * gauss.r2)
        fac /= numpy.sqrt(numpy.pi) * gauss.r2
        fac -= erf(numpy.sqrt(gauss.a) * gauss.r) / (gauss.r2 * gauss.r)
        fac *= 2.0 * 1.7724538509055159
        dx_phi_g = fac * x
        dy_phi_g = fac * y
        dz_phi_g = fac * z
        sp = dx_phi_g * dx_epsr + dy_phi_g * dy_epsr + dz_phi_g * dz_epsr
        rho = epsr * rho_g - 1. / (4. * numpy.pi) * sp
        invnorm = numpy.sqrt(4. * numpy.pi) / self.gd.integrate(rho)
        self.phi_gauss = phi_g * invnorm
        self.rho_gauss = rho * invnorm


class WeightedFDPoissonSolver(SolvationPoissonSolver):
    """
    Poisson solver including an electrostatic solvation model

    following Sanchez et al J. Chem. Phys. 131 (2009) 174108
    """

    def solve(self, phi, rho, charge=None, eps=None,
              maxcharge=1e-6,
              zero_initial_phi=False):
        self.restrict_op_weights()
        ret = PoissonSolver.solve(self, phi, rho, charge, eps, maxcharge,
                                  zero_initial_phi)
        return ret

    def restrict_op_weights(self):
        weights = [self.dielectric] + self.op_coarse_weights
        for i, res in enumerate(self.restrictors):
            for j in xrange(4):
                res.apply(weights[i][j], weights[i + 1][j])
        self.step = 0.66666666 / self.operators[0].get_diagonal_element()

    def set_grid_descriptor(self, gd):
        if gd.pbc_c.any():
            raise NotImplementedError(
                'WeightedFDPoissonSolver supports only '
                'non-periodic boundary conditions up to now.'
                )
        self.gd = gd
        self.gds = [gd]
        self.dv = gd.dv
        gd = self.gd
        self.B = None
        self.interpolators = []
        self.restrictors = []
        self.operators = []
        level = 0
        self.presmooths = [2]
        self.postsmooths = [1]
        self.weights = [2. / 3.]
        while level < 4:
            try:
                gd2 = gd.coarsen()
            except ValueError:
                break
            self.gds.append(gd2)
            self.interpolators.append(Transformer(gd2, gd))
            self.restrictors.append(Transformer(gd, gd2))
            self.presmooths.append(4)
            self.postsmooths.append(4)
            self.weights.append(1.0)
            level += 1
            gd = gd2
        self.levels = level

    def initialize(self, load_gauss=False):
        self.presmooths[self.levels] = 8
        self.postsmooths[self.levels] = 8
        self.phis = [None] + [gd.zeros() for gd in self.gds[1:]]
        self.residuals = [gd.zeros() for gd in self.gds]
        self.rhos = [gd.zeros() for gd in self.gds]
        self.op_coarse_weights = [[g.empty() for g in (gd, ) * 4] \
                               for gd in self.gds[1:]]
        scale = -0.25 / numpy.pi
        for i, gd in enumerate(self.gds):
            if i == 0:
                nn = self.nn
                weights = self.dielectric
            else:
                nn = 1
                weights = self.op_coarse_weights[i - 1]
            operators = [Laplace(gd, scale, nn)] + \
                        [Gradient(gd, j, scale, nn) for j in (0, 1, 2)]
            self.operators.append(WeightedFDOperator(weights, operators))
        if load_gauss:
            self.load_gauss()
        if self.relax_method == 1:
            self.description = 'Gauss-Seidel'
        else:
            self.description = 'Jacobi'
        self.description += ' solver with dielectric and ' \
                            '%d multi-grid levels' % (self.levels + 1, )
        self.description += '\nStencil: ' + self.operators[0].description


class PolarizationPoissonSolver(SolvationPoissonSolver):
    """
    Poisson solver including an electrostatic solvation model

    calculates the polarization charges first using only the
    vacuum poisson equation, then solves the vacuum equation
    with polarization charges
    """

    def __init__(self, nn=3, relax='J', eps=2e-10, dielectric=None):
        SolvationPoissonSolver.__init__(self, nn, relax, eps, dielectric)
        self.phi_tilde = None

    def solve(self, phi, rho, charge=None, eps=None,
              maxcharge=1e-6,
              zero_initial_phi=False):
        if self.phi_tilde is None:
            self.phi_tilde = self.gd.zeros()
        phi_tilde = self.phi_tilde
        niter_tilde = PoissonSolver.solve(
            self, phi_tilde, rho, None, self.eps,
            maxcharge, False
            )

        epsr, dx_epsr, dy_epsr, dz_epsr = self.dielectric
        dx_phi_tilde = self.gd.empty()
        dy_phi_tilde = self.gd.empty()
        dz_phi_tilde = self.gd.empty()
        Gradient(self.gd, 0, 1.0, 3).apply(phi_tilde, dx_phi_tilde)
        Gradient(self.gd, 1, 1.0, 3).apply(phi_tilde, dy_phi_tilde)
        Gradient(self.gd, 2, 1.0, 3).apply(phi_tilde, dz_phi_tilde)

        scalar_product = dx_epsr * dx_phi_tilde + \
                         dy_epsr * dy_phi_tilde + \
                         dz_epsr * dz_phi_tilde

        rho_and_pol = rho / epsr + \
                      scalar_product / (4. * numpy.pi * epsr ** 2)

        niter = PoissonSolver.solve(
            self, phi, rho_and_pol, None, eps,
            maxcharge, zero_initial_phi
            )
        return niter_tilde + niter

    def load_gauss(self):
        return PoissonSolver.load_gauss(self)


class IterativePoissonSolver(SolvationPoissonSolver):
    """
    Poisson solver including an electrostatic solvation model

    following Andreussi et al.
    The Journal of Chemical Physics 136, 064102 (2012)
    """
    # XXX broken convergence for eta != 1.0 (non-self-consistent) ???
    eta = 1.0

    def set_grid_descriptor(self, gd):
        SolvationPoissonSolver.set_grid_descriptor(self, gd)
        self.dx_phi = gd.empty()
        self.dy_phi = gd.empty()
        self.dz_phi = gd.empty()
        self.gradx = Gradient(gd, 0, 1.0, 3)
        self.grady = Gradient(gd, 1, 1.0, 3)
        self.gradz = Gradient(gd, 2, 1.0, 3)

    def solve_neutral(self, phi, rho, eps=2e-10):
        self.rho_iter = self.gd.zeros()
        self.rho = rho
        return SolvationPoissonSolver.solve_neutral(self, phi, rho, eps)

    def iterate2(self, step, level=0):
        if level == 0:
            epsr, dx_epsr, dy_epsr, dz_epsr = self.dielectric
            self.gradx.apply(self.phis[0], self.dx_phi)
            self.grady.apply(self.phis[0], self.dy_phi)
            self.gradz.apply(self.phis[0], self.dz_phi)
            sp = dx_epsr * self.dx_phi + \
                 dy_epsr * self.dy_phi + \
                 dz_epsr * self.dz_phi
            if self.eta == 1.0:
                self.rho_iter = 1. / (4. * numpy.pi) * sp
            else:
                self.rho_iter = self.eta / (4. * numpy.pi) * sp + \
                                (self.eta - 1.) * self.rho_iter
            self.rhos[0][:] = (self.rho_iter + self.rho) / epsr
        return SolvationPoissonSolver.iterate2(self, step, level)
