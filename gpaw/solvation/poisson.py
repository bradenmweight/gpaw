from gpaw.poisson import PoissonSolver
from gpaw.transformers import Transformer
from gpaw.fd_operators import Laplace, Gradient
from gpaw.wfd_operators import WeightedFDOperator
from gpaw.utilities.gauss import Gaussian
from gpaw.utilities import erf
import warnings
import numpy as np


class SolvationPoissonSolver(PoissonSolver):
    def set_dielectric(self, dielectric):
        self.dielectric = dielectric

    def load_gauss(self):
        # XXX Check if update is needed (dielectric changed)?
        epsr, dx_epsr, dy_epsr, dz_epsr = self.dielectric.eps_gradeps
        gauss = Gaussian(self.gd)
        rho_g = gauss.get_gauss(0)
        phi_g = gauss.get_gauss_pot(0)
        x, y, z = gauss.xyz
        fac = 2. * np.sqrt(gauss.a) * np.exp(-gauss.a * gauss.r2)
        fac /= np.sqrt(np.pi) * gauss.r2
        fac -= erf(np.sqrt(gauss.a) * gauss.r) / (gauss.r2 * gauss.r)
        fac *= 2.0 * 1.7724538509055159
        dx_phi_g = fac * x
        dy_phi_g = fac * y
        dz_phi_g = fac * z
        sp = dx_phi_g * dx_epsr + dy_phi_g * dy_epsr + dz_phi_g * dz_epsr
        rho = epsr * rho_g - 1. / (4. * np.pi) * sp
        invnorm = np.sqrt(4. * np.pi) / self.gd.integrate(rho)
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
        if self.gd.pbc_c.all():
            actual_charge = self.gd.integrate(rho)
            if abs(actual_charge) > maxcharge:
                raise NotImplementedError(
                    'charged periodic systems are not implemented'
                    )
        self.restrict_op_weights()
        ret = PoissonSolver.solve(self, phi, rho, charge, eps, maxcharge,
                                  zero_initial_phi)
        return ret

    def restrict_op_weights(self):
        weights = [self.dielectric.eps_gradeps] + self.op_coarse_weights
        for i, res in enumerate(self.restrictors):
            for j in xrange(4):
                res.apply(weights[i][j], weights[i + 1][j])
        self.step = 0.66666666 / self.operators[0].get_diagonal_element()

    def set_grid_descriptor(self, gd):
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
        scale = -0.25 / np.pi
        for i, gd in enumerate(self.gds):
            if i == 0:
                nn = self.nn
                weights = self.dielectric.eps_gradeps
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
    """Poisson solver including an electrostatic solvation model

    calculates the polarization charges first using only the
    vacuum poisson equation, then solves the vacuum equation
    with polarization charges

    warning: not exact enough, since the electric field is not exact enough

    """

    def __init__(self, nn=3, relax='J', eps=2e-10):
        polarization_warning = UserWarning(
            (
                'PolarizationPoissonSolver is not accurate enough'
                ' and therefore not recommended for production code!'
                )
            )
        warnings.warn(polarization_warning)
        SolvationPoissonSolver.__init__(self, nn, relax, eps)
        self.phi_tilde = None

    def set_grid_descriptor(self, gd):
        SolvationPoissonSolver.set_grid_descriptor(self, gd)
        if self.relax_method == 1:
            self.description = 'Polarization Gauss-Seidel'
        else:
            self.description = 'Polarization Jacobi'
        self.description += ' solver with %d multi-grid levels' % (
            self.levels + 1,
            )
        self.description += '\nStencil: ' + self.operators[0].description

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

        epsr, dx_epsr, dy_epsr, dz_epsr = self.dielectric.eps_gradeps
        dx_phi_tilde = self.gd.empty()
        dy_phi_tilde = self.gd.empty()
        dz_phi_tilde = self.gd.empty()
        Gradient(self.gd, 0, 1.0, self.nn).apply(phi_tilde, dx_phi_tilde)
        Gradient(self.gd, 1, 1.0, self.nn).apply(phi_tilde, dy_phi_tilde)
        Gradient(self.gd, 2, 1.0, self.nn).apply(phi_tilde, dz_phi_tilde)

        scalar_product = dx_epsr * dx_phi_tilde + \
                         dy_epsr * dy_phi_tilde + \
                         dz_epsr * dz_phi_tilde

        rho_and_pol = rho / epsr + \
                      scalar_product / (4. * np.pi * epsr ** 2)

        niter = PoissonSolver.solve(
            self, phi, rho_and_pol, None, eps,
            maxcharge, zero_initial_phi
            )
        return niter_tilde + niter

    def load_gauss(self):
        return PoissonSolver.load_gauss(self)


class ADM12PoissonSolver(SolvationPoissonSolver):
    """Poisson solver including an electrostatic solvation model

    following Andreussi et al.
    The Journal of Chemical Physics 136, 064102 (2012)

    XXX TODO : * Correction for charged systems???
               * Check: Can the polarization charge introduce a monopole?
               * Convergence problems depending on eta, apparently this
                 method works best with FFT as in the original Paper
               * Optimize numerics
    """

    def __init__(self, nn=3, relax='J', eps=2e-10, eta=.6):
        adm12_warning = UserWarning(
            (
                'ADM12PoissonSolver is not tested thoroughly'
                ' and therefore not recommended for production code!'
                )
            )
        warnings.warn(adm12_warning)
        self.eta = eta
        SolvationPoissonSolver.__init__(self, nn, relax, eps)

    def set_grid_descriptor(self, gd):
        SolvationPoissonSolver.set_grid_descriptor(self, gd)
        self.gradx = Gradient(gd, 0, 1.0, self.nn)
        self.grady = Gradient(gd, 1, 1.0, self.nn)
        self.gradz = Gradient(gd, 2, 1.0, self.nn)
        if self.relax_method == 1:
            self.description = 'ADM12 Gauss-Seidel'
        else:
            self.description = 'ADM12 Jacobi'
        self.description += ' solver with %d multi-grid levels' % (
            self.levels + 1,
            )
        self.description += '\nStencil: ' + self.operators[0].description

    def initialize(self, load_gauss=False):
        self.rho_iter = self.gd.zeros()
        self.d_phi = self.gd.empty()
        return SolvationPoissonSolver.initialize(self, load_gauss)

    def solve(self, phi, rho, charge=None, eps=None,
              maxcharge=1e-6,
              zero_initial_phi=False):
        if self.gd.pbc_c.all():
            actual_charge = self.gd.integrate(rho)
            if abs(actual_charge) > maxcharge:
                raise NotImplementedError(
                    'charged periodic systems are not implemented'
                    )
        return PoissonSolver.solve(
            self, phi, rho, charge, eps, maxcharge, zero_initial_phi
            )

    def solve_neutral(self, phi, rho, eps=2e-10):
        self.rho = rho
        return SolvationPoissonSolver.solve_neutral(self, phi, rho, eps)

    def iterate2(self, step, level=0):
        if level == 0:
            epsr, dx_epsr, dy_epsr, dz_epsr = self.dielectric.eps_gradeps
            self.gradx.apply(self.phis[0], self.d_phi)
            sp = dx_epsr * self.d_phi
            self.grady.apply(self.phis[0], self.d_phi)
            sp += dy_epsr * self.d_phi
            self.gradz.apply(self.phis[0], self.d_phi)
            sp += dz_epsr * self.d_phi
            self.rho_iter = self.eta / (4. * np.pi) * sp + \
                (1. - self.eta) * self.rho_iter
            self.rhos[0][:] = (self.rho_iter + self.rho) / epsr
        return SolvationPoissonSolver.iterate2(self, step, level)
