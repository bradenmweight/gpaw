import time
from math import log as ln

import numpy as np
from ase.units import Ha, Bohr

from gpaw import KohnShamConvergenceError
from gpaw.forces import calculate_forces


class SCFLoop:
    """Self-consistent field loop."""
    def __init__(self, eigenstates=0.1, energy=0.1, density=0.1, force=np.inf,
                 maxiter=100, niter_fixdensity=None, nvalence=None):
        self.max_errors = {'eigenstates': eigenstates,
                           'energy': energy,
                           'force': force,
                           'density': density}
        self.maxiter = maxiter
        self.niter_fixdensity = niter_fixdensity
        self.nvalence = nvalence

        self.old_energies = []
        self.old_F_av = None
        self.converged = False

        self.niter = None

        self.reset()

    def __str__(self):
        cc = self.max_errors
        s = 'Convergence criteria:\n'
        for name, val in [
            ('total energy change: {0:g} eV / electron',
             cc['energy'] * Ha / self.nvalence),
            ('integral of absolute density change: {0:g} electrons',
             cc['density'] / self.nvalence),
            ('integral of absolute eigenstate change: {0:g} eV^2',
             cc['eigenstates'] * Ha**2 / self.nvalence),
            ('change in atomic force: {0:g} eV / Ang',
             cc['force'] * Ha / Bohr),
            ('number of iterations: {0}', self.maxiter)]:
            if val < np.inf:
                s += '  Maximum {0}\n'.format(name.format(val))
        return s

    def write(self, writer):
        writer.write(converged=self.converged)

    def read(self, reader):
        self.converged = reader.scf.converged

    def reset(self):
        self.old_energies = []
        self.old_F_av = None
        self.converged = False

    def run(self, wfs, ham, dens, occ, log, callback):

        egs_name = getattr(wfs.eigensolver, "name", None)
        if egs_name == 'direct_min':
            self.run_dm(wfs, ham, dens, occ, log, callback)
            return

        self.niter = 1
        while self.niter <= self.maxiter:
            wfs.eigensolver.iterate(ham, wfs, occ)
            occ.calculate(wfs)
            energy = ham.get_energy(occ)
            self.old_energies.append(energy)
            errors = self.collect_errors(dens, ham, wfs)

            # Converged?
            for kind, error in errors.items():
                if error > self.max_errors[kind]:
                    self.converged = False
                    break
            else:
                self.converged = True

            callback(self.niter)
            self.log(log, self.niter, wfs, ham, dens, occ, errors)

            if self.converged and self.niter >= self.niter_fixdensity:
                break

            if self.niter > self.niter_fixdensity and not dens.fixed:
                dens.update(wfs)
                ham.update(dens)
            else:
                ham.npoisson = 0
            self.niter += 1

        # Don't fix the density in the next step:
        self.niter_fixdensity = 0

        if not self.converged:
            if not np.isfinite(errors['eigenstates']):
                msg = 'Not enough bands for ' + wfs.eigensolver.nbands_converge
                log(msg)
                raise KohnShamConvergenceError(msg)
            log(oops)
            raise KohnShamConvergenceError(
                'Did not converge!  See text output for help.')

    def collect_errors(self, dens, ham, wfs):
        """Check convergence of eigenstates, energy and density."""

        errors = {'eigenstates': wfs.eigensolver.error,
                  'density': dens.error,
                  'energy': np.inf}

        if dens.fixed:
            errors['density'] = 0.0

        if len(self.old_energies) >= 3:
            energies = self.old_energies[-3:]
            if np.isfinite(energies).all():
                errors['energy'] = np.ptp(energies)

        # We only want to calculate the (expensive) forces if we have to:
        check_forces = (self.max_errors['force'] < np.inf and
                        all(error <= self.max_errors[kind]
                            for kind, error in errors.items()))

        errors['force'] = np.inf
        if check_forces:
            with wfs.timer('Forces'):
                F_av = calculate_forces(wfs, dens, ham)
            if self.old_F_av is not None:
                errors['force'] = ((F_av - self.old_F_av)**2).sum(1).max()**0.5
            self.old_F_av = F_av

        return errors

    def log(self, log, niter, wfs, ham, dens, occ, errors):
        """Output from each iteration."""

        nvalence = wfs.nvalence
        if nvalence > 0:
            eigerr = errors['eigenstates'] * Ha**2 / nvalence
        else:
            eigerr = 0.0

        T = time.localtime()

        if niter == 1:
            header = """\
                     log10-error:    total        iterations:
           time      wfs    density  energy       fermi  poisson"""
            if wfs.nspins == 2:
                header += '  magmom'
            if self.max_errors['force'] < np.inf:
                l1 = header.find('total')
                header = header[:l1] + '       ' + header[l1:]
                l2 = header.find('energy')
                header = header[:l2] + 'force  ' + header[l2:]
            log(header)

        if eigerr == 0.0 or np.isinf(eigerr):
            eigerr = ''
        else:
            eigerr = '%+.2f' % (ln(eigerr) / ln(10))

        denserr = errors['density']
        assert denserr is not None
        if (denserr is None or np.isinf(denserr) or denserr == 0 or
            nvalence == 0):
            denserr = ''
        else:
            denserr = '%+.2f' % (ln(denserr / nvalence) / ln(10))

        niterocc = occ.niter
        if niterocc == -1:
            niterocc = ''
        else:
            niterocc = '%d' % niterocc

        if ham.npoisson == 0:
            niterpoisson = ''
        else:
            niterpoisson = str(ham.npoisson)

        log('iter: %3d  %02d:%02d:%02d %6s %6s  ' %
            (niter,
             T[3], T[4], T[5],
             eigerr,
             denserr), end='')

        if self.max_errors['force'] < np.inf:
            if errors['force'] == 0:
                log('    -oo', end='')
            elif errors['force'] < np.inf:
                log('  %+.2f' %
                    (ln(errors['force'] * Ha / Bohr) / ln(10)), end='')
            else:
                log('       ', end='')

        if np.isfinite(ham.e_total_extrapolated):
            energy = '{:11.6f}'.format(Ha * ham.e_total_extrapolated)
        else:
            energy = ' ' * 11

        log('%s    %-5s  %-7s' %
            (energy, niterocc, niterpoisson), end='')

        if wfs.nspins == 2:
            log('  %+.4f' % occ.magmom, end='')
        elif not wfs.collinear:
            totmom_v, magmom_av = dens.estimate_magnetic_moments()
            log(' {:+.1f},{:+.1f},{:+.1f}'.format(*totmom_v), end='')

        log(flush=True)

    def run_dm(self, wfs, ham, dens, occ, log, callback):

        self.niter = 1

        wfs.eigensolver.eg_count = 0

        while self.niter <= self.maxiter:
            wfs.eigensolver.iterate(ham, wfs, dens, occ, log)
            if hasattr(wfs.eigensolver, 'e_sic'):
                e_sic = wfs.eigensolver.e_sic
            else:
                e_sic = 0.0
            energy = ham.get_energy(occ, kin_en_using_band=False,
                                    e_sic=e_sic)

            self.old_energies.append(energy)
            errors = self.collect_errors(dens, ham, wfs)

            # Converged?
            for kind, error in errors.items():
                if error > self.max_errors[kind]:
                    self.converged = False
                    break
            else:
                self.converged = True

            callback(self.niter)
            self.log(log, self.niter, wfs, ham, dens, occ, errors)

            if self.converged:
                # lcao
                occ.calculate(wfs)
                wfs.eigensolver.get_canonical_representation(ham,
                                                             wfs)
                niter = wfs.eigensolver.eg_count
                log(
                    '\nOccupied states converged after {:d} e/g evaluations'.format(
                        niter))
                break

            ham.npoisson = 0
            self.niter += 1

        if not self.converged:
            log(oops)
            raise KohnShamConvergenceError(
                'Did not converge!  See text output for help.')

oops = """
Did not converge!

Here are some tips:

1) Make sure the geometry and spin-state is physically sound.
2) Use less aggressive density mixing.
3) Solve the eigenvalue problem more accurately at each scf-step.
4) Use a smoother distribution function for the occupation numbers.
5) Try adding more empty states.
6) Use enough k-points.
7) Don't let your structure optimization algorithm take too large steps.
8) Solve the Poisson equation more accurately.
9) Better initial guess for the wave functions.

See details here:

    https://wiki.fysik.dtu.dk/gpaw/documentation/convergence.html

"""
