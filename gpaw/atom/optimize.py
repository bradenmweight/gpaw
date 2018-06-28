from __future__ import print_function, division
import os
import re
import sys
import time

import numpy as np
from scipy.optimize import differential_evolution as DE
from ase import Atoms
from ase.data import covalent_radii, atomic_numbers
from ase.units import Bohr

from gpaw import GPAW, PW, setup_paths, Mixer, ConvergenceError, Davidson
from gpaw.atom.generator2 import _generate  # , DatasetGenerationError


my_covalent_radii = covalent_radii.copy()
for e in ['Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr']:  # missing radii
    my_covalent_radii[atomic_numbers[e]] = 1.7


class PAWDataError(Exception):
    """Error in PAW-data generation."""


class DatasetOptimizer:
    tolerances = np.array([0.2,
                           0.3,
                           40,
                           2 / 3 * 300 * 0.1**0.25,  # convergence
                           0.0005])  # eggbox error

    def __init__(self, symbol='H', nc=False, processes=None):
        self.old = False

        self.symbol = symbol
        self.nc = nc
        self.processes = processes

        with open('../start.txt') as fd:
            for line in fd:
                words = line.split()
                if words[1] == symbol:
                    projectors = words[3]
                    radii = [float(f) for f in words[5].split(',')]
                    r0 = float(words[7].split(',')[1])
                    break
            else:
                raise ValueError

        self.Z = atomic_numbers[symbol]
        self.rc = my_covalent_radii[self.Z]

        # Parse projectors string:
        pattern = r'(-?\d+\.\d)'
        energies = []
        for m in re.finditer(pattern, projectors):
            energies.append(float(projectors[m.start():m.end()]))
        self.projectors = re.sub(pattern, '{:.1f}', projectors)
        self.nenergies = len(energies)

        # Round to integers:
        self.x = energies + radii + [r0]
        self.bounds = ([(-1.0, 4.0)] * self.nenergies +
                       [(r * 0.6, r * 1.5) for r in radii] +
                       [(0.3, max(radii) * 1.4)])

        self.ecut1 = 450.0
        self.ecut2 = 800.0

        setup_paths[:0] = ['.']

        self.logfile = None
        self.tflush = time.time() + 60

    def run(self):
        print(self.symbol, self.rc / Bohr, self.projectors)
        print(self.x)
        print(self.bounds)
        self.logfile = open('data.csv', 'a')
        DE(self, self.bounds)

    def generate(self, fd, projectors, radii, r0, xc,
                 scalar_relativistic=True, tag=None, logderivs=True):

        if projectors[-1].isupper():
            nderiv0 = 5
        else:
            nderiv0 = 2

        type = 'poly'
        if self.nc:
            type = 'nc'

        try:
            gen = _generate(self.symbol, xc, None, projectors, radii,
                            scalar_relativistic, None, r0, nderiv0,
                            (type, 4), None, None, fd)
        except np.linalg.LinAlgError:
            raise PAWDataError('LinAlgError')

        if not scalar_relativistic:
            if not gen.check_all():
                raise PAWDataError('dataset check failed')

        if tag is not None:
            gen.make_paw_setup(tag or None).write_xml()

        r = 1.1 * gen.rcmax

        lmax = 2
        if 'f' in projectors:
            lmax = 3

        error = 0.0
        if logderivs:
            for l in range(lmax + 1):
                emin = -1.5
                emax = 2.0
                n0 = gen.number_of_core_states(l)
                if n0 > 0:
                    e0_n = gen.aea.channels[l].e_n
                    emin = max(emin, e0_n[n0 - 1] + 0.1)
                energies = np.linspace(emin, emax, 100)
                de = energies[1] - energies[0]
                ld1 = gen.aea.logarithmic_derivative(l, energies, r)
                ld2 = gen.logarithmic_derivative(l, energies, r)
                error += abs(ld1 - ld2).sum() * de

        return error

    def parameters(self, x):
        energies = x[:self.nenergies]
        radii = x[self.nenergies:-1]
        r0 = x[-1]
        projectors = self.projectors.format(*energies)
        return energies, radii, r0, projectors

    def __call__(self, x):
        energies, radii, r0, projectors = self.parameters(x)

        print('({}, {}, {:.2f}):'
              .format(', '.join('{:+.2f}'.format(e) for e in energies),
                      ', '.join('{:.2f}'.format(r) for r in radii),
                      r0),
              end='')

        fd = open('out.txt', 'w')
        errors, msg, convenergies, eggenergies = self.test(fd, projectors,
                                                           radii, r0)
        error = ((errors / self.tolerances)**2).sum()

        print('{:9.1f} ({:.2f}, {:.3f}, {:.1f}, {}, {:.5f}) {}'
              .format(error, *errors, msg))

        convenergies += [0] * (7 - len(convenergies))

        print(', '.join(repr(number) for number in
                        list(x) + [error] + errors +
                        convenergies + eggenergies),
              file=self.logfile)

        if time.time() > self.tflush:
            sys.stdout.flush()
            self.logfile.flush()
            self.tflush = time.time() + 60

        return error

    def test_old_paw_data(self):
        fd = open('old.txt', 'w')
        area, niter, convenergies = self.convergence(fd, 'paw')
        eggenergies = self.eggbox(fd, 'paw')
        print('RESULTS:',
              ', '.join(repr(number) for number in
                        [area, niter, max(eggenergies)] +
                        convenergies + eggenergies),
              file=fd)

    def test(self, fd, projectors, radii, r0):
        errors = [np.inf] * 5
        energies = []
        eggenergies = [0, 0, 0]
        msg = ''

        try:
            if any(r < r0 for r in radii):
                raise PAWDataError('Core radius too large')

            rc = self.rc / Bohr
            errors[0] = sum(r - rc for r in radii if r > rc)

            error = 0.0
            for kwargs in [dict(xc='PBE', tag='de'),
                           dict(xc='PBE', scalar_relativistic=False),
                           dict(xc='LDA'),
                           dict(xc='PBEsol'),
                           dict(xc='RPBE'),
                           dict(xc='PW91')]:
                error += self.generate(fd, projectors, radii, r0, **kwargs)
            errors[1] = error

            area, niter, energies = self.convergence(fd)
            errors[2] = area
            errors[3] = niter

            eggenergies = self.eggbox(fd)
            errors[4] = max(eggenergies)

        except (ConvergenceError, PAWDataError, RuntimeError,
                np.linalg.LinAlgError) as e:
            msg = str(e)

        return errors, msg, energies, eggenergies

    def eggbox(self, fd, setup='de'):
        energies = []
        for h in [0.16, 0.18, 0.2]:
            a0 = 16 * h
            atoms = Atoms(self.symbol, cell=(a0, a0, 2 * a0), pbc=True)
            if 58 <= self.Z <= 70 or 90 <= self.Z <= 102:
                M = 999
                mixer = {'mixer': Mixer(0.01, 5)}
            else:
                M = 333
                mixer = {}
            atoms.calc = GPAW(h=h,
                              eigensolver=Davidson(niter=2),
                              xc='PBE',
                              symmetry='off',
                              setups=setup,
                              maxiter=M,
                              txt=fd,
                              **mixer)
            atoms.positions += h / 2  # start with broken symmetry
            e0 = atoms.get_potential_energy()
            atoms.positions -= h / 6
            e1 = atoms.get_potential_energy()
            atoms.positions -= h / 6
            e2 = atoms.get_potential_energy()
            atoms.positions -= h / 6
            e3 = atoms.get_potential_energy()
            energies.append(np.ptp([e0, e1, e2, e3]))
        # print(energies)
        return energies

    def convergence(self, fd, setup='de'):
        a = 3.0
        atoms = Atoms(self.symbol, cell=(a, a, a), pbc=True)
        if 58 <= self.Z <= 70 or 90 <= self.Z <= 102:
            M = 999
            mixer = {'mixer': Mixer(0.01, 5)}
        else:
            M = 333
            mixer = {}
        atoms.calc = GPAW(mode=PW(1500),
                          xc='PBE',
                          setups=setup,
                          symmetry='off',
                          maxiter=M,
                          txt=fd,
                          **mixer)
        e0 = atoms.get_potential_energy()
        energies = [e0]
        iters = atoms.calc.get_number_of_iterations()
        oldfde = None
        area = 0.0

        def f(x):
            return x**0.25

        for ec in range(800, 200, -100):
            atoms.calc.set(mode=PW(ec))
            atoms.calc.set(eigensolver='rmm-diis')
            de = atoms.get_potential_energy() - e0
            energies.append(de)
            # print(ec, de)
            fde = f(abs(de))
            if fde > f(0.1):
                if oldfde is None:
                    return np.inf, iters, energies
                ec0 = ec + (fde - f(0.1)) / (fde - oldfde) * 100
                area += ((ec + 100) - ec0) * (f(0.1) + oldfde) / 2
                break

            if oldfde is not None:
                area += 100 * (fde + oldfde) / 2
            oldfde = fde

        return area, iters, energies

    def read(self):
        data = np.loadtxt('data.csv', delimiter=',')
        return data[data[:, len(self.x)].argsort()]

    def summary(self, N=10):
        n = len(self.x)
        for x in self.read()[:N]:
            print('{:3} {:2} {:6.1f} ({}) ({}, {})'
                  .format(self.Z,
                          self.symbol,
                          x[n],
                          ', '.join('{:4.1f}'.format(e) + s
                                    for e, s
                                    in zip(x[n + 1:n + 6] / self.tolerances,
                                           'rlcie')),
                          ', '.join('{:+.2f}'.format(e)
                                    for e in x[:self.nenergies]),
                          ', '.join('{:.2f}'.format(r)
                                    for r in x[self.nenergies:n])))

    def best(self):
        n = len(self.x)
        a = self.read()[0]
        x = a[:n]
        error = a[n]
        energies, radii, r0, projectors = self.parameters(x)
        if 1:
            if projectors[-1].isupper():
                nderiv0 = 5
            else:
                nderiv0 = 2
            fmt = '{0:2} {1:2} -P {2:31} -r {3:20} -0 {4},{5:.2f} # {6:10.3f}'
            print(fmt.format(self.Z,
                             self.symbol,
                             projectors,
                             ','.join('{0:.2f}'.format(r) for r in radii),
                             nderiv0,
                             r0,
                             error))
        if 1 and error != np.inf and error != np.nan:
            self.generate(None, projectors, radii, r0, 'PBE', True, 'a1',
                          logderivs=False)


if __name__ == '__main__':
    import argparse
    from pathlib import Path
    parser = argparse.ArgumentParser(usage='python -m gpaw.atom.optimize '
                                     '[options] [folder folder ...]',
                                     description='Optimize PAW data')
    parser.add_argument('-s', '--summary', type=int)
    parser.add_argument('-b', '--best', action='store_true')
    parser.add_argument('-n', '--norm-conserving', action='store_true')
    parser.add_argument('-o', '--old-setups', action='store_true')
    parser.add_argument('folder', nargs='*')
    args = parser.parse_args()
    folders = [Path(folder) for folder in args.folder or ['.']]
    home = Path.cwd()
    for folder in folders:
        try:
            os.chdir(folder)
            symbol = Path.cwd().name
            do = DatasetOptimizer(symbol)
            if args.summary:
                do.summary(args.summary)
            elif args.old_setups:
                do.test_old_paw_data()
            elif args.best:
                do.best()
            else:
                do.run()
        finally:
            os.chdir(home)
            