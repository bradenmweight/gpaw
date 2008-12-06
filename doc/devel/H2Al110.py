from gpaw import Calculator
from ase import *
from gpaw.mpi import rank

import time

a = 4.00
d = a / 2**0.5
z = 1.1
b = 1.5

slab = Atoms([Atom('Al', (0, 0, 0)),
                    Atom('Al', (a, 0, 0)),
                    Atom('Al', (a/2, d/2, -d/2)),
                    Atom('Al', (3*a/2, d/2, -d/2)),
                    Atom('Al', (0, 0, -d)),
                    Atom('Al', (a, 0, -d)),
                    Atom('Al', (a/2, d/2, -3*d/2)),
                    Atom('Al', (3*a/2, d/2, -3*d/2)),
                    Atom('Al', (0, 0, -2*d)),
                    Atom('Al', (a, 0, -2*d)),
                    Atom('H', (a/2-b/2, 0, z)),
                    Atom('H', (a/2+b/2, 0, z))],
                   cell=(2*a, d, 5*d), pbc=(1, 1, 1))

runs = 10
t = 0.0
for n in range(runs):
    t0 = time.time()
    for i in range(1):
        calc = Calculator(h=0.15, nbands=28, kpts=(2, 6, 1),
                  convergence={'eigenstates': 1e-5})
        slab.set_calculator(calc)
        e = slab.get_potential_energy()
        del calc
    t = t + time.time() - t0
    print 'Run: ', n, ' rank: ', str(rank), ' time: ', time.time() - t0
print 'Average time: ', t/runs, ' rank: ', str(rank)
