"""Test automagical calculation of wfs"""

import os
import sys
import time
from gpaw import GPAW
from ase import *
from gpaw.test import equal
from ase.parallel import rank, barrier, size

ending = 'gpw'
restart = 'gpaw-restart.' + ending

# H2
H = Atoms([Atom('H', (0, 0, 0)), Atom('H', (0, 0, 1))])
H.center(vacuum=2.0)

calc = GPAW(nbands=2, convergence={'eigenstates': 1e-3})
H.set_calculator(calc)
H.get_potential_energy()
calc.write(restart)

calc = GPAW(restart)
calc.set(nbands=5)
calc.converge_wave_functions()

if rank == 0:
    os.remove(restart)

