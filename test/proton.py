#!/usr/bin/env python
from ase import *
from gpaw import Calculator

a = 4.0
H = Atoms([Atom('H', (a/2, a/2, a/2), magmom=1)],
                pbc=0,
                cell=(a, a, a))
calc = Calculator(nbands=1, h=0.2, charge=1)
H.set_calculator(calc)
print H.get_potential_energy() + calc.get_reference_energy()
assert abs(H.get_potential_energy() + calc.get_reference_energy()) < 0.014
