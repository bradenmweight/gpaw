import os
import numpy as np

from ase import Atom, Atoms
from ase.optimize import BFGS
from ase.parallel import parprint
from ase.units import Hartree
from gpaw import GPAW
from gpaw.test import equal
from gpaw.lrtddft.kssingle import KSSingles

Be = Atoms('Be')
Be.center(vacuum=6)

eigensolver = None
#eigensolver = 'rmm-diis'

#modes = ['lcao', 'fd']
modes = ['fd']

for mode in modes:
    energy = {}
    osz = {}
    for pbc in [False, True]:
        Be.set_pbc(pbc)
        if pbc:
            name = 'periodic'
            calc = GPAW(h=0.25, nbands=4, kpts=(2,2,2), mode=mode, 
                        eigensolver=eigensolver, txt=None)
        else:
            name = 'zero bc'
            calc = GPAW(h=0.25, nbands=4, mode=mode, 
                        eigensolver=eigensolver, txt=None)
        Be.set_calculator(calc)
        Be.get_potential_energy()
        
        kss = KSSingles(calc)
        # all s->p transitions at the same energy [Ha] and 
        # oscillator_strength
        for ks in kss:
            equal(ks.get_energy(), kss[0].get_energy(), 1.e-5)
            equal(ks.get_oscillator_strength()[0],
                  kss[0].get_oscillator_strength()[0], 1.e-4)
        energy[name] = np.array(
            [ks.get_energy() * Hartree for ks in kss]).mean()
        osz[name] = np.array(
            [ks.get_oscillator_strength()[0] for ks in kss]).sum()

        parprint(kss)
        # I/O
        kss.write('kss.dat')
        kss = KSSingles('kss.dat')

    # periodic and non-periodic should be roughly equal
    equal(energy['zero bc'], energy['periodic'], 1.e-2)
    equal(osz['zero bc'], osz['periodic'], 1.e-3)
