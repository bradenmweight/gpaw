from __future__ import print_function
# This tests calculates the force on the atoms of a
# slightly distorted Silicon primitive cell.
#
# If the test fails, set the fd boolean below to enable a (costly) finite
# difference check.

import numpy as np
from ase import Atoms
from gpaw import GPAW
from gpaw.atom.basis import BasisMaker
from gpaw.test import equal


def test_generic_guc_force():
    sibasis = BasisMaker('Si').generate(2, 1, energysplit=0.3, tailnorm=0.03**.5)
    basis = {'Si' : sibasis}

    a = 5.475
    system = Atoms(symbols='Si2', pbc=True,
                   cell=0.5 * a * np.array([(1, 1, 0),
                                            (1, 0, 1),
                                            (0, 1, 1)]),
                   scaled_positions=[(0.0, 0.0, 0.0),
                                     (0.23, 0.23, 0.23)])

    calc = GPAW(h=0.2,
                mode='lcao',
                basis=basis,
                kpts=(2,2,2),
                convergence={'density':1e-5, 'energy': 1e-6}
                )
    system.set_calculator(calc)

    F_ac = system.get_forces()

    # Previous FD result, generated by disabled code below
    F_ac_ref = np.array(
        [[-1.3967114867039498, -1.3967115816022613, -1.396711581510779],
         [1.397400325299003, 1.3974003410455182, 1.3974003410801572]])

    err_ac = np.abs(F_ac - F_ac_ref)
    err = err_ac.max()

    print('Force')
    print(F_ac)
    print()
    print('Reference result')
    print(F_ac_ref)
    print()
    print('Error')
    print(err_ac)
    print()
    print('Max error')
    print(err)

    # ASE uses dx = [+|-] 0.001 by default,
    # error should be around 2e-3.  In fact 4e-3 would probably be acceptable


    equal(err, 0, 6e-3)

    # Set boolean to run new FD check
    fd = False

    if fd:
        from ase.calculators.test import numeric_force
        calc.set(symmetry='off')
        F_ac_fd = [[numeric_force(system, a, i) for i in range(3)]
                   for a in range(2)]
        print('Self-consistent forces')
        print(F_ac)
        print('FD')
        print(F_ac_fd)
        print(repr(F_ac_fd))
        print(F_ac - F_ac_fd, np.abs(F_ac - F_ac_fd).max())

        assert np.abs(F_ac - F_ac_fd).max() < 4e-3
