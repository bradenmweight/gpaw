import numpy as np
import pytest

from ase import Atoms
from ase.units import Ha

from gpaw import GPAW
from gpaw.scf import WorkFunction, Energy, Criterion


class FourIterations(Criterion):
    """A silly custom convergence criterion that ensures it runs
    at least four iterations."""
    name = 'four iterations'
    tablename = 'four'
    calc_last = False

    def __init__(self):
        self.description = 'At least four iterations must complete!'
        self.reset()

    def __call__(self, context):
        converged = self.iters >= 4
        entry = str(self.iters)
        self.iters += 1
        return converged, entry

    def reset(self):
        self.iters = 0


def test_scf_criterion(in_tmp_dir):
    """Tests different ways of setting SCF convergence criteria,
    and that it behaves consistenly with regard to the work function."""
    convergence = {'eigenstates': 1.0,
                   'density': 1.0,
                   'energy': 1.0,
                   'work function': 1.0}

    atoms = Atoms('HF', [(0., 0.5, 0.5),
                         (0., 0.4, -0.4)],
                  cell=(5., 5., 9.),
                  pbc=(True, True, False))
    atoms.center()
    calc = GPAW(h=0.3,
                nbands=-1,
                convergence=convergence,
                txt=None,
                poissonsolver={'dipolelayer': 'xy'})
    atoms.calc = calc
    atoms.get_potential_energy()
    fermilevel = calc.wfs.fermi_level
    workfunctions1 = Ha * calc.hamiltonian.get_workfunctions(fermilevel)
    calc.write('scf-criterion.gpw')

    # Flip and use saved calculator; work functions should be opposite.
    atoms = Atoms('HF', [(0., 0.5, -0.5),
                         (0., 0.4, +0.4)],
                  cell=(5., 5., 9.),
                  pbc=(True, True, False))
    atoms.center()
    calc = GPAW('scf-criterion.gpw', txt=None)  # checks loading
    atoms.calc = calc
    atoms.get_potential_energy()
    fermilevel = calc.wfs.fermi_level
    workfunctions2 = Ha * calc.hamiltonian.get_workfunctions(fermilevel)

    assert workfunctions1[0] == pytest.approx(workfunctions2[1])
    assert workfunctions1[1] == pytest.approx(workfunctions2[0])
    assert calc.scf.criteria['work function'].tol == pytest.approx(1.0)

    # Try import syntax.
    convergence = {'eigenstates': 1.0,
                   'density': 1.0,
                   'energy': 1.0,
                   'work function': WorkFunction(0.5)}
    calc.set(convergence=convergence)
    atoms.get_potential_energy()
    assert calc.scf.criteria['work function'].tol == pytest.approx(0.5)

    # Switch to H2 for faster calcs.
    for atom in atoms:
        atom.symbol = 'H'

    # Change a default.
    convergence = {'energy': Energy(2.0, n_old=4),
                   'density': np.inf,
                   'eigenstates': np.inf}
    calc.set(convergence=convergence)
    atoms.get_potential_energy()
    assert calc.scf.criteria['energy'].n_old == 4


def test_scf_custom_criterion(in_tmp_dir):
    """Simulate a user creating their own custom convergence criterion,
    saving the .gpw file, and re-loading it. It will warn the user at two
    points."""
    convergence = {'eigenstates': 1.0,
                   'density': 1.0,
                   'energy': 1.0,
                   'custom': FourIterations()}

    atoms = Atoms('HF', [(0., 0.5, 0.5),
                         (0., 0.4, -0.4)],
                  cell=(5., 5., 9.),
                  pbc=(True, True, False))
    atoms.center()
    atoms.rattle()
    calc = GPAW(h=0.3,
                nbands=-1,
                convergence=convergence,
                txt='out.txt',
                poissonsolver={'dipolelayer': 'xy'})
    atoms.calc = calc
    with pytest.warns(UserWarning):
        # Warns the user that their criterion must be a unique instance,
        # and can't be saved/loaded.
        atoms.get_potential_energy()
    calc.write('four.gpw')
    with pytest.warns(UserWarning):
        # Warns the user that their criterion did not load.
        calc = GPAW('four.gpw', txt='out2.txt')
    atoms[1].x += 0.1
    atoms.calc = calc
    atoms.get_potential_energy()
