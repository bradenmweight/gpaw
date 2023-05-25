import pytest

from gpaw import GPAW, PW
from gpaw.directmin.fdpw.directmin import DirectMin
from gpaw.mom import prepare_mom_calculation
from gpaw.directmin.tools import excite_and_sort
from ase import Atoms
import numpy as np


def test_xc_pw(gpw_files):
    atoms = Atoms('O2', [[0.0, 0.0, 0.0], [0.0, 0.0, 1.4]])
    atoms.center(vacuum=3.0)
    pos = atoms.get_positions()
    pos[0][2] += 0.01
    atoms.set_positions(pos)
    calc = GPAW(mode=PW(200), spinpol=True,
                symmetry='off',
                eigensolver=DirectMin(searchdir_algo={'name': 'LBFGS_P', 'memory': 1}),
                mixer={'name': 'dummy'},
                occupations={'name': 'fixed-uniform'}
                )
    atoms.calc = calc
    atoms.get_potential_energy()
    i, a = 0, 1
    excite_and_sort(calc.wfs, i, a, (0, 0), 'fdpw')
    calc.set(eigensolver=DirectMin(exstopt=True,
                                   searchdir_algo={'name': 'LBFGS_P', 'memory': 1}))
    f_sn = []
    for spin in range(calc.get_number_of_spins()):
        f_n = calc.get_occupation_numbers(spin=spin)
        f_sn.append(f_n)
    prepare_mom_calculation(calc, atoms, f_sn)
    e = atoms.get_potential_energy()
    assert e == pytest.approx(42.362787, abs=1.0e-4)
    f = atoms.get_forces()
    assert np.min(f) == pytest.approx(-10.60363, abs=1.0e-3)
    assert np.max(f) == pytest.approx(10.59878, abs=1.0e-3)

    pos = atoms.get_positions()
    pos[0][2] -= 0.01
    atoms.set_positions(pos)
    e2 = atoms.get_potential_energy()
    assert e2 == pytest.approx(42.257272, abs=1.0e-4)

test_xc_pw(0)