import pytest
from ase import Atoms
from ase.optimize import BFGS
from gpaw.new.ase_interface import GPAW
# from gpaw import GPAW


@pytest.mark.serial
def test_tb_mode_molecule():
    atoms = Atoms('LiH',
                  [[0, 0.1, 0.2],
                   [0, 0, 1.4]])
    atoms.calc = GPAW(
        mode='tb',
        symmetry='off',  # ??? should not be needed
        txt='-')  # None)
    atoms.get_potential_energy()
    # f1 = atoms.get_forces()
    # f2 = numeric_forces(atoms)
    # assert abs(f1 - f2).max() < 0.0005


@pytest.mark.serial
def test_tb_mode_bulk():
    a = 2.0
    atoms = Atoms('Li',
                  cell=[a, a, a],
                  pbc=True)
    atoms.calc = GPAW(
        mode='tb',
        kpts=(2, 2, 2))
    atoms.get_potential_energy()
    return
    from ase.filters import ExpCellFilter
    from ase.calculators.test import numeric_stress
    f = atoms.get_forces()
    assert abs(f).max() < 0.0001
    e = atoms.get_potential_energy()
    s = atoms.get_stress()
    print(a, e, s)
    s2 = numeric_stress(atoms)
    print(s2)
    assert abs(s - s2).max() < 0.0001
    BFGS(ExpCellFilter(atoms)).run(fmax=0.002)
    s = atoms.get_stress()
    print(s)
    assert abs(s).max() < 0.0001
    assert atoms.cell[0, 0] == pytest.approx(2.044, abs=0.0005)
