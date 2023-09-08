import pytest
from ase import Atoms

from gpaw.calculator import GPAW as OldGPAW
from gpaw.new.ase_interface import GPAW as NewGPAW


@pytest.mark.parametrize('random', [False, True])
def test_refactored_code(in_tmp_dir, random):
    new('n', random)


def new(x: str, random: bool) -> None:
    params = {'mode': {'name': 'fd'},
              'random': random,
              'kpts': (4, 1, 1),
              'spinpol': not True,
              'convergence': {'maximum iterations': 200}}

    if x == 'n':
        GPAW = NewGPAW
    else:
        GPAW = OldGPAW

    atoms = Atoms('H2', cell=[2, 2, 3], pbc=True)
    atoms.positions[1, 2] = 0.8
    atoms.calc = GPAW(**params,
                      txt=f'{x}s.txt')
    f = atoms.get_forces()
    f2 = atoms.calc.get_property('forces', allow_calculation=False)
    assert (f == f2).all()
    e = atoms.get_potential_energy()

    atoms.get_dipole_moment()
    print(e)
    print(f)
    if 1:
        atoms.positions[1, 2] = 0.75
        f2 = atoms.get_forces()
        e2 = atoms.get_potential_energy()
        print(f2)
        print(e2)


if __name__ == '__main__':
    import sys
    for x in sys.argv[1:]:
        new(x, False)
