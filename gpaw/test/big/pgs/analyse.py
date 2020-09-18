from pathlib import Path
from ase.io import read
from gpaw import GPAW
from gpaw.point_groups import SymmetryChecker

charges = {'Ih': -2,
           'Ico': -2,
           'Th': 2}

for path in Path().glob('Oh-*.xyz'):
    print(path, end='')
    pg = path.name.split('-')[0]
    atoms = read(path)
    atoms.center(vacuum=5)
    if 1:
        atoms.calc = GPAW(h=0.2,#mode='lcao',
                          charge=charges.get(pg, 0.0),
                          txt=path.with_suffix('.txt2'))
        atoms.get_potential_energy()
        atoms.calc.write(path.with_suffix('.gpw'), mode='all')
    elif 0:
        atoms.calc = GPAW(path.with_suffix('.gpw'))
    else:
        pass
    center = atoms.get_center_of_mass()
    R = atoms.positions
    if pg in {'Ico', 'Ih'}:
        z = R[17] - R[23]
        x = (R[5] + R[11]) / 2 - R[9]
    else:
        x = None
        z = None
    checker = SymmetryChecker(pg, center, 4.5,
                              x=x, z=z)
    c = checker.check_atoms(atoms)
    print(c)

    checker.check_calculation(atoms.calc,
                              0, atoms.calc.get_number_of_bands(),
                              output=path.with_suffix('.sym'))
    for n in range(atoms.calc.get_number_of_bands()):
        result = checker.check_band(atoms.calc, n)
        characters = result['characters']
        best = result['symmetry']
        for sym, value in characters.items():
            if sym != best:
                if abs(value) > 0.1 * characters[best]:
                    print(n, characters)
                # assert abs(value) < 0.1 * characters[best]

"""
Oh-F6S.xyz Oh True
"""
