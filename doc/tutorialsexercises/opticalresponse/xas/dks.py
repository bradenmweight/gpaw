from math import pi, cos, sin
from ase import Atoms
from ase.parallel import paropen
from gpaw import GPAW, setup_paths, FermiDirac
setup_paths.insert(0, '.')

a = 12.0  # use a large cell

d = 0.9575
t = pi / 180 * 104.51
atoms = Atoms('OH2',
              [(0, 0, 0),
               (d, 0, 0),
               (d * cos(t), d * sin(t), 0)],
              cell=(a, a, a))
atoms.center()

calc1 = GPAW(mode='fd',
             h=0.2,
             txt='h2o_gs.txt',
             xc='PBE')
atoms.calc = calc1
e1 = atoms.get_potential_energy() + calc1.get_reference_energy()

calc2 = GPAW(mode='fd',
             h=0.2,
             txt='h2o_exc.txt',
             xc='PBE',
             charge=-1,
             spinpol=True,
             occupations=FermiDirac(0.0, fixmagmom=True),
             setups={0: 'fch1s'})
atoms[0].magmom = 1
atoms.calc = calc2
e2 = atoms.get_potential_energy() + calc2.get_reference_energy()

with paropen('dks.result', 'w') as fd:
    print('Energy difference:', e2 - e1, file=fd)
