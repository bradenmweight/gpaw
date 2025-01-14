from ase.build import fcc100, add_adsorbate
from gpaw import GPAW

slab = fcc100('Al', (2, 2, 2), a=4.05, vacuum=7.5)
add_adsorbate(slab, 'Na', 4.0)
slab.center(axis=2)

slab.calc = GPAW(mode='fd',
                 txt='zero.txt',
                 xc='PBE',
                 setups={'Na': '1'},
                 kpts=(4, 4, 1))
e1 = slab.get_potential_energy()
slab.calc.write('zero.gpw')

slab.pbc = True
slab.calc = slab.calc.new(txt='periodic.txt')
e2 = slab.get_potential_energy()
slab.calc.write('periodic.gpw')

slab.pbc = (True, True, False)
slab.calc = slab.calc.new(poissonsolver={'dipolelayer': 'xy'},
                          txt='corrected.txt')
e3 = slab.get_potential_energy()
slab.calc.write('corrected.gpw')
