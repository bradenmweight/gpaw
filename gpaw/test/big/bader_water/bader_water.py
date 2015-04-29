from ase.structure import molecule
from ase.io import write
from ase.units import Bohr
from gpaw import GPAW

atoms = molecule('H2O')
atoms.center(vacuum=3.5)
atoms.calc = GPAW(h=0.17)#mode='lcao')
atoms.get_potential_energy()
rho = atoms.calc.get_all_electron_density(gridrefinement=4)
write('ae4h.cube', atoms, data=rho * Bohr**3)
#rho = atoms.calc.get_pseudo_density(gridrefinement=2)
#write('ps2.cube', atoms, data=rho * Bohr**3)
