"""First example for using RSF."""
from ase import Atoms
from gpaw.eigensolvers import RMMDIIS
from gpaw.occupations import FermiDirac

h = 0.30
co = Atoms('CO', positions=[(0, 0, 0), (0, 0, 1.15)])
co.center(5)
c = {'energy': 0.1, 'eigenstates': 3, 'density': 3}
calc = GPAW(txt='CO.txt', xc='LCY_PBE', convergence=c,
            eigensolver=RMMDIIS(), h=h,
            occupations=FermiDirac(width=0.0), spinpol=False)
co.set_calculator(calc)
co.get_potential_energy()
