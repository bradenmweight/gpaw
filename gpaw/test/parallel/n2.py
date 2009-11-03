from ase import *
from gpaw import GPAW, sl_diagonalize, extra_parameters, FermiDirac
from gpaw.eigensolvers.rmm_diis import RMM_DIIS
import gpaw.mpi as mpi
sl_diagonalize += [2, 2, 2, 'd']
extra_parameters['blacs'] = True
a = molecule('N2')
a.center(vacuum=3)
c = GPAW(nbands=8,
         parsize_bands=mpi.size,
         eigensolver=RMM_DIIS(0),
         occupations=FermiDirac(width=0.01))
a.set_calculator(c)
e1 = a.get_potential_energy()
c.write('N2')
c.write('N2', mode='all')
wf1 = c.get_pseudo_wave_function(4)
assert c.gd.comm.size == 1

extra_parameters['blacs'] = False
del sl_diagonalize[:]
c = GPAW(occupations=FermiDirac(width=0.01))
a.set_calculator(c)
e2 = a.get_potential_energy()
wf2 = c.get_pseudo_wave_function(4)
assert c.gd.comm.size == 4
assert abs(e1 -e2) < 1e-9
assert (wf1 - wf2).ptp() < 1e-9 or (wf1 + wf2).ptp() < 1e-9
