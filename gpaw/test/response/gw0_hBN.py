from gpaw.wavefunctions.pw import PW
from gpaw import GPAW, FermiDirac
import numpy as np
from ase.lattice import bulk
from gpaw.response.g0w0 import G0W0
from gpaw.test import equal

atoms = bulk('BN', 'zincblende', a=3.615)

calc = GPAW(mode = PW(400),
            kpts={'size': (2, 2, 2), 'gamma': True},
            dtype=complex,
            xc='LDA',
            eigensolver='rmm-diis',
            occupations=FermiDirac(0.001))

atoms.set_calculator(calc)
e0 = atoms.get_potential_energy()

calc.diagonalize_full_hamiltonian(scalapack=True)
calc.write('BN_bulk_k2_ecut400_allbands.gpw', mode='all')

gw = G0W0('BN_bulk_k2_ecut400_allbands.gpw',
          bands=(3,5),
          nbands=10,
          nblocks=1,
          method='GW0',
          maxiter=5,
          ecut=40)

result = gw.calculate()

gaps = [3.256, 4.599, 4.800, 4.812, 4.810, 4.808] 

for i in range(result['iqp'].shape[0]):
    equal(np.min(result['iqp'][i,0,:,1])-np.max(result['iqp'][i,0,:,0]), gaps[i], 0.01)
