from gpaw import GPAW, LCAO, FermiDirac
from ase import Atoms
import numpy as np
from gpaw.directmin.directmin_lcao import DirectMinLCAO

# Water molecule:
d = 0.9575
t = np.pi / 180 * 104.51
H2O = Atoms('OH2',
            positions=[(0, 0, 0),
                       (d, 0, 0),
                       (d * np.cos(t), d * np.sin(t), 0)])
H2O.center(vacuum=5.0)

calc = GPAW(mode=LCAO(force_complex_dtype=True),
            basis='dzp',
            occupations=FermiDirac(width=0.0, fixmagmom=True),
            eigensolver=DirectMinLCAO(
                odd_parameters={'name': 'PZ_SIC',  # half-SIC
                                'scaling_factor': (0.5, 0.5)}),
            mixer={'method': 'dummy'},
            nbands='nao'
            )
H2O.set_calculator(calc)
H2O.get_potential_energy()
H2O.get_forces()
