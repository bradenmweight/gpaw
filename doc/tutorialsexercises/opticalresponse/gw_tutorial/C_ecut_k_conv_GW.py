from ase.build import bulk

from gpaw import GPAW, FermiDirac
from gpaw import PW
from gpaw.response.g0w0 import G0W0

a = 3.567
atoms = bulk('C', 'diamond', a=a)

for k in [6, 8, 10, 12]:
    calc = GPAW(mode=PW(600),
                parallel={'domain': 1},
                kpts={'size': (k, k, k), 'gamma': True},
                xc='LDA',
                basis='dzp',
                occupations=FermiDirac(0.001),
                txt=f'C_groundstate_{k}.txt')

    atoms.calc = calc
    atoms.get_potential_energy()

    calc.diagonalize_full_hamiltonian()
    calc.write(f'C_groundstate_{k}.gpw', mode='all')

    for i, ecut in enumerate([100, 200, 300, 400]):
        gw = G0W0(calc=f'C_groundstate_{k}.gpw',
                  bands=(3, 5),
                  ecut=ecut,
                  kpts=[0],
                  integrate_gamma='WS',
                  filename=f'C-g0w0_k{k}_ecut{ecut}')

        result = gw.calculate()
