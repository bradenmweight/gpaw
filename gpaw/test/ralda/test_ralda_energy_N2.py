from gpaw.mpi import world
from ase import Atoms
from ase.build import molecule
from gpaw import GPAW, PW
from gpaw.xc.fxc import FXCCorrelation
from gpaw.test import equal


def test_ralda_ralda_energy_N2(in_tmp_dir, scalapack):
    if world.size == 1:
        scalapack1 = None
        scalapack2 = None
    elif world.size == 2:
        scalapack1 = (2, world.size // 2, 32)
        scalapack2 = None
    else:
        scalapack1 = (2, world.size // 2, 32)
        scalapack2 = (2, world.size // 2, 32)

    # N2
    N2 = molecule('N2')
    N2.set_cell((2.5, 2.5, 3.5))
    N2.center()
    calc = GPAW(mode=PW(force_complex_dtype=True),
                eigensolver='rmm-diis',
                xc='LDA',
                nbands=16,
                basis='dzp',
                parallel={'domain': 1},
                convergence={'density': 1.e-6})
    N2.calc = calc
    N2.get_potential_energy()
    calc.diagonalize_full_hamiltonian(nbands=80, scalapack=scalapack1)
    calc.write('N2.gpw', mode='all')

    ralda = FXCCorrelation('N2.gpw', xc='rALDA')
    Ec_N2 = ralda.calculate(ecut=[50])

    # N
    N = Atoms('N', [(0, 0, 0)])
    N.set_cell((2.5, 2.5, 3.5))
    N.center()
    calc = GPAW(mode=PW(force_complex_dtype=True),
                eigensolver='rmm-diis',
                xc='LDA',
                basis='dzp',
                nbands=8,
                hund=True,
                parallel={'domain': 1},
                convergence={'density': 1.e-6})
    N.calc = calc
    N.get_potential_energy()
    calc.diagonalize_full_hamiltonian(nbands=80, scalapack=scalapack2)
    calc.write('N.gpw', mode='all')

    ralda = FXCCorrelation('N.gpw', xc='rALDA')
    Ec_N = ralda.calculate(ecut=[50])

    equal(Ec_N2, -6.17, 0.01,)
    equal(Ec_N, -1.1085, 0.001)
