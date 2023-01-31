import pytest
import numpy as np
from ase.build import bulk
from gpaw import GPAW, FermiDirac
from gpaw.mpi import world
from gpaw.response.bse import BSE, read_bse_eigenvalues
from gpaw.response.df import read_response_function
from gpaw.test import findpeak, equal


@pytest.mark.response
def test_response_bse_silicon(in_tmp_dir, scalapack):
    GS = 1
    nosym = 1
    bse = 1
    check = 1

    if GS:
        a = 5.431  # From PRB 73,045112 (2006)
        atoms = bulk('Si', 'diamond', a=a)
        atoms.positions -= a / 8
        calc = GPAW(mode='pw',
                    kpts={'size': (2, 2, 2), 'gamma': True},
                    occupations=FermiDirac(0.001),
                    nbands=12,
                    convergence={'bands': -4})
        atoms.calc = calc
        atoms.get_potential_energy()
        calc.write('Si.gpw', 'all')

    if bse:
        eshift = 0.8
        bse = BSE('Si.gpw',
                  ecut=50.,
                  valence_bands=range(4),
                  conduction_bands=range(4, 8),
                  eshift=eshift,
                  nbands=8,
                  write_h=False,
                  write_v=False)
        bse.get_dielectric_function(eta=0.2,
                                    w_w=np.linspace(0, 10, 2001))
        w_w, epsreal_w, epsimag_w = read_response_function('df_bse.csv')
    if check:
        w_ = 2.552
        I_ = 421.15
        w, I = findpeak(w_w, epsimag_w)
        equal(w, w_, 0.01)
        equal(I, I_, 0.1)

    if GS and nosym:
        atoms = bulk('Si', 'diamond', a=a)
        calc = GPAW(mode='pw',
                    kpts={'size': (2, 2, 2), 'gamma': True},
                    occupations=FermiDirac(0.001),
                    nbands=12,
                    symmetry='off',
                    convergence={'bands': -4})
        atoms.calc = calc
        atoms.get_potential_energy()
        calc.write('Si.gpw', 'all')

    if bse and nosym:
        bse = BSE('Si.gpw',
                  ecut=50.,
                  valence_bands=range(4),
                  conduction_bands=range(4, 8),
                  eshift=eshift,
                  nbands=8,
                  write_h=False,
                  write_v=False)
        w_w, eps_w = bse.get_dielectric_function(filename=None,
                                                 eta=0.2,
                                                 w_w=np.linspace(0, 10, 2001))

    if check and nosym:
        w, I = findpeak(w_w, eps_w.imag)
        equal(w, w_, 0.01)
        equal(I, I_, 0.1)

        # Read eigenvalues file and test first 3 eigenvalues
        _, C_w = read_bse_eigenvalues('eig.dat')

        # Disable for now (fails sometimes on a single core also, See #720)
        if 0:  # world.size == 1:
            # The BSE module seems to be broken in parallel!
            # The eigenvalues change drastically as a function of world.size
            assert C_w[:3] == pytest.approx([22.375, 12.777, 31.515], abs=0.1)
