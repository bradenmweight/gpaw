import pytest
from ase.build import molecule
from gpaw import GPAW
import gpaw.mom as mom
from gpaw.test import equal


@pytest.mark.mom
def test_mom_lcao():
    f_n = [[1., 1., 1., 1., 0., 0.],
           [1., 1., 1., 1., 0., 0.]]
    dE_ref = [8.1458184518, 6.1579949366]

    for s in [0, 1]:
        # Ground state calculation

        calc = GPAW(mode='lcao',
                    basis='dzp',
                    nbands=6,
                    h=0.2,
                    xc='PBE',
                    spinpol=True,
                    symmetry='off',
                    convergence={'energy': 100,
                                 'density': 1e-3,
                                 'bands': -1})

        H2O = molecule('H2O')
        H2O.center(vacuum=3)
        H2O.calc = calc

        E_gs = H2O.get_potential_energy()

        # Excited-state MOM calculation
        f_n_exc = f_n.copy()
        f_n_exc[0][3] -= 1.
        f_n_exc[s][4] += 1.

        mom.mom_calculation(calc, H2O, f_n)

        E_es = H2O.get_potential_energy()

        dE = E_es - E_gs
        print(dE)
        equal(dE, 8.1458184518, 0.011)
