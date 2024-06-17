import pytest

from gpaw import GPAW, LCAO
from gpaw.directmin.etdm_lcao import LCAOETDM
from gpaw.directmin.tools import excite
from gpaw.directmin.derivatives import Davidson
from gpaw.mom import prepare_mom_calculation
from ase import Atoms
import numpy as np


@pytest.mark.later
@pytest.mark.sic
def test_gmf_lcaosic(in_tmp_dir):
    """
    test Perdew-Zunger Self-Interaction
    Correction  in LCAO mode using DirectMin
    :param in_tmp_dir:
    :return:
    """

    # Water molecule:
    d = 0.9575
    t = np.pi / 180 * 104.51
    H2O = Atoms('OH2',
                positions=[(0, 0, 0),
                           (d, 0, 0),
                           (d * np.cos(t), d * np.sin(t), 0)])
    H2O.center(vacuum=3.0)

    calc = GPAW(mode=LCAO(),
                basis='sz(dzp)',
                h=0.24,
                occupations={'name': 'fixed-uniform'},
                eigensolver='etdm-lcao',
                convergence={'eigenstates': 1e-4},
                mixer={'backend': 'no-mixing'},
                nbands='nao',
                spinpol=True,
                symmetry='off'
                )
    H2O.calc = calc
    H2O.get_potential_energy()

    calc.set(eigensolver=LCAOETDM(excited_state=True))
    f_sn = excite(calc, 0, 0, spin=(0, 0))
    prepare_mom_calculation(calc, H2O, f_sn)
    H2O.get_potential_energy()

    dave = Davidson(calc.wfs.eigensolver, None)
    appr_sp_order = dave.estimate_sp_order(calc)
    print(appr_sp_order)

    for kpt in calc.wfs.kpt_u:
        f_sn[kpt.s] = kpt.f_n
    calc.set(eigensolver=LCAOETDM(
        partial_diagonalizer={
            'name': 'Davidson', 'logfile': 'test.txt', 'seed': 42,
            'm': 20, 'eps': 5e-3, 'remember_sp_order': True,
            'sp_order': appr_sp_order},
        linesearch_algo={'name': 'max-step'},
        searchdir_algo={'name': 'LBFGS-P_GMF'},
        localizationtype='PM',
        functional={'name': 'PZ-SIC',
                    'scaling_factor': (0.5, 0.5)},
        need_init_orbs=False),
        occupations={'name': 'mom', 'numbers': f_sn,
                     'use_fixed_occupations': True})

    e = H2O.get_potential_energy()
    assert e == pytest.approx(-2.007241, abs=1.0e-3)

    f = H2O.get_forces()

    f_num = np.array([[-8.01206297e+00, -1.51553367e+01, 3.60670227e-03],
                      [1.42287594e+01, -9.81724693e-01, -5.09333905e-04],
                      [-4.92299436e+00, 1.55306540e+01, 2.12438557e-03]])

    numeric = False
    if numeric:
        from ase.calculators.test import numeric_force
        f_num = np.array([[numeric_force(H2O, a, i)
                          for i in range(3)]
                         for a in range(len(H2O))])
        print('Numerical forces')
        print(f_num)
        print(f - f_num, np.abs(f - f_num).max())

    assert f == pytest.approx(f_num, abs=0.75)
