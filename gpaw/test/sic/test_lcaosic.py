import pytest

from gpaw import GPAW, LCAO
from gpaw.directmin.etdm_lcao import LCAOETDM
from ase import Atoms
import numpy as np


@pytest.mark.later
@pytest.mark.sic
def test_lcaosic(in_tmp_dir):
    """
    Test Perdew-Zunger Self-Interaction
    Correction  in LCAO mode using ETDM
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
    H2O.center(vacuum=4.0)

    calc = GPAW(mode=LCAO(force_complex_dtype=True),
                h=0.22,
                occupations={'name': 'fixed-uniform'},
                eigensolver=LCAOETDM(localizationtype='PM_PZ',
                                     localizationseed=42,
                                     functional={'name': 'PZ-SIC',
                                                 'scaling_factor':
                                                     (0.5, 0.5)}),
                convergence={'eigenstates': 1e-4},
                mixer={'backend': 'no-mixing'},
                nbands='nao',
                symmetry='off'
                )
    H2O.calc = calc
    e = H2O.get_potential_energy()
    f = H2O.get_forces()

    assert e == pytest.approx(-12.16352, abs=1e-3)

    f2 = np.array([[-4.21747862, -4.63118948, 0.00303988],
                   [5.66636141, -0.51037693, -0.00049136],
                   [-1.96478031, 5.4043045, -0.0006107]])
    assert f2 == pytest.approx(f, abs=0.1)

    numeric = False
    if numeric:
        from ase.calculators.test import numeric_force
        f_num = np.array([[numeric_force(H2O, a, i)
                          for i in range(3)]
                         for a in range(len(H2O))])
        print('Numerical forces')
        print(f_num)
        print(f - f_num, np.abs(f - f_num).max())

    calc.write('h2o.gpw', mode='all')
    from gpaw import restart
    H2O, calc = restart('h2o.gpw', txt='-')
    H2O.positions += 1.0e-6
    f3 = H2O.get_forces()
    niter = calc.get_number_of_iterations()
    assert niter == pytest.approx(4, abs=3)
    assert f2 == pytest.approx(f3, abs=0.1)
