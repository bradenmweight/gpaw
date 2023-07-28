import pytest

from gpaw import GPAW, FD
from ase import Atoms
import numpy as np
from gpaw.directmin.etdm_fdpw import FDPWETDM
from ase.dft.bandgap import bandgap
from ase.units import Ha


@pytest.mark.sic
def test_fdsic(in_tmp_dir):
    """
    Test Perdew-Zunger Self-Interaction
    Correction in PW mode using DirectMin
    :param in_tmp_dir:
    :return:
    """

    # Water molecule:
    d = 0.9575
    t = np.pi / 180 * (104.51 + 2.0)
    eps = 0.02
    H2O = Atoms('OH2',
                positions=[(0, 0, 0),
                           (d + eps, 0, 0),
                           (d * np.cos(t), d * np.sin(t), 0)])
    H2O.center(vacuum=4.0)

    calc = GPAW(mode=FD(force_complex_dtype=True),
                h=0.25,
                occupations={'name': 'fixed-uniform'},
                eigensolver=FDPWETDM(
                    functional={'name': 'PZ-SIC',
                                'scaling_factor': (0.5, 0.5)},
                    localizationseed=42,
                    localizationtype='FB_ER',
                    grad_tol_pz_localization=5.0e-3,
                    maxiter_pz_localization=200),
                convergence={'eigenstates': 1e-4},
                mixer={'backend': 'no-mixing'},
                symmetry='off',
                spinpol=True
                )
    H2O.calc = calc
    e = H2O.get_potential_energy()
    f = H2O.get_forces()
    efermi = calc.wfs.fermi_levels[0] * Ha
    gap = bandgap(calc, efermi=efermi)[0]
    e_old = -18.144428

    assert e == pytest.approx(e_old, abs=1e-3)
    # Numeric forces, generated by disabled code below
    f_num = np.array([[2.11270273, 4.89616573e-01, -6.00482442e-04],
                      [-2.15829241, 3.54950512e-01, 1.04418211e-04],
                      [6.67703026e-01, -8.89596180e-01, 8.83126024e-05]])

    assert f == pytest.approx(f_num, abs=5e-2)
    assert gap == pytest.approx(10.215, abs=1e-2)

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
    e = H2O.get_potential_energy()
    f = H2O.get_forces()
    niter = calc.get_number_of_iterations()
    assert niter == pytest.approx(4, abs=3)
    assert e == pytest.approx(e_old, abs=1e-3)
    assert f == pytest.approx(f_num, abs=5e-2)
