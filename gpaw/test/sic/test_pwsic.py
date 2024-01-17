import pytest

from gpaw import GPAW, PW, restart
from ase import Atoms
import numpy as np
from gpaw.directmin.etdm_fdpw import FDPWETDM
from ase.dft.bandgap import bandgap
from ase.units import Ha


@pytest.mark.later
@pytest.mark.sic
def test_pwsic(in_tmp_dir, needs_ase_master):
    """
    test Perdew-Zunger Self-Interaction
    Correction  in PW mode using DirectMin
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

    calc = GPAW(mode=PW(300, force_complex_dtype=True),
                occupations={'name': 'fixed-uniform'},
                eigensolver=FDPWETDM(
                    functional={'name': 'pz-sic',
                                'scaling_factor': (0.5, 0.5)},
                    localizationseed=42,
                    localizationtype='FB_ER',
                    grad_tol_pz_localization=5.0e-3,
                    maxiter_pz_localization=200,
                    converge_unocc=True),
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

    assert e == pytest.approx(-9.98919, abs=1e-3)
    # Numeric forces, generated by disabled code below
    f2 = np.array([[0.22161312, -0.98564396, -0.00204214],
                   [-0.34986867, 0.17494903, 0.00029861],
                   [0.01085528, 0.56112341, 0.00129632]])
    assert f2 == pytest.approx(f, abs=3e-2)
    assert gap == pytest.approx(9.555, abs=1e-2)

    numeric = False
    if numeric:
        from ase.calculators.test import numeric_force
        f_num = np.array([[numeric_force(H2O, a, i)
                           for i in range(3)]
                          for a in range(len(H2O))])
        print('Numerical forces')
        print(f_num)
        print(f - f_num, np.abs(f - f_num).max())

    #
    calc.write('h2o.gpw', mode='all')

    H2O, calc = restart('h2o.gpw', txt='-')
    H2O.positions += 1.0e-6
    f3 = H2O.get_forces()
    niter = calc.get_number_of_iterations()
    assert niter == pytest.approx(4, abs=3)
    assert f2 == pytest.approx(f3, abs=3e-2)
