import pytest

from gpaw import GPAW, PW
from ase import Atoms
import numpy as np
from gpaw.directmin.etdm_fdpw import FDPWETDM
from gpaw.directmin.tools import excite
from gpaw.mom import prepare_mom_calculation


@pytest.mark.later
@pytest.mark.sic
def test_mom_pwsic(in_tmp_dir, needs_ase_master):
    # Water molecule:
    d = 0.9575
    t = np.pi / 180 * 104.51
    H2O = Atoms('OH2',
                positions=[(0, 0, 0),
                           (d, 0, 0),
                           (d * np.cos(t), d * np.sin(t), 0)])
    H2O.center(vacuum=3.0)

    calc = GPAW(mode=PW(300, force_complex_dtype=True),
                spinpol=True,
                symmetry='off',
                eigensolver=FDPWETDM(converge_unocc=True),
                mixer={'backend': 'no-mixing'},
                occupations={'name': 'fixed-uniform'},
                convergence={'eigenstates': 1e-4}
                )
    H2O.calc = calc
    H2O.get_potential_energy()

    calc.set(eigensolver=FDPWETDM(excited_state=True))
    f_sn = excite(calc, 0, 0, (0, 0))
    prepare_mom_calculation(calc, H2O, f_sn)
    H2O.get_potential_energy()

    calc.set(eigensolver=FDPWETDM(
        excited_state=True,
        need_init_orbs=False,
        functional={'name': 'PZ-SIC',
                    'scaling_factor': (0.5, 0.5)},  # SIC/2
        localizationseed=42,
        localizationtype='PM',
        grad_tol_pz_localization=1.0e-2,
        printinnerloop=False,
        grad_tol_inner_loop=1.0e-2),
        convergence={'eigenstates': 1e-3,
                     'density': 1e-3})

    e = H2O.get_potential_energy()
    assert e == pytest.approx(-3.302431, abs=0.2)

    f = H2O.get_forces()

    # Numeric forces, generated by disabled code below
    f_num = np.array([[-2.85022, -3.66863, -0.009059],
                      [3.916824, -0.204146, -0.000065],
                      [-1.204687, 3.822242, 0.000353]])
    assert f == pytest.approx(f_num, abs=0.3)

    numeric = False
    if numeric:
        from ase.calculators.test import numeric_force
        f_num = np.array([[numeric_force(H2O, a, i)
                           for i in range(3)]
                          for a in range(len(H2O))])
        print('Numerical forces')
        print(f_num)
        print(f - f_num, np.abs(f - f_num).max())
