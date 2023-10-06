import pytest

from gpaw import GPAW, LCAO, restart
from ase import Atoms
import numpy as np
from gpaw.directmin.etdm_lcao import LCAOETDM
from gpaw.directmin.tools import excite
from gpaw.mom import prepare_mom_calculation


@pytest.mark.sic
def test_mom_lcaosic(in_tmp_dir):
    # Water molecule:
    d = 0.9575
    t = np.pi / 180 * 104.51
    H2O = Atoms('OH2',
                positions=[(0, 0, 0),
                           (d, 0, 0),
                           (d * np.cos(t), d * np.sin(t), 0)])
    H2O.center(vacuum=3.0)

    calc = GPAW(mode=LCAO(force_complex_dtype=True),
                h=0.24,
                basis='sz(dzp)',
                spinpol=True,
                symmetry='off',
                eigensolver='etdm-lcao',
                mixer={'backend': 'no-mixing'},
                occupations={'name': 'fixed-uniform'},
                convergence={'eigenstates': 1e-4}
                )
    H2O.calc = calc
    H2O.get_potential_energy()

    calc.set(eigensolver=LCAOETDM(excited_state=True))
    f_sn = excite(calc, 0, 0, (0, 0))
    prepare_mom_calculation(calc, H2O, f_sn)
    H2O.get_potential_energy()

    test_restart = False
    if test_restart:
        calc.write('h2o.gpw', mode='all')
        H2O, calc = restart('h2o.gpw', txt='-')

    calc.set(eigensolver=LCAOETDM(searchdir_algo={'name': 'l-sr1p'},
                                  linesearch_algo={'name': 'max-step'},
                                  need_init_orbs=False,
                                  localizationtype='PM_PZ',
                                  localizationseed=42,
                                  functional={'name': 'pz-sic',
                                              'scaling_factor': (0.5, 0.5)}),
             convergence={'eigenstates': 1e-2})

    e = H2O.get_potential_energy()
    assert e == pytest.approx(-2.007092, abs=5.0e-3)

    f = H2O.get_forces()

    # Saved analytical forces
    f_old = np.array([[-8.76835930e+00, -1.44190257e+01, 2.53281902e-03],
                      [1.46304433e+01, -1.00089622e+00, -1.03430589e-02],
                      [-4.83182230e+00, 1.51551677e+01, -1.20513818e-02]])

    # Numeric forces, generated by disabled code below
    f_num = np.array([[-8.07002033e+00, -1.51005084e+01, -1.67903017e-03],
                      [1.42573640e+01, -9.80976464e-01, 8.86319871e-05],
                      [-4.93107360e+00, 1.54986283e+01, -3.45104541e-03]])

    numeric = False
    if numeric:
        from ase.calculators.test import numeric_force
        f_num = np.array([[numeric_force(H2O, a, i)
                           for i in range(3)]
                          for a in range(len(H2O))])
        print('Numerical forces')
        print(f_num)
        print(f - f_num, np.abs(f - f_num).max())

    assert f == pytest.approx(f_old, abs=0.5)
    assert f == pytest.approx(f_num, abs=1.0)
