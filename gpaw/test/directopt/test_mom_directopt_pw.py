import pytest

from gpaw import GPAW, PW, restart
from gpaw.directmin.etdm_fdpw import FDPWETDM
from gpaw.mom import prepare_mom_calculation
from gpaw.directmin.tools import excite
from ase import Atoms
import numpy as np


@pytest.mark.do
def test_mom_directopt_pw(in_tmp_dir):
    # Water molecule:
    d = 0.9575
    t = np.pi / 180 * 104.51
    atoms = Atoms('OH2',
                  positions=[(0, 0, 0),
                             (d, 0, 0),
                             (d * np.cos(t), d * np.sin(t), 0)])
    atoms.center(vacuum=4.0)

    calc = GPAW(mode=PW(300),
                spinpol=True,
                symmetry='off',
                eigensolver=FDPWETDM(converge_unocc=True),
                mixer={'backend': 'no-mixing'},
                occupations={'name': 'fixed-uniform'},
                convergence={'eigenstates': 1e-4}
                )
    atoms.calc = calc
    atoms.get_potential_energy()
    calc.write('h2o.gpw', mode='all')

    # Triplet excited state calculation
    calc.set(eigensolver=FDPWETDM(exstopt=True,
                                  need_init_orbs=False))
    f_sn = excite(calc, 0, 1, (0, 1))
    prepare_mom_calculation(calc, atoms, f_sn)

    e = atoms.get_potential_energy()
    assert e == pytest.approx(1.869659, abs=1.0e-3)

    # Mixed-spin excited state calculation
    atoms, calc = restart('h2o.gpw', txt='-')
    # Don't need to set need_init_orbs=False when restarting
    # from file
    calc.set(eigensolver=FDPWETDM(exstopt=True,
                                  printinnerloop=True))
    f_sn = excite(calc, 0, 0, (0, 0))
    prepare_mom_calculation(calc, atoms, f_sn)

    e = atoms.get_potential_energy()
    assert e == pytest.approx(0.027152, abs=1.0e-3)

    f = atoms.get_forces()

    # Numeric forces, generated by disabled code below
    f2 = np.array([[-4.070454, -5.464042, -0.000266],
                   [5.571928, -0.100377, 0.000202],
                   [-1.528699, 5.384741, 0.000204]])
    assert f2 == pytest.approx(f, abs=3e-2)

    numeric = False
    if numeric:
        from ase.calculators.test import numeric_force
        f_num = np.array([[numeric_force(atoms, a, i)
                          for i in range(3)]
                         for a in range(len(atoms))])
        print('Numerical forces')
        print(f_num)
        print(f - f_num, np.abs(f - f_num).max())

    calc.write('h2o.gpw', mode='all')

    # Test restart and fixed occupations
    atoms, calc = restart('h2o.gpw', txt='-')
    for kpt in calc.wfs.kpt_u:
        f_sn[kpt.s] = kpt.f_n
    prepare_mom_calculation(calc, atoms, f_sn, use_fixed_occupations='True')
    e2 = atoms.get_potential_energy()
    for spin in range(calc.get_number_of_spins()):
        f_n = calc.get_occupation_numbers(spin=spin)
        assert (np.allclose(f_sn[spin], f_n))
        assert (np.allclose(f_sn[spin], calc.wfs.occupations.numbers[spin]))
    niter = calc.get_number_of_iterations()
    assert niter == pytest.approx(4, abs=3)
    assert e == pytest.approx(e2, abs=1.0e-3)
