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
    mom_after_canonical = False  # Test MOM after canonical only
    if mom_after_canonical:
        momevery = np.inf
    else:
        momevery = 3
    calc.set(eigensolver=FDPWETDM(exstopt=True,
                                  momevery=momevery,
                                  restart_canonical=False,
                                  printinnerloop=False))
    f_sn = excite(calc, 0, 0, (0, 0))
    prepare_mom_calculation(calc, atoms, f_sn)

    def rotate_homo_lumo(calc=calc):
        angle = 70 * np.pi / 180.0
        iters = calc.get_number_of_iterations()
        if iters == 3:
            psit_nG_old = calc.wfs.kpt_u[0].psit_nG.copy()
            calc.wfs.kpt_u[0].psit_nG[3] = \
                np.cos(angle) * psit_nG_old[3] + np.sin(angle) * psit_nG_old[4]
            calc.wfs.kpt_u[0].psit_nG[4] = \
                np.cos(angle) * psit_nG_old[4] - np.sin(angle) * psit_nG_old[3]
            for kpt in calc.wfs.kpt_u:
                calc.wfs.pt.integrate(kpt.psit_nG, kpt.P_ani, kpt.q)

    calc.attach(rotate_homo_lumo, 1)
    e = atoms.get_potential_energy()
    assert e == pytest.approx(0.027152, abs=1.0e-3)

    f = atoms.get_forces()

    # Numeric forces, generated by disabled code below
    f2 = np.array([[-4.07053122, -5.46395110, -2.54491307e-04],
                   [5.57197214, -1.00373986e-01, 1.87598849e-04],
                   [-1.52873845, 5.38474614, 2.07965790e-04]])
    assert f2 == pytest.approx(f, abs=0.1)

    numeric = False
    if numeric:
        calc.observers = []
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
    atoms.calc.results.pop('energy')
    atoms.calc.scf.converged = False
    e2 = atoms.get_potential_energy()
    niter = calc.get_number_of_iterations()
    assert niter == pytest.approx(4, abs=3)
    assert e == pytest.approx(e2, abs=1.0e-3)

    prepare_mom_calculation(calc, atoms, f_sn, use_fixed_occupations='True')
    e2 = atoms.get_potential_energy()
    for spin in range(calc.get_number_of_spins()):
        f_n = calc.get_occupation_numbers(spin=spin)
        assert (np.allclose(f_sn[spin], f_n))
        assert (np.allclose(f_sn[spin], calc.wfs.occupations.numbers[spin]))
    niter = calc.get_number_of_iterations()
    assert niter == pytest.approx(4, abs=3)
    assert e == pytest.approx(e2, abs=1.0e-3)
