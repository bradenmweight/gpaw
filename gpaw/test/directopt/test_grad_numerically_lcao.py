import pytest

from gpaw import GPAW, LCAO
from gpaw.directmin.etdm import random_a
from ase import Atoms


def test_gradient_numerically_lcao(in_tmp_dir):
    """
    test exponential transformation
    direct minimization method for KS-DFT in LCAO
    :param in_tmp_dir:
    :return:
    """

    # Water molecule:
    atoms = Atoms('H3', positions=[(0, 0, 0),
                                   (0.59, 0, 0),
                                   (1.1, 0, 0)])
    atoms.center(vacuum=2.0)
    atoms.set_pbc(True)
    calc = GPAW(mode=LCAO(force_complex_dtype=True),
                basis='sz(dzp)',
                h=0.3,
                spinpol=False,
                convergence={'eigenstates': 10.0,
                             'density': 10.0,
                             'energy': 10.0},
                occupations={'name': 'fixed-uniform'},
                eigensolver={'name': 'etdm',
                             'matrix_exp': 'egdecomp'},
                mixer={'backend': 'no-mixing'},
                nbands='nao',
                symmetry='off'
                )
    atoms.calc = calc
   
    params = [{'name': 'etdm',
               'representation': 'sparse',
              'matrix_exp': 'egdecomp'},
              {'name': 'etdm',
               'representation': 'sparse',
               'matrix_exp': 'pade-approx'},
              {'name': 'etdm',
               'representation': 'u-invar',
               'matrix_exp': 'egdecomp'}]

    for eigsolver in params:
        calc.set(eigensolver=eigsolver)
        atoms.get_potential_energy()
        ham = calc.hamiltonian
        wfs = calc.wfs
        dens = calc.density
        a = random_a(wfs.eigensolver.a_mat_u[0].shape, wfs.dtype)
        wfs.gd.comm.broadcast(a, 0)
        amatu = {0: a}
        g_n = calc.wfs.eigensolver.get_numerical_derivatives(
            ham, wfs, dens, amatu=amatu, update_c_nm_ref=True)
        g_a = calc.wfs.eigensolver.get_analytical_derivatives(
            ham, wfs, dens, amatu=None, update_c_nm_ref=False)
        for x, y in zip(g_a[0], g_n[0]):
            assert x.real == pytest.approx(y.real, abs=1.0e-2)
            assert x.imag == pytest.approx(y.imag, abs=1.0e-2)
