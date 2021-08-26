import pytest

from gpaw import GPAW, LCAO
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
        g_a, g_n = calc.wfs.eigensolver.finite_diff_appr_of_derivative(
            ham, wfs, dens, random_amat=True, update_c_nm_ref=True)
        for x, y in zip(g_a[0], g_n[0]):
            assert x.real == pytest.approx(y.real, abs=1.0e-2)
            assert x.imag == pytest.approx(y.imag, abs=1.0e-2)
