import pytest

from ase import Atoms
from gpaw import GPAW, LCAO
import numpy as np


def test_directmin_lcao_numerical_hessian(in_tmp_dir):
    """
    Test complex numerical Hessian
    w.r.t rotation parameters in LCAO

    :param in_tmp_dir:
    :return:
    """

    calc = GPAW(xc='PBE',
                mode=LCAO(force_complex_dtype=True),
                h=0.25,
                basis='dz(dzp)',
                spinpol=False,
                eigensolver={'name': 'etdm',
                             'representation': 'u-invar'},
                occupations={'name': 'fixed-uniform'},
                mixer={'backend': 'no-mixing'},
                nbands='nao',
                symmetry='off',
                )

    atoms = Atoms('H', positions=[[0, 0, 0]])
    atoms.center(vacuum=5.0)
    atoms.set_pbc(False)
    atoms.calc = calc
    atoms.get_potential_energy()

    hess_a, hess_n = calc.wfs.eigensolver.finite_diff_appr_of_derivative(
        calc.hamiltonian,
        calc.wfs,
        calc.density,
        what2calc='hessian'
    )
    hess_nt = np.asarray([[1.32720630e+00, -1.93947467e-11],
                         [3.95786680e-09, 1.14599176e+00]])
    assert hess_n == pytest.approx(hess_nt, abs=1e-4)

    a_mat_u = {0: [np.sqrt(2) * np.pi / 4.0 + 1.0j * np.sqrt(2) * np.pi / 4.0]}
    calc.wfs.eigensolver.rotate_wavefunctions(calc.wfs,
                                              a_mat_u,
                                              {0: calc.wfs.bd.nbands},
                                              calc.wfs.eigensolver.c_nm_ref)
    c_nm = {x: calc.wfs.kpt_u[x].C_nM.copy()
            for x in range(len(calc.wfs.kpt_u))}
    hess_a, hess_n = calc.wfs.eigensolver.finite_diff_appr_of_derivative(
        calc.hamiltonian,
        calc.wfs,
        calc.density,
        c_nm_ref=c_nm,
        what2calc='hessian'
    )
    hess_nt = np.asarray([[-1.08209601e+00, -1.11022302e-09],
                         [8.50014503e-10, -9.37664521e-01]])
    assert hess_n == pytest.approx(hess_nt, abs=1e-4)
