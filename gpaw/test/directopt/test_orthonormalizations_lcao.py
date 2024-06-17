import pytest

import numpy as np

from gpaw import GPAW
from ase import Atoms


@pytest.mark.do
def test_orthonormalizations_lcao(in_tmp_dir):
    """
    Test Loewdin and Gram-Schmidt orthonormalization
    of orbitals in LCAO
    :param in_tmp_dir:
    :return:
    """

    atoms = Atoms('H3', positions=[(0, 0, 0),
                                   (0.59, 0, 0),
                                   (1.1, 0, 0)])
    atoms.set_initial_magnetic_moments([1, 0, 0])

    atoms.center(vacuum=2.0)
    atoms.set_pbc(False)
    calc = GPAW(mode='lcao',
                basis='sz(dzp)',
                h=0.3,
                spinpol=True,
                convergence={'energy': np.inf,
                             'eigenstates': np.inf,
                             'density': np.inf,
                             'minimum iterations': 1},
                eigensolver={'name': 'etdm-lcao'},
                occupations={'name': 'fixed-uniform'},
                mixer={'backend': 'no-mixing'},
                nbands='nao',
                symmetry='off',
                txt=None)
    atoms.calc = calc
    atoms.get_potential_energy()

    for type in ['loewdin', 'gramschmidt']:
        atoms.positions[0] += 0.1
        calc.initialize_positions(atoms)
        for kpt in calc.wfs.kpt_u:
            calc.wfs.orthonormalize(kpt, type=type)
            overlaps = np.dot(kpt.C_nM.conj(),
                              np.dot(kpt.S_MM, kpt.C_nM.T))
            assert overlaps == pytest.approx(np.identity(3), abs=1e-10)
