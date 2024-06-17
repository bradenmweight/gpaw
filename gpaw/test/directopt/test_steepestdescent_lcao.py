import pytest

from gpaw import GPAW
from ase import Atoms


@pytest.mark.do
def test_steepestdescent_lcao(in_tmp_dir):
    """
    Test steepest descent and conjugate gradients
    search direction algorithms
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
                h=0.3,
                spinpol=True,
                convergence={'energy': 0.1,
                             'eigenstates': 1e-4,
                             'density': 1e-4},
                eigensolver={'name': 'etdm-lcao'},
                occupations={'name': 'fixed-uniform'},
                mixer={'backend': 'no-mixing'},
                nbands='nao',
                symmetry='off')
    atoms.calc = calc

    for sd_algo in ['sd', 'fr-cg']:
        calc = calc.new(eigensolver={'name': 'etdm-lcao',
                                     'searchdir_algo': sd_algo})
        atoms.calc = calc
        e = atoms.get_potential_energy()
        assert e == pytest.approx(6.021948, abs=1.0e-5)
