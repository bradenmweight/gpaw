"""Test of ase wannier using gpaw."""
import numpy as np
import pytest
from ase.build import molecule
from ase.dft.wannier import Wannier

from gpaw.mpi import world
from gpaw import GPAW

pytestmark = pytest.mark.skipif(world.size > 1,
                                reason='world.size > 1')


@pytest.mark.wannier
def test_ase_features_asewannier(in_tmp_dir):
    calc = GPAW(mode='fd', gpts=(32, 32, 32), nbands=4)
    atoms = molecule('H2', calculator=calc)
    atoms.center(vacuum=3.)
    e = atoms.get_potential_energy()

    pos = atoms.positions + np.array([[0, 0, .2339], [0, 0, -.2339]])
    com = atoms.get_center_of_mass()

    wan = Wannier(nwannier=2, calc=calc, initialwannier='bloch')
    assert wan.get_functional_value() == pytest.approx(2.964, abs=1e-3)
    assert np.linalg.norm(wan.get_centers() - [com, com]) == pytest.approx(
        0, abs=1e-4)

    wan = Wannier(nwannier=2, calc=calc, initialwannier='projectors')
    assert wan.get_functional_value() == pytest.approx(3.100, abs=2e-3)
    assert np.linalg.norm(wan.get_centers() - pos) == pytest.approx(
        0, abs=1e-3)

    wan = Wannier(nwannier=2,
                  calc=calc,
                  initialwannier=[[0, 0, .5], [1, 0, .5]])
    assert wan.get_functional_value() == pytest.approx(3.100, abs=2e-3)
    assert np.linalg.norm(wan.get_centers() - pos) == pytest.approx(
        0, abs=1e-3)

    wan.localize()
    assert wan.get_functional_value() == pytest.approx(3.100, abs=2e-3)
    assert np.linalg.norm(wan.get_centers() - pos) == pytest.approx(
        0, abs=1e-3)
    assert np.linalg.norm(wan.get_radii() - 1.2393) == pytest.approx(
        0, abs=2e-3)
    eig = np.sort(np.linalg.eigvals(wan.get_hamiltonian(k=0).real))
    assert np.linalg.norm(eig - calc.get_eigenvalues()[:2]) == pytest.approx(
        0, abs=1e-4)

    wan.write_cube(0, 'H2.cube')

    energy_tolerance = 0.002
    assert e == pytest.approx(-6.652, abs=energy_tolerance)


@pytest.mark.wannier
def test_wannier_pw(in_tmp_dir, gpw_files, needs_ase_master):
    calc = GPAW(gpw_files['fancy_si_pw_nosym'])
    wan = Wannier(nwannier=4, calc=calc, fixedstates=4,
                  initialwannier='orbitals')
    wan.localize()
    wan.translate_all_to_cell(cell=(0, 0, 0))
    assert wan.get_functional_value() == pytest.approx(9.853, abs=2e-3)

    # test that one of the Wannier functions is localized between the
    # two Si atoms in the unit cell, as it should be
    center = calc.atoms.positions.mean(axis=0)
    min_dist_to_center = np.min(np.linalg.norm(wan.get_centers()
                                               - center, axis=1))
    assert min_dist_to_center < 1e-3
