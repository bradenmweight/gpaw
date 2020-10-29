
import numpy as np
from gpaw import GPAW, PW
from ase import Atoms
from gpaw.nlopt.shift import get_shift
from gpaw.nlopt.matrixel import make_nlodata
from gpaw.mpi import world
from gpaw.nlopt.basic import is_file_exist


def test_shift(in_tmp_dir):

    # Check for Hydrogen atom
    atoms = Atoms('H', cell=(3 * np.eye(3)), pbc=True)

    # Do a GS and save it
    if is_file_exist('gs.gpw'):
        calc = GPAW(
            mode=PW(600), symmetry={'point_group': False},
            kpts={'size': (2, 2, 2)}, nbands=5, txt=None)
        atoms.calc = calc
        atoms.get_potential_energy()
        calc.write('gs.gpw', 'all')

    # Get the mml
    if is_file_exist('mml.npz'):
        make_nlodata()

    # Do a shift current caclulation
    get_shift(freqs=np.linspace(0, 5, 100))

    # Check it
    if world.rank == 0:
        shift = np.load('shift.npy')

        # Check for nan's
        assert not np.isnan(shift).any()
        # It should be zero (small) since H2 is centro-symm.
        assert np.all(np.abs(shift[1]) < 1e-8)
