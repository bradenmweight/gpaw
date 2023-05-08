"""Test for elph/gmatrix"""
import numpy as np
import pytest

from ase.build import bulk
# from ase.phonons import Phonons

from gpaw import GPAW
from gpaw.mpi import world

from gpaw.elph import ElectronPhononMatrix


@pytest.mark.elph
def test_gmatrix(module_tmp_path, supercell_cache):
    atoms = bulk('Li', crystalstructure='bcc', a=3.51, cubic=True)
    supercell_cache
    elph = ElectronPhononMatrix(atoms, 'supercell', 'elph')
    q = [[0, 0, 0], [0.5, 0.5, 0.5]]

    # phonon = Phonons(atoms, name='elph', supercell=(2,1,1),
    #                  center_refcell=True)
    # phonon.read()
    # w_ql = phonon.band_structure(q, modes=False)
    # print(w_ql)
    # [[0.00124087  0.00044091  0.00132004  0.0421112   0.04212737  0.04218485]
    # [0.03031018  0.03031948  0.03041029  0.03041035  0.04326759  0.04327498]]

    calc = GPAW(mode='lcao',
                basis='sz(dzp)',
                kpts={'size': (2, 2, 2), 'gamma': False},
                symmetry='off',
                txt='li_gs_nosym.txt')
    atoms.calc = calc
    atoms.get_potential_energy()

    g_sqklnn = elph.bloch_matrix(calc, k_qc=q,
                                 savetofile=False, prefactor=False)

    assert g_sqklnn.shape == (1, 2, 8, 6, 2, 2)

    # NOTE: It seems g is Hermitian if q=0 and symmetric otherwise. CHECK THIS!

    if world.rank == 0:
        # Li has lots of degeneratephonon modes, average/sum those

        # q = 0
        # accoustic sum rule
        assert np.allclose(g_sqklnn[0, 0, :, 0:3], 0.)

        g_knn = np.sum(g_sqklnn[0, 0], axis=1, initial=3)  # modes 4-6
        assert g_knn.shape == (8, 2, 2)
        # Hermitian
        assert np.allclose(g_knn[:, 0, 1], g_knn[:, 1, 0].conj())
        # and check one specific value
        print(g_knn[0, 0, 1])
        assert g_knn[0, 0, 1] == pytest.approx(3.00000 - 0.162686j, rel=4e-2)
        # print(g_knn)

        # q = 1
        g_knn = np.sum(g_sqklnn[0, 1], axis=1, initial=4)  # modes 5-6
        assert g_knn.shape == (8, 2, 2)
        # Hermitian, actually symmetric
        assert np.allclose(g_knn[:, 0, 1], g_knn[:, 1, 0])
        # and check one specific value
        # print(g_knn)
        print(g_knn[0, 0, 1])
        assert g_knn[0, 0, 1] == pytest.approx(3.50743 - 0.237765j, rel=4e-2)
