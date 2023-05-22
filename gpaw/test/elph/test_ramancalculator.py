import numpy as np
import pytest

from ase.utils.filecache import MultiFileJSONCache
from gpaw.calculator import GPAW
from gpaw.lcao.dipoletransition import get_momentum_transitions
from gpaw.elph import ResonantRamanCalculator
from gpaw.mpi import world


def get_random_g(nk, nb):
    g_sqklnn = np.zeros((1, 1, 4, 3, 4, 4), dtype=complex)
    rng = np.random.default_rng()
    tmp = rng.random((4, 4)) + 1j * rng.random((4, 4))
    # make hermitian
    for i in range(4):
        for j in range(i + 1, 4):
            tmp[i, j] = tmp[j, i].conj()
    g_sqklnn[0, 0, 0, 2] = tmp
    return g_sqklnn


@pytest.mark.serial
def test_ramancalculator(gpw_files, tmp_path_factory):
    """Test of ResonantRamanCalculator object"""
    calc = GPAW(gpw_files['bcc_li_lcao_wfs'])
    atoms = calc.atoms
    # Initialize calculator if necessary
    if not hasattr(calc.wfs, 'C_nM'):
        calc.initialize_positions(atoms)
    # need to fiddle with some occupation numnbers as this exampe is
    # not properly converged
    for kpt in calc.wfs.kpt_u:
        kpt.f_n[0] = kpt.weight

    # prepare some required data
    wph_w = np.array([0., 0., 0.1])
    get_momentum_transitions(calc.wfs)
    if world.rank == 0:
        g_sqklnn = get_random_g(4, 4)
        np.save("gsqklnn.npy", g_sqklnn)

    rrc = ResonantRamanCalculator(calc, wph_w)

    # check reading of file cache
    check_cache = MultiFileJSONCache("Rlab")
    assert check_cache["phonon_frequencies"] == pytest.approx(wph_w)
    assert check_cache["frequency_grid"] is None

    rrc.calculate_raman_tensor(1.0)
    for i in range(3):
        for j in range(3):
            R_l = check_cache[f"{'xyz'[i]}{'xyz'[j]}"]
            assert R_l is not None
            assert R_l[0] == pytest.approx(0.0 + 1j * 0.)
            assert R_l[1] == pytest.approx(0.0 + 1j * 0.)

        if j > i:
            # need to make sure momentum matrix is perfectly hermitian too
            Rother_l = check_cache[f"{'xyz'[j]}{'xyz'[i]}"]
            assert R_l[2].real == pytest.approx(Rother_l[2].real, rel=0.1)
            assert R_l[2].imag == pytest.approx(-Rother_l[2].imag, rel=1.0)
