import pytest
import os
import numpy as np
from gpaw import GPAW
from gpaw.wannier90 import Wannier90
from gpaw.wannier.w90 import read_wout_all
from pathlib import Path
from gpaw.spinorbit import soc_eigenstates
from subprocess import PIPE, run


def out():
    result = run('wannier90.x --version',
                 stdout=PIPE,
                 stderr=PIPE,
                 universal_newlines=True,
                 shell=True)
    return result.stdout


@pytest.mark.wannier
@pytest.mark.serial
@pytest.mark.skipif(': 3.' not in out(),
                    reason="requires at least Wannier90 version 3.0")
@pytest.mark.parametrize('mode', ['sym', 'nosym'])
def test_wannier90(gpw_files, mode, in_tmp_dir):
    o_ai = [[], [0, 1, 2, 3]]
    bands = range(4)
    
    if mode == 'sym':
        calc = GPAW(gpw_files['gaas_pw'])
        assert calc.wfs.kd.nbzkpts > calc.wfs.kd.nibzkpts
    else:
        calc = GPAW(gpw_files['gaas_pw_nosym'])
        assert calc.wfs.kd.nbzkpts == calc.wfs.kd.nibzkpts

    seed = f'GaAs_{mode}'

    wannier = Wannier90(calc,
                        seed=seed,
                        bands=bands,
                        orbitals_ai=o_ai)
    
    wannier.write_input(num_iter=1000,
                        plot=False)

    os.system('wannier90.x -pp ' + seed)
    wannier.write_projections()
    wannier.write_eigenvalues()
    wannier.write_overlaps()
    os.system('wannier90.x ' + seed)
    with (Path(f'{seed}.wout')).open() as fd:
        w = read_wout_all(fd)
    centers = np.sum(np.array(w['centers']), axis=0)
    print('centers:', centers)
    centers_correct = np.array([5.68, 5.68, 5.68])
    assert np.allclose(centers, centers_correct, atol=1e-3)
    spreads = np.sum(np.array(w['spreads']))
    assert spreads == pytest.approx(9.9733, abs=0.002)


@pytest.mark.wannier
@pytest.mark.serial
@pytest.mark.skipif(': 3.' not in out(),
                    reason="requires at least Wannier90 version 3.0")
def test_wannier90_soc(gpw_files, in_tmp_dir):
    calc = GPAW(gpw_files['fe_pw_nosym'])
    soc = soc_eigenstates(calc)
    seed = 'Fe'
    assert calc.wfs.kd.nbzkpts == calc.wfs.kd.nibzkpts

    wannier = Wannier90(calc,
                        seed=seed,
                        bands=range(9),
                        spinors=True)

    wannier.write_input(num_iter=200,
                    dis_num_iter=500,
                    dis_mix_ratio=1.0)
    os.system('wannier90.x -pp ' + seed)
    wannier.write_projections()
    wannier.write_eigenvalues()
    wannier.write_overlaps()

    os.system('wannier90.x ' + seed)

    with (Path('Fe.wout')).open() as fd:
        w = read_wout_all(fd)
    centers = np.sum(np.array(w['centers']), axis=0)
    centers_correct = [12.9, 12.9, 12.9]
    assert np.allclose(centers, centers_correct, atol=0.1)
    spreads = np.sum(np.array(w['spreads']))
    assert spreads == pytest.approx(20.1, abs=0.6)
