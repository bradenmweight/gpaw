import pytest
import os
import numpy as np
from gpaw import GPAW
import gpaw.wannier90 as w90
from gpaw.wannier.w90 import read_wout_all
from pathlib import Path
from gpaw.spinorbit import soc_eigenstates

@pytest.mark.wannier
@pytest.mark.serial
def test_wannier90(gpw_files, in_tmp_dir):
    o_ai = [[], [0, 1, 2, 3]]
    bands = range(4)
    calc = GPAW(gpw_files['gaas_pw_nosym_wfs'])
    seed = 'GaAs'
    assert calc.wfs.kd.nbzkpts == calc.wfs.kd.nibzkpts

    w90.write_input(calc, orbitals_ai=o_ai,
                    bands=bands,
                    seed=seed,
                    num_iter=1000,
                    plot=False)
    try:
        os.system('wannier90.x -pp ' + seed)
    except FileNotFoundError:
        return  # no wannier90.x executable
    w90.write_projections(calc, orbitals_ai=o_ai, seed=seed)
    w90.write_eigenvalues(calc, seed=seed)
    w90.write_overlaps(calc, seed=seed)
    os.system('wannier90.x ' + seed)
    with (Path('GaAs.wout')).open() as fd:
        w = read_wout_all(fd)
    centers = np.sum(np.array(w['centers']), axis=0)
    print('centers:', centers)
    centers_correct = np.array([5.68, 5.68, 5.68])
    assert np.allclose(centers, centers_correct)
    spreads = np.sum(np.array(w['spreads']))
    assert spreads == pytest.approx(9.9733, abs=0.002)


@pytest.mark.wannier
@pytest.mark.serial
def test_wannier90_soc(gpw_files, in_tmp_dir):
    calc = GPAW(gpw_files['fe_pw_nosym_wfs'])
    soc = soc_eigenstates(calc)
    seed = 'Fe'
    assert calc.wfs.kd.nbzkpts == calc.wfs.kd.nibzkpts
    
    w90.write_input(calc,
                bands=range(9),
                spinors=True,
                num_iter=200,
                dis_num_iter=500,
                dis_mix_ratio=1.0,
                seed=seed)
    try:
        os.system('wannier90.x -pp ' + seed)
    except FileNotFoundError:
        return  # no wannier90.x executable
    w90.write_projections(calc,
                      seed=seed, soc=soc)
    w90.write_eigenvalues(calc, seed=seed, soc=soc)
    w90.write_overlaps(calc, seed=seed, soc=soc)

    os.system('wannier90.x ' + seed)
    
    with (Path('Fe.wout')).open() as fd:
        w = read_wout_all(fd)
    centers = np.sum(np.array(w['centers']), axis=0)
    centers_correct = [12.9034  , 12.981219, 12.932828]
    assert np.allclose(centers, centers_correct)
    spreads = np.sum(np.array(w['spreads']))
    assert spreads == pytest.approx(19.69782, abs=0.002)

