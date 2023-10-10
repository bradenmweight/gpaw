import numpy as np
import pytest
from gpaw import GPAW, PW, FermiDirac, Mixer
from gpaw.response.modelinteraction import initialize_w_model
from gpaw.response.chi0 import Chi0
from gpaw.wannier90 import Wannier90
import os
from ase.build import bulk
from ase import Atoms

@pytest.mark.serial  # W90 does not work in parallel
@pytest.mark.parametrize('symm', [True, False])
def test_w90(in_tmp_dir, gpw_files, symm):
    k = 4
    cell = bulk('Ga', 'fcc', a=5.68).cell
    a = Atoms('GaAs', cell=cell, pbc=True,
              scaled_positions=((0, 0, 0), (0.25, 0.25, 0.25)))

    if symm:
        calc = GPAW(gpw_files['gaas_pw'])
    else:
        calc = GPAW(gpw_files['gaas_pw_nosym'])

    
    seed = 'GaAs'
    w90 = Wannier90(calc, orbitals_ai=[[], [0, 1, 2, 3]],
                    bands=range(4),
                    seed=seed)
    w90.write_input(num_iter=1000,
                    plot=True,
                    write_u_matrices=True)

    w90.write_wavefunctions()
    os.system('wannier90.x -pp ' + seed)

    w90.write_projections()
    w90.write_eigenvalues()
    w90.write_overlaps()

    os.system('wannier90.x ' + seed)

    omega = np.array([0])
    chi0calc = Chi0(calc, frequencies=omega, hilbert=False,ecut=100, txt='test.log',intraband=False)
    txt='out.txt'
    truncation = None
    Wm = initialize_w_model(calc, chi0calc)
    Wwann = Wm.calc_in_Wannier(chi0calc,Uwan=seed,bandrange=[0,4])
    (1.059069004580567+2.873325682818419e-17j)
    (0.7744817912919332+1.699625225785978e-17j)
    (0.7744817912872252-2.4259190263322844e-18j)
    (0.7744817912816624-1.5508496535341566e-17j)
    (0.4058975488523744+3.711483001934005e-18j)
    (0.7548791410801549-0.0017578979720108043j)

    assert Wwann[0, 0, 0, 0, 0] == pytest.approx(1.059, abs=0.003)
    assert Wwann[0, 1, 1, 1, 1] == pytest.approx(0.774, abs=0.003)
    assert Wwann[0, 2, 2, 2, 2] == pytest.approx(0.774, abs=0.003)
    assert Wwann[0, 3, 3, 3, 3] == pytest.approx(0.774, abs=0.003)
    assert Wwann[0, 3, 3, 0, 0] == pytest.approx(0.406, abs=0.003)
    assert Wwann[0, 3, 0, 3, 0] == pytest.approx(0.755, abs=0.003)
