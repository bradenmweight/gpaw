import pytest
from gpaw.mpi import world
import numpy as np
from gpaw.response.g0w0 import G0W0
import pickle


@pytest.mark.response
def test_do_GW_too(in_tmp_dir, gpw_files, scalapack):
    ref_gap = 4.7747
    gw = G0W0(gpw_files['bn_pw_wfs'],
              bands=(3, 5),
              nbands=9,
              nblocks=1,
              xc='rALDA',
              method='G0W0',
              ecut=40,
              fxc_mode='GWP',
              do_GW_too=True)

    gw.calculate()

    world.barrier()

    with open('gw_results_GW.pckl', 'rb') as handle:
        results_GW = pickle.load(handle)
    calculated_gap = np.min(results_GW['qp'][0, :, 1])\
        - np.max(results_GW['qp'][0, :, 0])
    assert calculated_gap == pytest.approx(ref_gap, abs=0.001)
