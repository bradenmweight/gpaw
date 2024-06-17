import pytest
import numpy as np
from gpaw.response.g0w0 import G0W0
from gpaw.mpi import world


@pytest.mark.response
def test_ppa(in_tmp_dir, gpw_files, scalapack, gpaw_new):
    if gpaw_new and world.size > 1:
        pytest.skip('Hybrids not working in parallel with GPAW_NEW=1')
    ref_result = np.asarray([[[11.30094393, 21.62842077],
                              [5.33751513, 16.06905725],
                              [8.75269938, 22.46579489]]])
    gw = G0W0(gpw_files['bn_pw'],
              bands=(3, 5),
              nbands=9,
              nblocks=1,
              ecut=40,
              ppa=True)

    results = gw.calculate()
    np.testing.assert_allclose(results['qp'], ref_result, rtol=1e-03)
