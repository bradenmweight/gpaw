import pytest
from gpaw.mpi import world
from gpaw.utilities import compiled_with_sl
import numpy as np
from gpaw.response.g0w0 import G0W0

pytestmark = pytest.mark.skipif(
    world.size != 1 and not compiled_with_sl(),
    reason='world.size != 1 and not compiled_with_sl()')


@pytest.mark.response
def test_do_GW_too(in_tmp_dir, gpw_files):
    ref_result = np.asarray([[[11.30094393, 21.62842077],
                              [5.33751513, 16.06905725],
                              [8.75269938, 22.46579489]]])
    gw = G0W0(gpw_files['bn_pw_wfs'],
              bands=(3, 5),
              nbands=9,
              nblocks=1,
              method='G0W0',
              ecut=40,
              ppa=True)

    results = gw.calculate()
    np.testing.assert_allclose(results['qp'], ref_result, rtol=1e-03)
