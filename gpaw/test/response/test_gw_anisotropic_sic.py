import pytest
from gpaw.response.g0w0 import G0W0
import numpy as np
from gpaw.mpi import world


@pytest.mark.response
def test_gw_anisotropic(in_tmp_dir, gpw_files, gpaw_new):
    if gpaw_new and world.size > 1:
        pytest.skip('Hybrids not working in parallel with GPAW_NEW=1')
    print(gpw_files)

    gw = G0W0(gpw_files['sic_pw'],
              'gw-test',
              nbands=5,
              ecut=20,
              eta=0.2,
              frequencies={'type': 'nonlinear', 'domega0': 0.3},
              truncation=None,
              kpts=[(-0.125, 0.125, 0.125), (-0.125, -0.125, -0.125)],
              bands=(3, 5))

    e_qp = gw.calculate()['qp']

    print(e_qp)
    assert np.allclose(e_qp, [[[8.48523631, 14.83302236],
                               [7.12242815, 15.76829637]]], atol=0.001)
