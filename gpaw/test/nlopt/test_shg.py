import numpy as np
import pytest

from gpaw.mpi import world
from gpaw.nlopt.basic import NLOData
from gpaw.nlopt.shg import get_shg


def test_shg_spinpol(mme_files):
    shg_values = np.array([-0.77053399 - 0.37041593j,
                           -0.87903174 - 0.4177294j,
                           -1.00791251 - 0.51051291j,
                           -1.15465962 - 0.66642326j,
                           -1.30812094 - 0.92114822j,
                           -1.42138133 - 1.33424513j,
                           -1.34649601 - 1.96084827j,
                           -0.78891819 - 2.66240386j,
                           0.27801137 - 2.8572836j,
                           1.12315952 - 2.30446868j,
                           1.38569995 - 1.59698796j])

    freqs = np.linspace(2.3, 2.4, 11)
    shg = {}
    for spinpol in ['spinpaired', 'spinpol']:
        tag = '_spinpol' if spinpol == 'spinpol' else ''

        # Get pre-calculated nlodata from SiC fixtures
        nlodata = NLOData.load(mme_files[f'sic_pw{tag}'], comm=world)

        # Calculate 'xyz' tensor element of SHG spectra
        get_shg(nlodata, freqs=freqs, eta=0.025, pol='xyz',
                out_name=f'shg_xyz{tag}.npy')
        world.barrier()

        # Load the calculated SHG spectra (in units of nm/V)
        shg[spinpol] = np.load(f'shg_xyz{tag}.npy')[1] * 1e9
        assert shg[spinpol] == pytest.approx(shg_values, abs=1e-3), \
            np.max(np.abs(shg[spinpol] - shg_values))

    # Assert that the difference between spectra from spinpaired and
    # spinpolarised calculations is small

    # Absolute error
    assert shg['spinpol'].real == pytest.approx(
        shg['spinpaired'].real, abs=5e-4)
    assert shg['spinpol'].imag == pytest.approx(
        shg['spinpaired'].imag, abs=5e-4)

    # Relative error
    assert shg['spinpol'].real == pytest.approx(
        shg['spinpaired'].real, rel=2e-3)
    assert shg['spinpol'].imag == pytest.approx(
        shg['spinpaired'].imag, rel=2e-3)
