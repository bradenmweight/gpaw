import pytest
import numpy as np

from gpaw.mpi import world
from gpaw.nlopt.basic import NLOData
from gpaw.nlopt.shift import get_shift


def test_shift_spinpol(mme_files):
    shift_values = np.array([98.0118706 + 0.j,
                             98.89282526 + 0.j,
                             99.78333556 + 0.j,
                             100.68356521 + 0.j,
                             101.59368155 + 0.j,
                             102.51385566 + 0.j,
                             103.44426247 + 0.j,
                             104.38508086 + 0.j,
                             105.33649376 + 0.j,
                             106.29868833 + 0.j,
                             107.271856 + 0.j])

    freqs = np.linspace(2.3, 2.4, 11)
    shift = {}
    for spinpol in ['spinpaired', 'spinpol']:
        tag = '_spinpol' if spinpol == 'spinpol' else ''

        # Get pre-calculated nlodata from SiC fixtures
        nlodata = NLOData.load(mme_files[f'sic_pw{tag}'], comm=world)

        # Calculate 'xyz' tensor element of shift spectra
        get_shift(nlodata, freqs=freqs, eta=0.025, pol='xyz',
                  out_name=f'shift_xyz{tag}.npy')
        world.barrier()

        # Load the calculated SHG spectra (in units of nm/V)
        shift[spinpol] = np.load(f'shift_xyz{tag}.npy')[1] * 1e9
        assert shift[spinpol] == pytest.approx(shift_values, abs=5e-2)

    # Assert that the difference between spectra from spinpaired and
    # spinpolarised calculations is small

    # Absolute error
    assert shift['spinpol'].real == pytest.approx(
        shift['spinpaired'].real, abs=2e-2)

    # Relative error
    assert shift['spinpol'].real == pytest.approx(
        shift['spinpaired'].real, rel=1e-3)

    # Imaginary value should've been removed
    assert shift['spinpol'].imag == pytest.approx(
        shift['spinpaired'].imag, abs=1e-10)
