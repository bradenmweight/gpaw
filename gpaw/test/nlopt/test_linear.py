import pytest
import numpy as np

from gpaw.mpi import world
from gpaw.nlopt.basic import NLOData
from gpaw.nlopt.linear import get_chi_tensor


def test_chi_spinpol(mme_files):
    chi_values = np.array([7.97619464 + 0.06564128j,
                           7.98936984 + 0.06620975j,
                           8.00265916 + 0.06678416j,
                           8.01606379 + 0.06736459j,
                           8.02958494 + 0.06795116j,
                           8.04322385 + 0.06854397j,
                           8.05698176 + 0.06914314j,
                           8.07085997 + 0.06974877j,
                           8.08485976 + 0.07036098j,
                           8.09898247 + 0.07097988j,
                           8.11322944 + 0.07160561j])

    freqs = np.linspace(2.3, 2.4, 11)
    chi_xx = {}
    for spinpol in ['spinpaired', 'spinpol']:
        tag = '_spinpol' if spinpol == 'spinpol' else ''

        # Get pre-calculated nlodata from SiC fixtures
        nlodata = NLOData.load(mme_files[f'sic_pw{tag}'], comm=world)

        # Calculate tensor elements of susceptibility spectra
        get_chi_tensor(nlodata, freqs=freqs,
                       eta=0.05, out_name=f'chi{tag}.npy')
        world.barrier()

        # Load the calculated susceptibility
        chi_xx[spinpol] = np.load(f'chi{tag}.npy')[1]

        assert chi_xx[spinpol] == pytest.approx(chi_values, abs=5e-2)

    # Assert that the difference between spectra from spinpaired and
    # spinpolarised calculations is small

    # Absolute error
    assert chi_xx['spinpol'].real == pytest.approx(
        chi_xx['spinpaired'].real, abs=2e-2)

    # Relative error
    assert chi_xx['spinpol'].real == pytest.approx(
        chi_xx['spinpaired'].real, rel=1e-3)
    assert chi_xx['spinpol'].imag == pytest.approx(
        chi_xx['spinpaired'].imag, abs=1e-4)
