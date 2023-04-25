"""
Calculate the magnetic response in iron using ALDA.

Tests whether the magnon energies and scattering intensities
have changed for:
 * Different kernel calculation strategies
 * Different chi0 transitions summation strategies
"""

import numpy as np
import pytest

from gpaw import GPAW
from gpaw.test import findpeak
from gpaw.mpi import world

from gpaw.response import ResponseGroundStateAdapter
from gpaw.response.chiks import ChiKSCalculator
from gpaw.response.susceptibility import ChiFactory
from gpaw.response.fxc_kernels import AdiabaticFXCCalculator
from gpaw.response.pair_functions import read_pair_function


pytestmark = pytest.mark.skipif(world.size < 4, reason='world.size < 4')


@pytest.mark.kspair
@pytest.mark.response
def test_response_iron_sf_ALDA(in_tmp_dir, gpw_files, scalapack):
    # ---------- Inputs ---------- #

    # Magnetic response calculation
    q_c = [0.0, 0.0, 1 / 4.]
    fxc = 'ALDA'
    ecut = 300
    eta = 0.1

    # Test different kernel, summation and symmetry strategies
    # rshelmax, rshewmin, bandsummation, bundle_integrals, disable_syms
    strat_sd = [(None, None, 'pairwise', True, False),
                (-1, 0.001, 'pairwise', True, False),
                (-1, 0.001, 'pairwise', False, False),
                (-1, 0.000001, 'pairwise', True, False),
                (-1, 0.000001, 'double', True, False),
                (-1, 0.000001, 'double', True, True),
                (-1, None, 'pairwise', True, False),
                (3, None, 'pairwise', True, False)]
    frq_sw = [np.linspace(0.320, 0.480, 21),
              np.linspace(0.420, 0.580, 21),
              np.linspace(0.420, 0.580, 21),
              np.linspace(0.420, 0.580, 21),
              np.linspace(0.420, 0.580, 21),
              np.linspace(0.420, 0.580, 21),
              np.linspace(0.420, 0.580, 21),
              np.linspace(0.420, 0.580, 21),
              np.linspace(0.420, 0.580, 21)]

    # ---------- Script ---------- #

    calc = GPAW(gpw_files['fe_pw_wfs'], parallel=dict(domain=1))
    nbands = calc.parameters.convergence['bands']
    gs = ResponseGroundStateAdapter(calc)

    for s, ((rshelmax, rshewmin, bandsummation, bundle_integrals,
             disable_syms), frq_w) in enumerate(zip(strat_sd, frq_sw)):
        complex_frequencies = frq_w + 1.j * eta
        chiks_calc = ChiKSCalculator(gs,
                                     nbands=nbands,
                                     ecut=ecut,
                                     bandsummation=bandsummation,
                                     bundle_integrals=bundle_integrals,
                                     disable_point_group=disable_syms,
                                     disable_time_reversal=disable_syms,
                                     nblocks=2)
        fxc_calculator = AdiabaticFXCCalculator.from_rshe_parameters(
            gs, chiks_calc.context,
            rshelmax=rshelmax,
            rshewmin=rshewmin)
        chi_factory = ChiFactory(chiks_calc, fxc_calculator)
        chi = chi_factory('+-', q_c, complex_frequencies, fxc=fxc)
        chi.write_macroscopic_component('iron_dsus' + '_G%d.csv' % (s + 1))
        chi_factory.context.write_timer()

    world.barrier()

    # Identify magnon peak in scattering functions
    w1_w, chi1_w = read_pair_function('iron_dsus_G1.csv')
    w2_w, chi2_w = read_pair_function('iron_dsus_G2.csv')
    w3_w, chi3_w = read_pair_function('iron_dsus_G3.csv')
    w4_w, chi4_w = read_pair_function('iron_dsus_G4.csv')
    w5_w, chi5_w = read_pair_function('iron_dsus_G5.csv')
    w6_w, chi6_w = read_pair_function('iron_dsus_G6.csv')
    w7_w, chi7_w = read_pair_function('iron_dsus_G7.csv')
    w8_w, chi8_w = read_pair_function('iron_dsus_G8.csv')
    print(-chi1_w.imag)
    print(-chi2_w.imag)

    wpeak1, Ipeak1 = findpeak(w1_w, -chi1_w.imag)
    wpeak2, Ipeak2 = findpeak(w2_w, -chi2_w.imag)
    wpeak3, Ipeak3 = findpeak(w3_w, -chi3_w.imag)
    wpeak4, Ipeak4 = findpeak(w4_w, -chi4_w.imag)
    wpeak5, Ipeak5 = findpeak(w5_w, -chi5_w.imag)
    wpeak6, Ipeak6 = findpeak(w6_w, -chi6_w.imag)
    wpeak7, Ipeak7 = findpeak(w7_w, -chi7_w.imag)
    wpeak8, Ipeak8 = findpeak(w8_w, -chi8_w.imag)

    mw1 = wpeak1 * 1000
    mw2 = wpeak2 * 1000
    mw3 = wpeak3 * 1000
    mw4 = wpeak4 * 1000
    mw5 = wpeak5 * 1000
    mw6 = wpeak6 * 1000
    mw7 = wpeak7 * 1000
    mw8 = wpeak8 * 1000

    # Compare new results to test values
    print(mw1, mw2, mw4, Ipeak1, Ipeak2, Ipeak4)
    test_mw1 = 402.  # meV
    test_mw2 = 490.  # meV
    test_mw4 = 490.  # meV
    test_Ipeak1 = 4.10  # a.u.
    test_Ipeak2 = 5.05  # a.u.
    test_Ipeak4 = 5.04  # a.u.

    # Different kernel strategies should remain the same
    # Magnon peak:
    assert mw1 == pytest.approx(test_mw1, abs=22.)
    assert mw2 == pytest.approx(test_mw2, abs=22.)
    assert mw4 == pytest.approx(test_mw4, abs=22.)

    # Scattering function intensity:
    assert Ipeak1 == pytest.approx(test_Ipeak1, abs=0.5)
    assert Ipeak2 == pytest.approx(test_Ipeak2, abs=0.5)
    assert Ipeak4 == pytest.approx(test_Ipeak4, abs=0.5)

    # The bundled and unbundled integration methods should give the same
    assert mw2 == pytest.approx(mw3, abs=2.)
    assert Ipeak2 == pytest.approx(Ipeak3, abs=0.1)

    # The two transitions summation strategies should give identical results
    assert mw4 == pytest.approx(mw5, abs=20.)
    assert Ipeak4 == pytest.approx(Ipeak5, abs=0.5)

    # Toggling symmetry should preserve the result
    assert mw6 == pytest.approx(mw4, abs=20.)
    assert Ipeak6 == pytest.approx(Ipeak4, abs=0.5)

    # Including vanishing coefficients should not matter for the result
    assert mw7 == pytest.approx(mw4, abs=2.)
    assert Ipeak7 == pytest.approx(Ipeak4, abs=0.1)
    assert mw8 == pytest.approx(mw2, abs=2.)
    assert Ipeak8 == pytest.approx(Ipeak2, abs=0.1)
