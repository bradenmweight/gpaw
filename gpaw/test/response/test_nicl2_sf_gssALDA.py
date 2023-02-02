import numpy as np
import pytest

from gpaw.mpi import world
from gpaw.test import findpeak

from gpaw.response import ResponseContext, ResponseGroundStateAdapter
from gpaw.response.frequencies import ComplexFrequencyDescriptor
from gpaw.response.chiks import ChiKSCalculator
from gpaw.response.localft import LocalFTCalculator
from gpaw.response.fxc_kernels import FXCScaling
from gpaw.response.susceptibility import ChiFactory
from gpaw.response.df import read_response_function


@pytest.mark.kspair
@pytest.mark.response
def test_nicl2_magnetic_response(in_tmp_dir, gpw_files):
    # ---------- Inputs ---------- #

    q_qc = [[0., 0., 0.],
            [1. / 3., 1. / 3., 0.]]
    fxc = 'ALDA'
    fxc_scaling = FXCScaling('fm')
    rshelmax = -1
    rshewmin = 1e-8
    ecut = 200
    frq_w = np.linspace(-0.09, 0.06, 21)
    eta = 0.05
    zd = ComplexFrequencyDescriptor.from_array(frq_w + 1.j * eta)
    if world.size > 1:
        nblocks = 2
    else:
        nblocks = 1

    # ---------- Script ---------- #

    # Magnetic response calculation
    context = ResponseContext()
    gs = ResponseGroundStateAdapter.from_gpw_file(gpw_files['nicl2_pw_wfs'],
                                                  context=context)
    chiks_calc = ChiKSCalculator(gs,
                                 ecut=ecut,
                                 gammacentered=True,
                                 nblocks=nblocks)
    localft_calc = LocalFTCalculator.from_rshe_parameters(
        gs, chiks_calc.context,
        rshelmax=rshelmax,
        rshewmin=rshewmin)
    chi_factory = ChiFactory(chiks_calc)

    for q, q_c in enumerate(q_qc):
        filename = 'nicl2_macro_tms_q%d.csv' % q
        txt = 'nicl2_macro_tms_q%d.txt' % q
        chi = chi_factory('+-', q_c, zd,
                          fxc=fxc,
                          localft_calc=localft_calc,
                          fxc_scaling=fxc_scaling,
                          txt=txt)
        chi.write_macroscopic_component(filename)

    chi_factory.context.write_timer()
    world.barrier()

    # Identify magnon peaks and extract kernel scaling
    w0_w, _, chi0_w = read_response_function('nicl2_macro_tms_q0.csv')
    w1_w, _, chi1_w = read_response_function('nicl2_macro_tms_q1.csv')

    wpeak0, Ipeak0 = findpeak(w0_w, -chi0_w.imag / np.pi)
    wpeak1, Ipeak1 = findpeak(w1_w, -chi1_w.imag / np.pi)
    mw0 = wpeak0 * 1e3  # meV
    mw1 = wpeak1 * 1e3  # meV

    assert fxc_scaling.has_scaling
    fxcs = fxc_scaling.get_scaling()

    if world.rank == 0:
        # import matplotlib.pyplot as plt
        # plt.plot(w0_w, -chi0_w.imag / np.pi)
        # plt.plot(w1_w, -chi1_w.imag / np.pi)
        # plt.show()

        print(fxcs, mw0, mw1, Ipeak0, Ipeak1)
    world.barrier()

    # Compare new results to test values
    test_fxcs = 0.68648
    test_mw0 = -1.8  # meV
    test_mw1 = -44.0  # meV
    test_Ipeak0 = 0.5528  # a.u.
    test_Ipeak1 = 0.2332  # a.u.

    # Test fxc scaling
    assert fxcs == pytest.approx(test_fxcs, abs=0.005)

    # Test magnon peaks
    assert test_mw0 == pytest.approx(test_mw0, abs=2.0)
    assert test_mw1 == pytest.approx(test_mw1, abs=2.0)

    # Test peak intensities
    assert test_Ipeak0 == pytest.approx(test_Ipeak0, abs=0.01)
    assert test_Ipeak1 == pytest.approx(test_Ipeak1, abs=0.01)
