"""
Calculate the magnetic response in iron using ALDA.

Fast test, where the kernel is scaled to fulfill the Goldstone theorem.
"""

# Workflow modules
import pytest
import numpy as np

# Script modules
import time

from ase.build import bulk
from ase.dft.kpoints import monkhorst_pack
from ase.parallel import parprint

from gpaw import GPAW, PW
from gpaw.response.context import ResponseContext
from gpaw.response.groundstate import ResponseGroundStateAdapter
from gpaw.response.tms import TransverseMagneticSusceptibility
from gpaw.response.susceptibility import read_macroscopic_component
from gpaw.test import findpeak, equal
from gpaw.mpi import world


@pytest.mark.kspair
@pytest.mark.response
def test_response_iron_sf_gssALDA(in_tmp_dir):
    # ------------------- Inputs ------------------- #

    # Part 1: ground state calculation
    xc = 'LDA'
    kpts = 4
    nb = 6
    pw = 300
    a = 2.867
    mm = 2.21

    # Part 2: magnetic response calculation
    q_qc = [[0.0, 0.0, 0.0], [0.0, 0.0, 1. / 4.]]  # Two q-points along G-N
    frq_qw = [np.linspace(-0.080, 0.120, 26), np.linspace(0.100, 0.300, 26)]
    fxc = 'ALDA'
    fxc_scaling = [True, None, 'fm']
    ecut = 300
    eta = 0.01
    if world.size > 1:
        nblocks = 2
    else:
        nblocks = 1

    # ------------------- Script ------------------- #

    # Part 1: ground state calculation

    t1 = time.time()

    Febcc = bulk('Fe', 'bcc', a=a)
    Febcc.set_initial_magnetic_moments([mm])

    calc = GPAW(xc=xc,
                mode=PW(pw),
                kpts=monkhorst_pack((kpts, kpts, kpts)),
                nbands=nb,
                symmetry={'point_group': False},
                parallel={'domain': 1})

    Febcc.calc = calc
    Febcc.get_potential_energy()
    calc.write('Fe', 'all')
    t2 = time.time()

    # Part 2: magnetic response calculation
    context = ResponseContext()
    gs = ResponseGroundStateAdapter.from_gpw_file('Fe', context=context)
    fxckwargs = {'rshelmax': None, 'fxc_scaling': fxc_scaling}
    tms = TransverseMagneticSusceptibility(gs,
                                           context=context,
                                           fxc=fxc,
                                           eta=eta,
                                           ecut=ecut,
                                           fxckwargs=fxckwargs,
                                           gammacentered=True,
                                           nblocks=nblocks)

    for q in range(2):
        tms.get_macroscopic_component(
            '+-', q_qc[q], frq_qw[q],
            filename='iron_dsus' + '_%d.csv' % (q + 1))
        tms.context.write_timer()

    t3 = time.time()

    parprint('Ground state calculation took', (t2 - t1) / 60, 'minutes')
    parprint('Excited state calculation took', (t3 - t2) / 60, 'minutes')

    world.barrier()

    # Part 3: identify magnon peaks in scattering function
    w1_w, chiks1_w, chi1_w = read_macroscopic_component('iron_dsus_1.csv')
    w2_w, chiks2_w, chi2_w = read_macroscopic_component('iron_dsus_2.csv')

    print(w1_w, -chi1_w.imag)
    print(w2_w, -chi2_w.imag)

    wpeak1, Ipeak1 = findpeak(w1_w, -chi1_w.imag)
    wpeak2, Ipeak2 = findpeak(w2_w, -chi2_w.imag)

    mw1 = wpeak1 * 1000
    mw2 = wpeak2 * 1000

    # Part 4: compare new results to test values
    test_fxcs = 1.033
    test_mw1 = -0.03  # meV
    test_mw2 = 176.91  # meV
    test_Ipeak1 = 71.20  # a.u.
    test_Ipeak2 = 44.46  # a.u.

    # fxc_scaling:
    equal(fxc_scaling[1], test_fxcs, 0.005)

    # Magnon peak:
    equal(mw1, test_mw1, 0.1)
    equal(mw2, test_mw2, eta * 650)

    # Scattering function intensity:
    equal(Ipeak1, test_Ipeak1, 5)
    equal(Ipeak2, test_Ipeak2, 5)
