"""Test GW band-gaps for Si."""

import pytest
from ase.build import bulk

from gpaw import GPAW
from gpaw.mpi import world
from gpaw.response.g0w0 import G0W0


def generate_si_systems():
    a = 5.43
    si1 = bulk('Si', 'diamond', a=a)
    si2 = si1.copy()
    si2.positions -= a / 8

    return [si1, si2]


def run(gpw_filename, nblocks, integrate_gamma):
    # This tests checks the actual numerical accuracy which is asserted below
    calc = GPAW(gpw_filename)
    e = calc.get_potential_energy()

    gw = G0W0(gpw_filename, 'gw_None',
              nbands=8, integrate_gamma=integrate_gamma,
              kpts=[(0, 0, 0), (0.5, 0.5, 0)],  # Gamma, X
              ecut=40, nblocks=nblocks,
              frequencies={'type': 'nonlinear',
                           'domega0': 0.1, 'omegamax': None},
              eta=0.2, relbands=(-1, 2))
    results = gw.calculate()

    G, X = results['eps'][0]
    output = [e, G[0], G[1] - G[0], X[1] - G[0], X[2] - X[1]]
    G, X = results['qp'][0]
    output += [G[0], G[1] - G[0], X[1] - G[0], X[2] - X[1]]

    return output


reference = {0: pytest.approx([-9.253, 5.442, 2.389, 0.403, 0.000,
                               6.261, 3.570, 1.323, 0.001], abs=0.0035),
            'WS': pytest.approx([-9.253, 5.442, 2.389, 0.403, 0.000,
                                 6.284, 3.551, 1.285, 0.001], abs=0.0035)}

@pytest.mark.response
@pytest.mark.slow
@pytest.mark.parametrize('si', [0, 1])
@pytest.mark.parametrize('integrate_gamma', [0, 'WS'])
@pytest.mark.parametrize('symm', ['all', 'no', 'tr', 'pg'])
@pytest.mark.parametrize('nblocks',
                         [x for x in [1, 2, 4, 8] if x <= world.size])
def test_response_gwsi(in_tmp_dir, si, symm, nblocks, integrate_gamma, 
                       scalapack, gpw_files, gpaw_new):
    if gpaw_new and world.size > 1:
        pytest.skip('Hybrids not working in parallel with GPAW_NEW=1')
    filename = gpw_files[f'si_gw_a{si}_{symm}']
    assert run(filename, nblocks, integrate_gamma) == reference[integrate_gamma]


@pytest.mark.response
@pytest.mark.ci
@pytest.mark.parametrize('si', [0, 1])
@pytest.mark.parametrize('symm', ['all'])
def test_small_response_gwsi(in_tmp_dir, si, symm, scalapack,
                             gpw_files, gpaw_new):
    if gpaw_new and world.size > 1:
        pytest.skip('Hybrids not working in parallel with GPAW_NEW=1')
    filename = gpw_files[f'si_gw_a{si}_{symm}']
    assert run(filename, 1, 0) == reference[0]


@pytest.mark.response
@pytest.mark.ci
def test_few_freq_response_gwsi(in_tmp_dir, scalapack,
                                gpw_files, gpaw_new):
    if gpaw_new and world.size > 1:
        pytest.skip('Hybrids not working in parallel with GPAW_NEW=1')

    if world.size > 1:
        nblocks = 2
    else:
        nblocks = 1

    # This test has very few frequencies and tests that the code doesn't crash.
    filename = gpw_files['si_gw_a0_all']
    gw = G0W0(filename, 'gw_0.2',
              nbands=8, integrate_gamma=0,
              kpts=[(0, 0, 0), (0.5, 0.5, 0)],  # Gamma, X
              ecut=40, nblocks=nblocks,
              frequencies={'type': 'nonlinear',
                           'domega0': 0.1, 'omegamax': 0.2},
              eta=0.2, relbands=(-1, 2))
    gw.calculate()
