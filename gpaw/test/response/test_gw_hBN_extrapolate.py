""" Tests extrapolation to infinite energy cutoff + block paralellization.
It takes ~109 s on one core"""

import pytest
from gpaw.mpi import world
from gpaw.utilities import compiled_with_sl
import numpy as np

from ase.lattice.hexagonal import Graphene

from gpaw import GPAW, FermiDirac
from gpaw.test import equal
from gpaw.response.g0w0 import G0W0

pytestmark = pytest.mark.skipif(
    world.size != 1 and not compiled_with_sl(),
    reason='world.size != 1 and not compiled_with_sl()')


@pytest.mark.response
def test_response_gw_hBN_extrapolate(in_tmp_dir):
    if 1:
        calc = GPAW(mode='pw',
                    xc='PBE',
                    nbands=16,
                    setups={'Mo': '6'},
                    eigensolver='rmm-diis',
                    occupations=FermiDirac(0.001),
                    parallel={'domain': 1},
                    kpts={'size': (6, 6, 1), 'gamma': True})

        layer = Graphene(symbol='B',
                         latticeconstant={'a': 2.5, 'c': 1.0},
                         size=(1, 1, 1))
        layer[0].symbol = 'N'
        layer.pbc = (1, 1, 0)
        layer.center(axis=2, vacuum=4.0)

        layer.calc = calc
        layer.get_potential_energy()

        nbecut = 50
        from ase.units import Bohr, Hartree
        vol = layer.get_volume() / Bohr**3
        nbands = int(vol * (nbecut / Hartree)**1.5 * 2**0.5 / 3 / np.pi**2)
        calc.diagonalize_full_hamiltonian(nbands)
        calc.write('hBN.gpw', mode='all')

    gw = G0W0('hBN.gpw',
              'gw-hBN',
              ecut=50,
              frequencies={'type': 'nonlinear',
                           'domega0': 0.1},
              eta=0.2,
              truncation='2D',
              kpts=[0],
              bands=(3, 5),
              ecut_extrapolation=[30, 40, 50],
              nblocksmax=True)

    e_qp = gw.calculate()['qp'][0, 0]

    ev = -4.381
    ec = 3.71013806
    equal(e_qp[0], ev, 0.01)
    equal(e_qp[1], ec, 0.01)
