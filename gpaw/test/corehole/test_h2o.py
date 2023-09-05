import numpy as np
import pytest

import gpaw.mpi as mpi
from gpaw import GPAW
from gpaw.atom.generator2 import generate
from gpaw.test import equal
from gpaw.xas import XAS


@pytest.mark.later
def test_corehole_h2o(in_tmp_dir, add_cwd_to_setup_paths, gpw_files):
    # Generate setup for oxygen with half a core-hole:
    gen = generate('O', 8, '2s,s,2p,p,d', [1.2], 1.0, None, 2,
                   core_hole='1s,0.5')
    setup = gen.make_paw_setup('hch1s')
    setup.write_xml()

    calc = GPAW(gpw_files['h2o_xas'])
    if mpi.size == 1:
        xas = XAS(calc)
        x, y = xas.get_spectra()
        e1_n = xas.eps_n
        de1 = e1_n[1] - e1_n[0]

    if mpi.size == 1:
        # calc = GPAW('h2o-xas.gpw')
        # poissonsolver=FDPoissonSolver(use_charge_center=True))
        # calc.initialize()
        xas = XAS(calc)
        x, y = xas.get_spectra()
        e2_n = xas.eps_n
        w_n = np.sum(xas.sigma_cn.real**2, axis=0)
        de2 = e2_n[1] - e2_n[0]

        equal(de2, 2.064, 0.005)
        equal(w_n[1] / w_n[0], 2.22, 0.01)

        assert de1 == de2

    if 0:
        import matplotlib.pyplot as plt
        plt.plot(x, y[0])
        plt.plot(x, sum(y))
        plt.show()
