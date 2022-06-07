"""Test the site kernel calculation functionality of the response code"""
from gpaw import GPAW, PW
from gpaw.response.mft import StaticChiKSFactory
from gpaw.response.site_kernels import site_kernel_interface
from gpaw.response.susceptibility import get_pw_coordinates
from ase.build import bulk
import numpy as np


# ---------- Main test ---------- #


def test_site_kernels():
    Co_hcp_test()


# ---------- Actual tests ---------- #


def Co_hcp_test():
    # ---------- Inputs ---------- #

    # Part 1: Perform ground state calculation
    # This step should be made redudant! XXX
    xc = 'LDA'
    kpts = 4
    pw = 200
    a = 2.5071
    c = 4.0695
    mm = 1.6

    # Part 2: Generate PWDescriptor
    ecut = 50

    # Part 3: Calculate site kernels
    rc_m = np.array([2.0, 3.0])  # radii in Aa
    zc_m = ['diameter', 'unit cell']

    # Part 4: Check the calculated kernels

    # ---------- Script ---------- #

    # Part 1: Perform ground state calculation to get PWDescriptor

    atoms = bulk('Co', 'hcp', a=a, c=c)
    atoms.set_initial_magnetic_moments([mm, mm])

    calc = GPAW(xc=xc,
                spinpol=True,
                mode=PW(pw),
                kpts={'size': (kpts, kpts, kpts),
                      'gamma': True}
                )

    atoms.calc = calc
    atoms.get_potential_energy()

    # Part 2: Generate PWDescriptor
    chiksf = StaticChiKSFactory(calc, ecut=ecut)
    pd0 = chiksf.get_PWDescriptor([0, 0, 0])

    # Part 3: Calculate site kernels
    siteposition_mv = atoms.get_positions()
    nG = len(get_pw_coordinates(pd0))

    # Compute site-kernels
    Kuc_GGm = site_kernel_interface(pd0, siteposition_mv,
                                    shapes_m='unit cell')
    Ksph_GGm = site_kernel_interface(pd0, siteposition_mv,
                                     shapes_m='sphere', rc_m=rc_m)
    Kcyl_GGm = site_kernel_interface(pd0, siteposition_mv,
                                     shapes_m='cylinder', rc_m=rc_m, zc_m=zc_m)

    # Part 4: Check the calculated kernels
    V0 = atoms.get_volume()  # Volume of unit cell in Aa^3

    # Calculate integration volumes in Aa^3
    Vsphere_m = 4 / 3 * np.pi * rc_m**3
    height_m = np.array([2 * rc_m[0], np.sum(atoms.cell[:, -1])])
    Vcylinder_m = np.pi * rc_m**2 * height_m

    # Check shape of K-arrays
    assert Kuc_GGm.shape == (nG, nG, 1)
    assert Ksph_GGm.shape == (nG, nG, 2)
    assert Kcyl_GGm.shape == (nG, nG, 2)

    # Check that K_00(q=0) gives Vint / V0 (fractional integration volume)
    assert abs(Kuc_GGm[0, 0, 0] - 1.) < 1.e-8
    assert np.allclose(Ksph_GGm[0, 0, :], Vsphere_m / V0)
    assert np.allclose(Kcyl_GGm[0, 0, :], Vcylinder_m / V0)
