"""Test the site kernel calculation functionality of the response code"""
import numpy as np
import scipy.special as sc

from ase.build import bulk

from gpaw import GPAW, PW
from gpaw.response.site_kernels import (SphericalSiteKernels,
                                        CylindricalSiteKernels,
                                        ParallelepipedicSiteKernels,
                                        sinc,
                                        spherical_geometry_factor,
                                        cylindrical_geometry_factor,
                                        parallelepipedic_geometry_factor)
from gpaw.response.susceptibility import get_pw_coordinates


# ---------- Main test ---------- #


def test_site_kernels():
    spherical_kernel_test()
    cylindrical_kernel_test()
    parallelepipedic_kernel_test()
    Co_hcp_test()


# ---------- Actual tests ---------- #


def spherical_kernel_test():
    """Check the numerics of the spherical kernel"""
    # ---------- Inputs ---------- #

    # Relative wave vector lengths to check (relative to 1/rc)
    Qrel_Q = np.array([0.,
                       np.pi / 2,
                       np.pi])

    # Expected results (assuming |Q| * rc = 1.)
    Vsphere = 4. * np.pi / 3.
    test_K_Q = Vsphere * np.array([1.,
                                   3 / (np.pi / 2)**3.,
                                   3 / (np.pi)**2.])

    # Spherical radii to check
    nr = 5
    rc_r = np.random.rand(nr)

    # Wave vector directions to check
    nd = 41
    Q_dv = 2. * np.random.rand(nd, 3) - 1.
    Q_dv /= np.linalg.norm(Q_dv, axis=1)[:, np.newaxis]  # normalize

    # ---------- Script ---------- #

    # Set up wave vectors
    Q_Qdv = Qrel_Q[:, np.newaxis, np.newaxis] * Q_dv[np.newaxis, ...]

    for rc in rc_r:
        # Calculate site centered geometry factor with rescaled wave vector
        K_Qd = spherical_geometry_factor(Q_Qdv / rc, rc)
        # Check against expected result
        assert np.allclose(K_Qd, rc**3. * test_K_Q[:, np.newaxis])


def cylindrical_kernel_test():
    """Check the numerics of the spherical kernel"""
    # ---------- Inputs ---------- #

    # Relative wave vectors (relative to 1/rc) in radial direction
    Qrhorel_Q1 = np.array([0.] + list(sc.jn_zeros(1, 4)))  # Roots of J1(x)

    # Relative wave vectors (relative to 2/hc) in cylindrical direction
    Qzrel_Q2 = list(np.pi * np.arange(5))  # Roots of sin(x)
    Qzrel_Q2 += list(np.pi * np.arange(4) + np.pi / 2)  # Extrema of sin(x)
    Qzrel_Q2 = np.array(Qzrel_Q2)

    # Expected results for roots of J1 (assuming rc=1. and hc=2.)
    Vcylinder = 2. * np.pi
    nQ2 = 13  # Choose random Q_z r_z
    test_Krho_Q1 = np.array([1., 0., 0., 0., 0.])
    Qzrand_Q2 = 10. * np.random.rand(nQ2)
    sinc_zrand_Q2 = np.sin(Qzrand_Q2) / Qzrand_Q2
    test_Krho_Q1Q2 = Vcylinder * test_Krho_Q1[:, np.newaxis]\
        * sinc_zrand_Q2[np.newaxis, :]

    # Expected results for roots and extrema of sin (assuming rc=1. and hc=2.)
    nQ1 = 15  # Choose random Q_ρ h_c
    test_Kz_Q2 = [1., 0., 0., 0., 0.]  # Nodes in sinc(Q_z h_c)
    test_Kz_Q2 += list(np.array([1., -1., 1., -1.]) / Qzrel_Q2[5:])  # Extrema
    test_Kz_Q2 = np.array(test_Kz_Q2)
    Qrhorand_Q1 = 10. * np.random.rand(nQ1)
    J1term_rhorand_Q1 = 2. * sc.jv(1, Qrhorand_Q1) / Qrhorand_Q1
    test_Kz_Q1Q2 = Vcylinder * J1term_rhorand_Q1[:, np.newaxis]\
        * test_Kz_Q2[np.newaxis, :]

    # Cylinder radii to check
    nr = 5
    rc_r = 3. * np.random.rand(nr)

    # Cylinder height to check
    nh = 3
    hc_h = 4. * np.random.rand(nh)

    # Cylindrical axes to check
    nc = 7
    ez_cv = 2. * np.random.rand(nc, 3) - 1.
    ez_cv /= np.linalg.norm(ez_cv, axis=1)[:, np.newaxis]

    # Wave vector directions in-plane to check. Generated through the cross
    # product of a random direction with the cylindrical axis
    nd = 11
    Qrho_dv = 2. * np.random.rand(nd, 3) - 1.
    Qrho_cdv = np.cross(Qrho_dv[np.newaxis, ...], ez_cv[:, np.newaxis, :])
    Qrho_cdv /= np.linalg.norm(Qrho_cdv, axis=-1)[..., np.newaxis]  # normalize

    # ---------- Script ---------- #

    for rc in rc_r:
        for hc in hc_h:
            # Set up wave vectors for radial tests
            Qrho_cdQ1v = Qrhorel_Q1[np.newaxis, np.newaxis, :, np.newaxis]\
                * Qrho_cdv[..., np.newaxis, :] / rc
            Qrho_cQ2v = Qzrand_Q2[np.newaxis, :, np.newaxis]\
                * ez_cv[:, np.newaxis, :] / (hc / 2.)
            Qrho_cdQ1Q2v = Qrho_cdQ1v[..., np.newaxis, :]\
                + Qrho_cQ2v[:, np.newaxis, np.newaxis, ...]

            # Set up wave vectors for cylindrical tests
            Qz_cdQ1v = Qrhorand_Q1[np.newaxis, np.newaxis, :, np.newaxis]\
                * Qrho_cdv[..., np.newaxis, :] / rc
            Qz_cQ2v = Qzrel_Q2[np.newaxis, :, np.newaxis]\
                * ez_cv[:, np.newaxis, :] / (hc / 2.)
            Qz_cdQ1Q2v = Qz_cdQ1v[..., np.newaxis, :]\
                + Qz_cQ2v[:, np.newaxis, np.newaxis, ...]

            # Test one cylindrical direction at a time
            for ez_v, Qrho_dQ1Q2v, Qz_dQ1Q2v in zip(ez_cv,
                                                    Qrho_cdQ1Q2v, Qz_cdQ1Q2v):
                # Calculate geometry factors
                Krho_dQ1Q2 = cylindrical_geometry_factor(Qrho_dQ1Q2v,
                                                         ez_v, rc, hc)
                Kz_dQ1Q2 = cylindrical_geometry_factor(Qz_dQ1Q2v,
                                                       ez_v, rc, hc)

                # Check against expected result
                assert np.allclose(Krho_dQ1Q2, rc**2. * hc / 2.
                                   * test_Krho_Q1Q2[np.newaxis, ...],
                                   atol=1.e-8)
                assert np.allclose(Kz_dQ1Q2, rc**2. * hc / 2.
                                   * test_Kz_Q1Q2[np.newaxis, ...],
                                   atol=1.e-8)


def parallelepipedic_kernel_test():
    """Check the numerics of the parallelepipedic site kernel."""
    # ---------- Inputs ---------- #

    # Relative wave vectors to check and corresponding sinc(x/2)
    Qrel_Q = np.pi * np.arange(5)
    sinchalf_Q = np.array([1., 2. / np.pi, 0., - 2. / (3. * np.pi), 0.])

    # Random parallelepipedic cell vectors to check
    nC = 9
    cell_Ccv = 2. * np.random.rand(nC, 3, 3) - 1.
    volume_C = np.abs(np.linalg.det(cell_Ccv))
    # Normalize the cell volume
    cell_Ccv /= (volume_C**(1 / 3))[:, np.newaxis, np.newaxis]

    # Transverse wave vector components to check. Generated through the cross
    # product of a random direction with the first cell axis.
    v0_Cv = cell_Ccv[:, 0, :].copy()
    v0_C = np.linalg.norm(v0_Cv, axis=-1)  # Length of primary vector
    v0n_Cv = v0_Cv / v0_C[:, np.newaxis]  # Normalize
    nd = 11
    Q_dv = 2. * np.random.rand(nd, 3) - 1.
    Q_dv[0, :] = np.array([0., 0., 0.])  # Check also parallel Q-vector
    Q_Cdv = np.cross(Q_dv[np.newaxis, ...], v0_Cv[:, np.newaxis, :])

    # Volumes to test
    nV = 7
    Vparlp_V = 10. * np.random.rand(nV)

    # ---------- Script ---------- #

    # Rescale cell
    cell_CVcv = cell_Ccv[:, np.newaxis, ...]\
        * (Vparlp_V**(1 / 3.))[np.newaxis, :, np.newaxis, np.newaxis]

    # Rescale primary vector to let Q.a follow Qrel
    Qrel_CQ = Qrel_Q[np.newaxis, :] / v0_C[:, np.newaxis]
    Qrel_CVQ = Qrel_CQ[:, np.newaxis, :]\
        / (Vparlp_V**(1 / 3.))[np.newaxis, :, np.newaxis]
    # Generate Q-vectors
    Q_CVdQv = Qrel_CVQ[..., np.newaxis, :, np.newaxis]\
        * v0n_Cv[:, np.newaxis, np.newaxis, np.newaxis, :]\
        + Q_Cdv[:, np.newaxis, :, np.newaxis, :]

    # Generate test values
    sinchalf_CVdQ = sinc(np.sum(cell_CVcv[..., np.newaxis, np.newaxis, 1, :]
                                * Q_CVdQv, axis=-1) / 2)\
        * sinc(np.sum(cell_CVcv[..., np.newaxis, np.newaxis, 2, :]
                      * Q_CVdQv, axis=-1) / 2)
    test_Theta_CVdQ = Vparlp_V[np.newaxis, :, np.newaxis, np.newaxis]\
        * sinchalf_Q[np.newaxis, np.newaxis, np.newaxis, :]\
        * sinchalf_CVdQ

    for Q_VdQv, test_Theta_VdQ, cell_Vcv in zip(Q_CVdQv, test_Theta_CVdQ,
                                                cell_CVcv):
        for Q_dQv, test_Theta_dQ, cell_cv in zip(Q_VdQv, test_Theta_VdQ,
                                                 cell_Vcv):
            for _ in range(3):  # Check that primary axis can be anywhere
                # Slide the cell axis indices
                cell_cv[:, :] = cell_cv[[2, 0, 1], :]
                # Calculate geometry factors
                Theta_dQ = parallelepipedic_geometry_factor(Q_dQv, cell_cv)

                # Check against expected results
                assert np.allclose(Theta_dQ, test_Theta_dQ, atol=1.e-8)


def Co_hcp_test():
    """Check that the site kernel interface works on run-time inputs.

    To do: Remember to check full diagonal as well as hermitian properties XXX
    """
    # ---------- Inputs ---------- #

    # Part 1: Generate plane wave representation (PWDescriptor)
    # Atomic configuration
    a = 2.5071
    c = 4.0695
    mm = 1.6
    # Ground state settings
    xc = 'LDA'
    kpts = 4
    pw = 200
    # Response settings
    ecut = 50.
    gammacentered = False
    q_c = [0., 0., 0.]

    # Part 2: Calculate site kernels
    # Define partitions to try
    rc_pa = np.array([[1., 2.], [2., 3.]])  # radii in Å
    # hc_pa = np.array([[2., 3.], [3., 4.]])  # heights in Å
    # ez_pav = np.array([[[0., 0., 1.], [1., 0., 0.]],
    #                    [[0., 0., 1.], [0., 0., 1.]]])
    # cell_acv = np.array([[[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
    #                      [[1., 0., 0.], [0., 1., 0.], [0., 0., 2.]]])
    # cell_pacv = np.append(cell_acv, cell_acv * 2).reshape(2, 2, 3, 3)

    # Part 3: Check the calculated kernels

    # ---------- Script ---------- #

    # Part 1: Generate plane wave representation (PWDescriptor)
    atoms = bulk('Co', 'hcp', a=a, c=c)
    atoms.set_initial_magnetic_moments([mm, mm])

    calc = GPAW(xc=xc,
                spinpol=True,
                mode=PW(pw),
                kpts={'size': (kpts, kpts, kpts),
                      'gamma': True}
                )

    pd0 = get_PWDescriptor(atoms, calc, q_c,
                           ecut=ecut,
                           gammacentered=gammacentered)

    # Part 2: Calculate site kernels
    positions = atoms.get_positions()

    # Generate spherical site kernels instances
    # Normally
    sph_sitekernels = SphericalSiteKernels(positions, rc_pa)
    # Separately as sum of site kernels
    sph_sitekernels0 = SphericalSiteKernels(positions[:1, :], rc_pa[:, :1])
    sph_sitekernels1 = SphericalSiteKernels(positions[1:, :], rc_pa[:, 1:])
    sph_sitekernels_sum = sph_sitekernels0 + sph_sitekernels1
    # Seperately as appended partitions
    sph_sitekernelsp0 = SphericalSiteKernels(positions, rc_pa[:1, :])
    sph_sitekernelsp1 = SphericalSiteKernels(positions, rc_pa[1:, :])
    sph_sitekernels_app = sph_sitekernelsp0.copy()
    sph_sitekernels_app.append(sph_sitekernelsp1)

    """
    # Generate cylindrical site kernels instances
    # Normally
    cyl_sitekernels = CylindricalSiteKernels(positions, [ez_mv],
                                             [rc_m], [height_m])
    # Separately as a sum of site kernels
    cyl_sitekernels0 = CylindricalSiteKernels([positions[0]], [[ez_mv[0]]],
                                              [[rc_m[0]]], [[height_m[0]]])
    cyl_sitekernels1 = CylindricalSiteKernels([positions[1]], [[ez_mv[1]]],
                                              [[rc_m[1]]], [[height_m[1]]])
    cyl_sitekernels_sum = cyl_sitekernels0 + cyl_sitekernels1

    # Generate unit cell site kernels instances
    # To get a single site kernel spanning the entire unit cell, we use
    # the unit cell center as the site position and the unit cell vectors
    # to define a parallelepipedic site kernel
    cell_cv = atoms.get_cell()
    cell_center_v = np.sum(cell_cv, axis=0) / 2.
    uc_sitekernels = ParallelepipedicSiteKernels([cell_center_v],
                                                 [[cell_cv]])
    # Do a sum with itself
    uc_sitekernels_sum = uc_sitekernels + uc_sitekernels

    # Do a sum of all unique site kernels to show, that we can compute them
    # in parallel
    all_sitekernels = sph_sitekernels + cyl_sitekernels + uc_sitekernels
    """

    # Calculate spherical site kernels
    Ksph_paGG = np.array([K_aGG for K_aGG in sph_sitekernels.calculate(pd0)])
    Ksph0_paGG = np.array([K_aGG for K_aGG in sph_sitekernels0.calculate(pd0)])
    Ksph1_paGG = np.array([K_aGG for K_aGG in sph_sitekernels1.calculate(pd0)])
    Ksph_sum_paGG = np.array([K_aGG
                              for K_aGG in sph_sitekernels_sum.calculate(pd0)])
    Ksphp0_paGG = np.array([K_aGG
                            for K_aGG in sph_sitekernelsp0.calculate(pd0)])
    Ksphp1_paGG = np.array([K_aGG
                            for K_aGG in sph_sitekernelsp1.calculate(pd0)])
    Ksph_app_paGG = np.array([K_aGG
                              for K_aGG in sph_sitekernels_app.calculate(pd0)])

    """
    # Calculate cylindrical site kernels
    Kcyl_mGG = list([K_mGG for K_mGG in cyl_sitekernels.calculate(pd0)])[0]
    Kcyl0_mGG = list([K_mGG for K_mGG in cyl_sitekernels0.calculate(pd0)])[0]
    Kcyl1_mGG = list([K_mGG for K_mGG in cyl_sitekernels1.calculate(pd0)])[0]
    Kcyl_sum_mGG = list([K_mGG for K_mGG in cyl_sitekernels_sum.calculate(pd0)])[0]

    # Calculate unit cell site kernels
    Kuc_mGG = list([K_mGG for K_mGG in uc_sitekernels.calculate(pd0)])[0]
    Kuc_sum_mGG = list([K_mGG for K_mGG in uc_sitekernels_sum.calculate(pd0)])[0]

    # Calculate all site kernels together
    Kall_mGG = list([K_mGG for K_mGG in all_sitekernels.calculate(pd0)])[0]
    """

    # Part 4: Check the calculated kernels

    # Check geometry shapes of basic arrays
    assert all([gs == 'sphere' for gs in sph_sitekernels.geometry_shapes])
    # assert cyl_sitekernels.geometry_shapes[0] == 'cylinder'
    # assert uc_sitekernels.geometry_shapes[0] == 'parallelepiped'

    # Check geometry shapes of summed arrays
    assert all([gs == 'sphere' for gs in sph_sitekernels_sum.geometry_shapes])
    # assert cyl_sitekernels_sum.geometry_shapes[0] == 'cylinder'
    # assert uc_sitekernels_sum.geometry_shapes[0] == 'parallelepiped'
    # assert all_sitekernels.geometry_shapes[0] is None

    # Check geometry shapes of appended arrays
    assert all([gs == 'sphere' for gs in sph_sitekernels_app.geometry_shapes])
    # To do XXX

    # Check shape of spherical kernel arrays
    nG = len(get_pw_coordinates(pd0))
    assert sph_sitekernels.shape == Ksph_paGG.shape[:2]
    assert Ksph_paGG.shape == rc_pa.shape + (nG, nG)
    assert Ksph0_paGG.shape == (rc_pa.shape[0], 1) + (nG, nG)
    assert Ksph0_paGG.shape == Ksph1_paGG.shape
    assert Ksph_sum_paGG.shape == rc_pa.shape + (nG, nG)
    assert Ksphp0_paGG.shape == (1, rc_pa.shape[1]) + (nG, nG)
    assert Ksphp0_paGG.shape == Ksphp1_paGG.shape
    assert Ksph_app_paGG.shape == rc_pa.shape + (nG, nG)

    """
    # Check shape of cylindrical kernel arrays
    assert Kcyl_mGG.shape == (2, nG, nG)
    assert Kcyl0_mGG.shape == (1, nG, nG)
    assert Kcyl1_mGG.shape == (1, nG, nG)
    assert Kcyl_sum_mGG.shape == (2, nG, nG)

    # Check shape of parallelepipedic kernel arrays
    assert Kuc_mGG.shape == (1, nG, nG)
    assert Kuc_sum_mGG.shape == (2, nG, nG)

    # Check shape of array calculated in parallel
    assert Kall_mGG.shape == (5, nG, nG)
    """

    # Check self-consitency of spherical arrays
    assert np.allclose(Ksph0_paGG[:, 0, ...], Ksph_sum_paGG[:, 0, ...])
    assert np.allclose(Ksph1_paGG[:, 0, ...], Ksph_sum_paGG[:, 1, ...])
    assert np.allclose(Ksph_paGG, Ksph_sum_paGG)
    assert np.allclose(Ksphp0_paGG[0], Ksph_app_paGG[0])
    assert np.allclose(Ksphp1_paGG[0], Ksph_app_paGG[1])
    assert np.allclose(Ksph_paGG, Ksph_app_paGG)

    """
    # Check self-consistency of cylindrical arrays
    assert np.allclose(Kcyl0_mGG[0], Kcyl_sum_mGG[0])
    assert np.allclose(Kcyl1_mGG[0], Kcyl_sum_mGG[1])
    assert np.allclose(Kcyl_mGG, Kcyl_sum_mGG)

    # Check self-consistency of unit cell arrays
    assert np.allclose(Kuc_mGG[0], Kuc_sum_mGG[0])
    assert np.allclose(Kuc_mGG[0], Kuc_sum_mGG[1])

    # Check self-consistency of kernel calculated in parallel
    assert np.allclose(Ksph_mGG, Kall_mGG[:2])
    assert np.allclose(Kcyl_mGG, Kall_mGG[2:4])
    assert np.allclose(Kuc_mGG, Kall_mGG[-1:])
    """

    # Check that K_00(q=0) gives Vint / V0 (fractional integration volume)
    # Volume of unit cell in Å^3
    V0 = atoms.get_volume()
    # Calculate integration volumes in Å^3
    Vsphere_pa = 4 / 3 * np.pi * rc_pa**3
    # Vcylinder_pa = np.pi * rc_pa**2 * hc_pa
    assert np.allclose(Ksph_paGG[..., 0, 0], Vsphere_pa / V0)
    # assert np.allclose(Kcyl_mGG[:, 0, 0], Vcylinder_m / V0)
    # assert abs(Kuc_paGG[0, 0, 0] - 1.) < 1.e-8


# ---------- Test functionality ---------- #


def get_PWDescriptor(atoms, calc, q_c, ecut=50., gammacentered=False):
    """Mock-up of PlaneWaveKSLRF.get_PWDescriptor working on a bare calculator
    instance without any actual data in it."""
    from ase.units import Ha
    from gpaw.pw.descriptor import PWDescriptor
    from gpaw.kpt_descriptor import KPointDescriptor

    # Perform inexpensive calculator initialization
    calc.initialize(atoms)

    # Create the plane wave descriptor
    q_c = np.asarray(q_c, dtype=float)
    qd = KPointDescriptor([q_c])
    pd = PWDescriptor(ecut / Ha, calc.wfs.gd,
                      complex, qd, gammacentered=gammacentered)

    return pd
