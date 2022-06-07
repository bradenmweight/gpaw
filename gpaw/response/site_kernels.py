"""Compute site-kernels. Used for computing Heisenberg exchange.
Specifically, one maps DFT calculations onto a Heisenberg lattice model,
where the site-kernels define the lattice sites and magnetic moments."""

import numpy as np
from scipy.special import jv
from gpaw.response.susceptibility import get_pw_coordinates
from ase.units import Bohr


def sinc(x):
    """np.sinc(x) = sin(pi*x) / (pi*x), hence the division by pi"""
    return np.sinc(x / np.pi)


def site_kernel_interface(pd, sitePos_mv, shapes_m='sphere',
                          rc_m=1.0, zc_m='diameter'):
    """Compute site kernels using an arbitrary combination of shapes for the
        integration region at different magnetic sites.
        Accepts spheres, cylinders and unit cell.

    Parameters
    ----------
    pd : Planewave Descriptor. Contains mixed information about
            plane-wave basis
    sitePos_mv : Array with positions of magnetic sites within unit cell
    shapes_m : List of str or str. Which shapes to use for the
            different integration regions.
            Options are 'sphere', 'cylinder', 'unit cell'.
            If 'unit cell' then rc_rm, zc_rm, sitePos_mv do nothing and
                a single site in the middle of the unit cell is used.
    rc_m : Radius of integration region
    zc_m : Height of integration cylinders (if any)
            If 'diameter' then zc=2*rc
            If 'unit cell' then use height of unit cell (makes sense in 2D)
    """

    # Number of sites
    if shapes_m == 'unit cell':
        nsites = 1
    else:
        nsites = len(sitePos_mv)

    # Reformat shape parameters
    if type(shapes_m) is str:
        shapes_m = np.array([shapes_m] * nsites)
    if type(rc_m) in [int, float]:
        rc_m = np.array([rc_m] * nsites)
    if type(zc_m) in [int, float]:
        zc_m = np.array([zc_m] * nsites)
    if type(zc_m) == str:
        zc_m = [zc_m] * nsites
    for m in range(nsites):
        # Reformat zc if needed
        if zc_m[m] == 'unit cell':
            zc_m[m] = np.sum(pd.gd.cell_cv[:, -1]) * Bohr   # Units of Å.
        elif zc_m[m] == 'diameter':
            zc_m[m] = 2. * rc_m[m]
    zc_m = np.array(zc_m, dtype=float)

    # Convert input units (Å) to atomic units (Bohr)
    sitePos_mv = sitePos_mv / Bohr
    rc_m = rc_m / Bohr
    zc_m = zc_m / Bohr

    # Construct Fourier components
    G_Gv, q_v, Omega_cell = _extract_pd_info(pd)
    Q_GGv = _construct_wave_vectors(G_Gv, q_v)

    # Array to fill
    K_GGm = np.zeros(Q_GGv.shape[:2] + (nsites,), dtype=np.complex128)

    # --- The Calculation itself --- #

    # Loop through magnetic sites
    # Should be vectorized? XXX
    for m in range(nsites):
        # Get site specific values
        shape, rc, zc = shapes_m[m], rc_m[m], zc_m[m]
        sitePos_v = sitePos_mv[m, :]

        # Compute complex prefactor
        prefactor = _makePrefactor(sitePos_v, Q_GGv, Omega_cell)

        # Do computation for relevant shape
        if shape == 'sphere':
            K_GG = K_sphere(Q_GGv, rc=rc)

        elif shape == 'cylinder':
            K_GG = K_cylinder(Q_GGv, rc=rc, zc=zc)

        elif shape == 'unit cell':
            # Get real-space basis vectors
            # Give the user control over these XXX
            a1, a2, a3 = pd.gd.cell_cv

            # # Default site position is center of unit cell
            # # This should not be up to some secret functionality to decide XXX
            # if sitePos_v is None:
            #     sitePos_v = 1 / 2 * (a1 + a2 + a3)
            K_GG = K_unit_cell(Q_GGv, a1, a2, a3)

        else:
            print('Not a recognised shape')

        # Update data
        K_GGm[:, :, m] = prefactor * K_GG

    return K_GGm


def K_sphere(Q_GGv, rc=1.0):
    """Compute site-kernel for a spherical integration region """

    # Combine arrays
    magsq_GG = np.sum(Q_GGv ** 2, axis=-1)  # |G_1 + G_2 + q|^2
    mag_GG = np.sqrt(magsq_GG)  # |G_1 + G_2 + q|

    # Find singular and regular points
    is_sing = mag_GG * rc < 1.e-8  # cutoff is arbitrary
    is_reg = np.logical_not(is_sing)

    # Separate calculation into regular and singular part
    magReg_GG = mag_GG[is_reg]
    magSing_GG = mag_GG[is_sing]
    magsqReg_GG = magsq_GG[is_reg]

    # Compute integral part of kernel
    K_GG = np.zeros(Q_GGv.shape[:2], dtype=np.complex128)
    # Full formula
    K_GG[is_reg] = 4 * np.pi / magsqReg_GG * \
        (-rc * np.cos(magReg_GG * rc) + np.sin(magReg_GG * rc) / magReg_GG)
    # Taylor expansion around singularity
    K_GG[is_sing] = 4 * np.pi * rc**3 / 3\
        - 2 * np.pi / 15 * magSing_GG**2 * rc**5
    
    return K_GG


def K_cylinder(Q_GGv, rc=1.0, zc=1.0):
    """Compute site-kernel for a cylindrical integration region"""

    # Combine arrays
    # sqrt([G1_x + G2_x + q_x]^2 + [G1_y + G2_y + q_y]^2)
    Qrho_GG = np.sqrt(Q_GGv[:, :, 0]**2 + Q_GGv[:, :, 1]**2)
    Qz_GG = Q_GGv[:, :, 2]  # G1_z + G2_z + q_z

    # Set values of |G_1 + G_2 + q|*r_c below sing_cutoff equal to
    #   sing_cutoff (deals with division by 0)
    # Note : np.sinc does this on it's own, so Qz_GGq needs no adjustment
    sing_cutoff = 1.0e-15
    Qrho_GG = np.where(np.abs(Qrho_GG) * rc < sing_cutoff,
                       sing_cutoff / rc, Qrho_GG)

    # Compute site kernel
    K_GG = 2 * np.pi * zc * rc**2 * sinc(Qz_GG * zc / 2)\
        * jv(1, rc * Qrho_GG) / (rc * Qrho_GG)

    return K_GG


def K_unit_cell(Q_GGv, a1, a2, a3):
    """Compute site-kernel for a spherical integration region"""

    # Calculate the paralleliped volume
    cell_cv = np.array([a1, a2, a3])
    Vparlp = abs(np.linalg.det(cell_cv))

    # Calculate the site-kernel
    K_GG = Vparlp * sinc(Q_GGv @ a1 / 2) * sinc(Q_GGv @ a2 / 2) * \
        sinc(Q_GGv @ a3 / 2)

    return K_GG


def _makePrefactor(sitePos_v, sum_GGv, Omega_cell):
    """Make the complex prefactor which occurs for all site-kernels,
    irrespective of shape of integration region"""
    # Phase factor
    phaseFactor_GG = np.exp(1j * sum_GGv @ sitePos_v)

    # Scale factor
    scaleFactor = 1. / Omega_cell

    return scaleFactor * phaseFactor_GG


def _extract_pd_info(pd):
    """Get relevant quantities from pd object (plane-wave descriptor)
    In particular reciprocal space vectors and unit cell volume
    Note : all in Bohr and absolute coordinates."""
    q_qc = pd.kd.bzk_kc
    assert len(q_qc) == 1
    q_c = q_qc[0, :]  # Assume single q
    G_Gc = get_pw_coordinates(pd)

    # Convert to cartesian coordinates
    B_cv = 2.0 * np.pi * pd.gd.icell_cv  # Coordinate transform matrix
    q_v = np.dot(q_c, B_cv)  # Unit = Bohr^(-1)
    G_Gv = np.dot(G_Gc, B_cv)

    # Get unit cell volume in bohr^3
    Omega_cell = pd.gd.volume

    return G_Gv, q_v, Omega_cell


def _construct_wave_vectors(G_Gv, q_v):
    """Construct wave vector array with shape (NG, NG, 3)"""
    NG = len(G_Gv)
    G1_GGv = np.tile(G_Gv[:, np.newaxis, :], [1, NG, 1])
    G2_GGv = np.tile(G_Gv[np.newaxis, :, :], [NG, 1, 1])
    q_GGv = np.tile(q_v[np.newaxis, np.newaxis, :], [NG, NG, 1])

    Q_GGv = G1_GGv + G2_GGv + q_GGv  # G_1 + G_2 + q

    return Q_GGv
