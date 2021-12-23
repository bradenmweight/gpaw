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

    # Get number of reciprocal lattice vectors
    G_Gc = get_pw_coordinates(pd)
    NG = len(G_Gc)

    # Number of sites
    if shapes_m == 'unit cell':
        N_sites = 1
    else:
        N_sites = len(sitePos_mv)

    # Reformat shape parameters
    if type(shapes_m) is str:
        shapes_m = np.array([shapes_m]*N_sites)
    if type(rc_m) in {int, float}:
        rc_m = np.array([rc_m]*N_sites)
    if type(zc_m) in {int, float, str}:
        zc_m = np.array([zc_m]*N_sites)

    # Array to fill
    K_GGm = np.zeros([NG, NG, N_sites], dtype=np.complex128)

    # --- The Calculation itself --- #

    # Loop through magnetic sites
    for m in range(N_sites):
        # Get site specific values
        shape, rc, zc, sitePos_v = shapes_m[m], rc_m[m], zc_m[m], \
                                   sitePos_mv[m, :]

        # Do computation for relevant shape
        if shape == 'sphere':
            K_GG = K_sphere(pd, sitePos_v=sitePos_v, rc=rc)

        elif shape == 'cylinder':
            K_GG = K_cylinder(pd, sitePos_v=sitePos_v, rc=rc, zc=zc)

        elif shape == 'unit cell':
            K_GG = K_unit_cell(pd, sitePos_v=sitePos_v)

        else:
            print('Not a recognised shape')

        # Update data
        K_GGm[:, :, m] = K_GG

    return K_GGm


def K_sphere(pd, sitePos_v, rc=1.0):
    """Compute site-kernel for a spherical integration region """

    # Get relevant quantities from pd object
    G_Gv, q_v, Omega_cell = _extract_pd_info(pd)
    NG = len(G_Gv)

    # Convert from Å to Bohr
    rc = rc / Bohr
    sitePos_v = sitePos_v / Bohr

    # Construct arrays
    G1_GGv, G2_GGv, q_GGv, sum_GGv = _constructArrays(G_Gv, q_v)

    # Combine arrays
    magsq_GG = np.sum(sum_GGv ** 2, axis=-1)  # |G_1 + G_2 + q|^2
    mag_GG = np.sqrt(magsq_GG)  # |G_1 + G_2 + q|

    # Find singular and regular points
    is_sing = mag_GG * rc < 1.e-8  # cutoff is arbitrary
    is_reg = np.logical_not(is_sing)

    # Separate calculation into regular and singular part
    magReg_GG = mag_GG[is_reg]
    magSing_GG = mag_GG[is_sing]
    magsqReg_GG = magsq_GG[is_reg]

    # Compute integral part of kernel
    K_GG = np.zeros([NG, NG], dtype=np.complex128)
    # Full formula
    K_GG[is_reg] = 4*np.pi / magsqReg_GG * \
        (-rc*np.cos(magReg_GG*rc) + np.sin(magReg_GG*rc)/magReg_GG)
    # Taylor expansion around singularity
    K_GG[is_sing] = 4*np.pi*rc**3/3 - 2*np.pi/15 * magSing_GG**2 * rc**5

    # Compute complex prefactor
    prefactor = _makePrefactor(sitePos_v, sum_GGv, Omega_cell)
    K_GG *= prefactor
    
    return K_GG


def K_cylinder(pd, sitePos_v, rc=1.0, zc='unit cell'):
    """Compute site-kernel for a cylindrical integration region"""

    # Get reciprocal lattice vectors and q-vector from pd
    q_qc = pd.kd.bzk_kc
    assert len(q_qc) == 1
    q_c = q_qc[0, :]     # Assume single q
    G_Gc = get_pw_coordinates(pd)
    NG = len(G_Gc)

    # Convert to cartesian coordinates
    B_cv = 2.0 * np.pi * pd.gd.icell_cv  # Coordinate transform matrix
    q_v = np.dot(q_c, B_cv)  # Unit = Bohr^(-1)
    G_Gv = np.dot(G_Gc, B_cv)

    # Set height to that of unit cell (only makes sense in 2D)
    if zc == 'unit cell':
        zc = np.sum(pd.gd.cell_cv[:, -1]) * Bohr   # Units of Å.
    elif zc == 'diameter':
        zc = 2*rc

    # Get unit cell volume in bohr^3
    Omega_cell = pd.gd.volume

    # Convert from Å to Bohr
    rc = rc / Bohr
    zc = zc / Bohr
    sitePos_v = sitePos_v / Bohr

    # Construct arrays
    G1_GGv = np.tile(G_Gv[:, np.newaxis, :], [1, NG, 1])
    G2_GGv = np.tile(G_Gv[np.newaxis, :, :], [NG, 1, 1])
    q_GGv = np.tile(q_v[np.newaxis, np.newaxis, :], [NG, NG, 1])

    # Combine arrays
    sum_GGv = G1_GGv + G2_GGv + q_GGv  # G_1 + G_2 + q
    # sqrt([G1_x + G2_x + q_x]^2 + [G1_y + G2_y + q_y]^2)
    Qrho_GG = np.sqrt(sum_GGv[:, :, 0]**2 + sum_GGv[:, :, 1]**2)
    Qz_GG = sum_GGv[:, :, 2]  # G1_z + G2_z + q_z

    # Set values of |G_1 + G_2 + q|*r_c below sing_cutoff equal to
    #   sing_cutoff (deals with division by 0)
    # Note : np.sinc does this on it's own, so Qz_GGq needs no adjustment
    sing_cutoff = 1.0e-15
    Qrho_GG = np.where(np.abs(Qrho_GG)*rc < sing_cutoff,
                       sing_cutoff/rc, Qrho_GG)

    # Compute site kernel
    K_GG = 2*np.pi*zc*rc**2 * sinc(Qz_GG*zc/2) * jv(1, rc*Qrho_GG)/(rc*Qrho_GG)

    # Compute complex prefactor
    tau_GGv = np.tile(sitePos_v[np.newaxis, np.newaxis, :], [NG, NG, 1])
    # e^{i tau_mu . (G_1 + G_2 + q)}
    phaseFactor_GG = np.exp(1j * np.sum(tau_GGv * sum_GGv, axis=-1))
    K_GG = K_GG.astype(np.complex128)
    K_GG *= phaseFactor_GG  # Phase factor
    K_GG *= np.sqrt(2) / Omega_cell ** (3 / 2)  # Real-valued prefactor

    return K_GG


def K_unit_cell(pd, sitePos_v=None):
    """Compute site-kernel for a spherical integration region"""

    # Get reciprocal lattice vectors and q-vector from pd
    q_qc = pd.kd.bzk_kc
    assert len(q_qc) == 1
    q_c = q_qc[0, :]     # Assume single q
    G_Gc = get_pw_coordinates(pd)
    NG = len(G_Gc)

    # Convert to cartesian coordinates
    B_cv = 2.0 * np.pi * pd.gd.icell_cv  # Coordinate transform matrix
    q_v = np.dot(q_c, B_cv)  # Unit = Bohr^(-1)
    G_Gv = np.dot(G_Gc, B_cv)

    # Get unit cell vectors and volume in Bohr
    Omega_cell = pd.gd.volume
    a1, a2, a3 = pd.gd.cell_cv

    # Default is center of unit cell
    if sitePos_v is None:
        sitePos_v = 1/2 * (a1 + a2 + a3)

    # Construct arrays
    G1_GGv = np.tile(G_Gv[:, np.newaxis, :], [1, NG, 1])
    G2_GGv = np.tile(G_Gv[np.newaxis, :, :], [NG, 1, 1])
    q_GGv = np.tile(q_v[np.newaxis, np.newaxis, :], [NG, NG, 1])

    # Combine arrays
    sum_GGv = G1_GGv + G2_GGv + q_GGv  # G_1 + G_2 + q

    # Compute site-kernel
    K_GG = sinc(sum_GGv@a1 / 2) * sinc(sum_GGv@a2 / 2) * sinc(sum_GGv@a3 / 2)

    # Multiply by prefactor
    # e^{i (G_1 + G_2 + q) . (a_1 + a_2 + a_3)/2}
    phase_factor_GG = np.exp(1j * sum_GGv @ sitePos_v)
    K_GG = K_GG * np.sqrt(2/Omega_cell) * phase_factor_GG

    return K_GG


def _makePrefactor(sitePos_v, sum_GGv, Omega_cell):
    """Make the complex prefactor which occurs for all site-kernels,
    irrespective of shape of integration region"""
    # Phase factor
    phaseFactor_GG = np.exp(1j * sum_GGv @ sitePos_v)

    # Scale factor
    scaleFactor = np.sqrt(2) / Omega_cell ** (3 / 2)

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


def _constructArrays(G_Gv, q_v):
    """Construct arrays with shape (NG, NG, 3)"""
    NG = len(G_Gv)
    G1_GGv = np.tile(G_Gv[:, np.newaxis, :], [1, NG, 1])
    G2_GGv = np.tile(G_Gv[np.newaxis, :, :], [NG, 1, 1])
    q_GGv = np.tile(q_v[np.newaxis, np.newaxis, :], [NG, NG, 1])
    sum_GGv = G1_GGv + G2_GGv + q_GGv  # G_1 + G_2 + q

    return G1_GGv, G2_GGv, q_GGv, sum_GGv
