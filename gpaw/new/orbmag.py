r"""This module calculates the orbital magnetic moment vector for each atom.

The orbital magnetic moment is calculated in the atom-centred approximation
where only the PAW correction to the wave function is assumed to contribute.
This leads to the equation (presented in SI units):

::

                 ===  ===
   a         e   \    \         / a   \*  a      a
  m      = - --  /    /    f   | P     | P      L
   orb,v     2m  ===  ===   kn  \ knsi/   knsi'  vii'
                 kn   sii'

with L^a_vii' containing the matrix elements of the angular momentum operator
between two partial waves centred at atom a.

The orbital magnetic moments are returned in units of μ_B without the sign of
the negative electronic charge, q = - e.
"""

import numpy as np
from gpaw.new import zips
from gpaw.spinorbit import get_L_vlmm

L_vlmm = get_L_vlmm()


def get_orbmag_from_density(D_asii, n_aj, l_aj):
    """Returns orbital magnetic moment vectors for each atom a
    calculated from its respective atomic density matrix."""

    orbmag_av = np.zeros([len(n_aj), 3])
    for (a, D_sii), n_j, l_j in zips(D_asii.items(), n_aj, l_aj):
        D_ii = D_sii[0]  # Only the electron density

        Ni = 0
        for n, l in zips(n_j, l_j):
            Nm = 2 * l + 1
            if n < 0:
                Ni += Nm
                continue
            for v in range(3):
                orbmag_av[a, v] += np.einsum('ij,ij->',
                                             D_ii[Ni:Ni + Nm, Ni:Ni + Nm],
                                             L_vlmm[v][l]).real
            Ni += Nm

    return orbmag_av


# Unused, but it might be good to keep around for the purposes of analysis.
def calculate_orbmag_1k(f_n, P_nsi, nl_j):
    """Calculate contribution to orbital magnetic moment for a single k-point.

    Parameters
    ----------
    f_n : list or ndarray
        Occupations for each state n
        (Fermi-Dirac occupation multiplied by k-point weight)
    P_nsi : ndarray
        Projector overlaps for each state n, spin s, and partial wave i
    nl_j : sequence of tuples
        Principal quantum number and angular momentum quantum number
        for each radial function j

    NB: i is an index for all partial waves for one atom and j is an index for
    only the radial wave function which is used to build all of the partial
    waves. i and j do not refer to the same kind of index.

    Only pairs of partial waves with the same radial function may yield
    nonzero contributions. The sum can therefore be limited to diagonal blocks
    of shape [2 * l_j + 1, 2 * l_j +1] where l_j is the angular momentum
    quantum number of the j'th radial function.

    Partials with unbounded radial functions (negative n_j) are skipped.
    """

    orbmag_v = np.zeros(3)
    Ni = 0
    for n, l in nl_j:
        Nm = 2 * l + 1
        if n < 0:
            Ni += Nm
            continue
        for v in range(3):
            orbmag_v[v] += np.einsum('nsi,nsj,n,ij->',
                                     P_nsi[:, :, Ni:Ni + Nm].conj(),
                                     P_nsi[:, :, Ni:Ni + Nm],
                                     f_n, L_vlmm[v][l]).real
        Ni += Nm

    return orbmag_v
