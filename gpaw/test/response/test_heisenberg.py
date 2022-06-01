"""Test the Heisenberg model based methodology of the response code."""

# General modules
import numpy as np

# Script modules
from gpaw.response.heisenberg import calculate_magnon_energy_single_site,\
    compute_magnon_energy_FM


# ---------- Main test functionality ---------- #


def test_heisenberg():
    magnon_dispersion_tests()


def magnon_dispersion_tests():
    single_site_magnons_test()
    single_site_magnons_consistency_test()
    multiple_sites_magnons_test()


def single_site_magnons_test():
    """Check the single site magnon dispersion functionality."""
    # ---------- Inputs ---------- #

    # Magnetic moment
    mm = 1.
    # q-point grid
    q_qc = np.zeros((11, 3), dtype=np.float)
    q_qc[:, 2] = np.linspace(0., np.pi, 11)
    np.random.shuffle(q_qc[:, 2])

    # Random J_q, with J=0 at q=0
    J_q = np.random.rand(q_qc.shape[0])
    J_q[list(q_qc[:, 2]).index(0.)] = 0.

    # Cosine J_qD with different spin wave stiffnesses D
    D_D = np.linspace(400., 800., 5)
    J_qD = D_D[np.newaxis, :] * np.cos(q_qc[:, 2])[:, np.newaxis]

    # ---------- Script ---------- #

    # Calculate magnon energies
    E_q = calculate_magnon_energy_single_site(J_q, q_qc, mm)
    E_qD = calculate_magnon_energy_single_site(J_qD, q_qc, mm)

    # Check dimensions of arrays
    assert E_q.shape == (q_qc.shape[0],)
    assert E_qD.shape == J_qD.shape

    # Check versus formulas
    assert np.allclose(E_q, -2. / mm * J_q)  # Remember: J(0) = 0
    assert np.allclose(E_qD, 2. / mm * D_D[np.newaxis, :]
                       * (1. - np.cos(q_qc[:, 2]))[:, np.newaxis])


def single_site_magnons_consistency_test():
    """Check that the generalized magnon dispersion calculation is consistent
    for a single site system with the simple analytical formula valid in that
    case."""
    # ---------- Inputs ---------- #

    # Magnetic moment
    mm = 1.
    # q-point grid
    q_qc = np.zeros((11, 3), dtype=np.float)
    q_qc[:, 2] = np.linspace(0., np.pi, 11)
    np.random.shuffle(q_qc[:, 2])

    # Random isotropic exchange constants
    J_q = np.random.rand(q_qc.shape[0])

    # ---------- Script ---------- #

    # Calculate assuming a single site
    E_q = calculate_magnon_energy_single_site(J_q, q_qc, mm)

    # Calcualte using generalized functionality
    E_nq = compute_magnon_energy_FM(J_q[np.newaxis, np.newaxis, :], q_qc, mm)

    # Test self-consistency
    assert E_nq.shape[0] == 1
    assert E_nq.shape[1] == len(E_q)
    assert np.allclose(E_nq[0, :], E_q, atol=1e-8)


def multiple_sites_magnons_test():
    pass
