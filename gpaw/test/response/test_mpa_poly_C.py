import numpy as np
from ase.units import Hartree as Ha
from _gpaw import evaluate_mpa_poly as mpa_C


def mpa_py(omega, f, omegat_nGG, W_nGG, eta, factor):
    x1_nGG = f / (omega + omegat_nGG - 1j * eta)
    x2_nGG = (1.0 - f) / (omega - omegat_nGG + 1j * eta)

    x_GG = 2 * factor * np.sum(W_nGG * (x1_nGG + x2_nGG),
                               axis=0)  # Why 2 here

    eps = 0.0001 / Ha
    xp_nGG = f / (omega + eps + omegat_nGG - 1j * eta)
    xp_nGG += (1.0 - f) / (omega + eps - omegat_nGG + 1j * eta)
    xm_nGG = f / (omega - eps + omegat_nGG - 1j * eta)
    xm_nGG += (1.0 - f) / (omega - eps - omegat_nGG + 1j * eta)
    dx_GG = 2 * factor * np.sum(W_nGG * (xp_nGG - xm_nGG) / (2 * eps),
                                axis=0)  # Why 2 here
    return x_GG, dx_GG


def test_residues(in_tmp_dir):
    f = 0.5
    factor = 2.0
    eta = 0.1 * Ha
    nG = 5
    npols = 10
    omegat_nGG = np.empty((npols, nG, nG), dtype=np.complex128)
    W_nGG = np.empty((npols, nG, nG), dtype=np.complex128)
    omega = 0.5

    rng = np.random.default_rng(seed=1)
    omegat_nGG = rng.random((npols, nG, nG)) * 0.05 + 5.5 - 0.01j
    W_nGG[:] = rng.random((npols, nG, nG))
    W_nGG = np.ascontiguousarray(W_nGG)

    x_GG_py, dx_GG_py = mpa_py(omega, f, omegat_nGG, W_nGG, eta, factor)

    x_GG_C = np.empty(omegat_nGG.shape[1:], dtype=complex)
    dx_GG_C = np.empty(omegat_nGG.shape[1:], dtype=complex)
    mpa_C(x_GG_C, dx_GG_C, omega, f, omegat_nGG, W_nGG, eta, factor)

    print(x_GG_py)
    print(x_GG_C)
    print(x_GG_py / x_GG_C)

    assert np.allclose(x_GG_py, x_GG_C, atol=1e-6)
