import numpy as np
from gpaw.response.mpa_interpolation import fit_residue, RESolver
from .mpa_interpolation_from_fortran import mpa_R_fit as fit_residue_fortran, mpa_RE_solver


def Xeval(Omega_GGp, residues_GGp, omega_w):
    X_GGpw = residues_GGp[...,:, np.newaxis] * 2 * Omega_GGp[...,:, np.newaxis] / (omega_w[None,None,None,:]**2 - Omega_GGp[...,:, np.newaxis]**2)
    return np.sum(X_GGpw, axis=2)

def test_residues():
    nG = 5
    npols = 10
    Omega_GGp = np.empty((nG,nG,npols), dtype=np.complex128)
    residues_GGp = np.empty((nG,nG,npols), dtype=np.complex128)
    X_GGw = np.empty((nG,nG,2*npols), dtype=np.complex128)
    R_GGp = np.empty((nG,nG,npols), dtype=np.complex128)
    R_fortran_GGp = np.empty((nG,nG,npols), dtype=np.complex128)
    omega_w = np.linspace(0., 5., 2*npols) + 0.1j

    rng = np.random.default_rng(seed=1)
    for g1 in range(nG):
        for g2 in range(nG):
            Omega_GGp[g1,g2] = rng.random(npols)*0.05 + 5.5 - 0.01j
            residues_GGp[g1,g2] = rng.random(npols)
            X_GGw[g1,g2] = Xeval(Omega_GGp[g1,g2], residues_GGp[g1,g2], omega_w)
            R_fortran_GGp[g1,g2] = fit_residue_fortran(npols, npols, omega_w, X_GGw[g1,g2], Omega_GGp[g1,g2])

    R_pGG = fit_residue(np.ones((nG,nG))*npols, omega_w, X_GGw.transpose(2,0,1), Omega_GGp.transpose(2,0,1))

    R_GGp = R_pGG.transpose(1,2,0)

    X_fit_GGw = Xeval(Omega_GGp, R_GGp,omega_w)
    X_fortran_fit_GGw = Xeval(Omega_GGp, R_fortran_GGp,omega_w)
    assert np.allclose(X_fit_GGw, X_fortran_fit_GGw, atol=1e-6)


    if 0:
        g1, g2 = 0,0
        from matplotlib import pyplot as plt
        plt.plot(omega_w, X_GGw[g1,g2].real, 'k',ls='--')
        plt.plot(omega_w, X_GGw[g1,g2].imag, 'gray',ls='--')

        plt.plot(omega_w, X_fit_GGw[g1,g2].real)
        plt.plot(omega_w, X_fit_GGw[g1,g2].imag)
        plt.plot(omega_w, X_fortran_fit_GGw[g1,g2].real,ls=':')
        plt.plot(omega_w, X_fortran_fit_GGw[g1,g2].imag,ls=':')
        plt.show()



def test_poles():
    nG = 2
    npols = 100
    Omega_GGp = np.empty((nG,nG,npols), dtype=np.complex128)
    residues_GGp = np.empty((nG,nG,npols), dtype=np.complex128)
    wmax = 2 

    npols_mpa = 6
    omega_p = np.linspace(0,wmax, npols_mpa)
    omega_w = np.concatenate((omega_p + 0.1j, omega_p +1.j))

    X_GGw = np.empty((nG,nG,2*npols_mpa), dtype=np.complex128)
    E_GGp = np.empty((nG,nG,npols_mpa), dtype=np.complex128)
    R_GGp = np.empty((nG,nG,npols_mpa), dtype=np.complex128)
    E_fortran_GGp = np.empty((nG,nG,npols_mpa), dtype=np.complex128)
    R_fortran_GGp = np.empty((nG,nG,npols_mpa), dtype=np.complex128)

    for g1 in range(nG):
        for g2 in range(nG):
            Omega_GGp[g1,g2] = np.random.normal(1, 0.5, npols)  - 0.05j
            residues_GGp[g1,g2] = np.random.rand(npols)
            X_GGw[g1,g2] = Xeval(Omega_GGp[g1,g2], residues_GGp[g1,g2], omega_w)

            R_fortran_GGp[g1,g2], E_fortran_GGp[g1,g2], _, _ = mpa_RE_solver(npols_mpa, omega_w, X_GGw[g1,g2])

    E_pGG, R_pGG = RESolver(omega_w).solve(X_GGw.transpose(2,0,1))

    E_GGp = E_pGG.transpose(1,2,0)
    R_GGp = R_pGG.transpose(1,2,0)

    omega_grid = np.linspace(0., wmax, 100)  + 0.01j
    X_fit_GGw = Xeval(E_GGp, R_GGp,omega_grid)
    X_fortran_fit_GGw = Xeval(E_fortran_GGp, R_fortran_GGp,omega_grid)
    assert np.allclose(X_fit_GGw, X_fortran_fit_GGw, atol=1e-6)

    if 1:
        X_num_GGw = Xeval(Omega_GGp, residues_GGp, omega_grid)
        g1, g2 = 0,0
        from matplotlib import pyplot as plt
        plt.plot(omega_grid.real, X_num_GGw[g1,g2].real, 'k',ls='--')
        plt.plot(omega_grid.real, X_num_GGw[g1,g2].imag, 'gray',ls='--')

        plt.plot(omega_grid.real, X_fit_GGw[g1,g2].real)
        plt.plot(omega_grid.real, X_fit_GGw[g1,g2].imag)
        plt.plot(omega_grid.real, X_fortran_fit_GGw[g1,g2].real,ls=':')
        plt.plot(omega_grid.real, X_fortran_fit_GGw[g1,g2].imag,ls=':')
        plt.show()
