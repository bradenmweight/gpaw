from gpaw.utilities.extrapolate import extrapolate
import numpy as np

E_pbe, E_hf = np.genfromtxt('diamond.ralda.PBE_HF_diamond.dat')[:, 1]
assert abs(E_pbe - -7.75) < 0.01
assert abs(E_hf - -5.17) < 0.01

a = np.loadtxt('diamond.ralda.rpa.dat')
b = np.loadtxt('CO.ralda_rpa_C.dat')
ext, A, B, sigma = extrapolate(a[:, 0], a[:, 1] / 2 - b[:, 1],
                               reg=3, plot=False)
assert abs(A - -1.89) < 0.01, A
assert abs(E_hf + A - -7.06) < 0.01

if 0:
    # RAPBE currently not working!
    a = np.loadtxt('diamond.ralda.rapbe.dat')
    b = np.loadtxt('CO.ralda_rapbe_C.dat')
    ext, A, B, sigma = extrapolate(a[:, 0], a[:, 1] / 2 - b[:, 1],
                                   reg=3, plot=False)
    assert abs(A - -1.48) < 0.01
    assert abs(E_hf + A - -6.65) < 0.01
