import numpy as np
from gpaw.utilities.extrapolate import extrapolate

E_pbe, E_hf = np.genfromtxt('CO.ralda.PBE_HF_CO.dat')[:, 1]
assert abs(E_pbe - -11.74) < 0.01
assert abs(E_hf - -7.37) < 0.01

CO_rpa = np.loadtxt('CO.ralda_rpa_CO.dat')
C_rpa = np.loadtxt('CO.ralda_rpa_C.dat')
O_rpa = np.loadtxt('CO.ralda_rpa_O.dat')

a = CO_rpa
a[:, 1] -= (C_rpa[:, 1] + O_rpa[:, 1])
ext, A, B, sigma = extrapolate(a[:, 0], a[:, 1], reg=3, plot=False)
assert abs(A - -3.27) < 0.02, A
assert abs(E_hf + A - -10.64) < 0.01

if 0:
    # RAPBE not currently working!
    CO_rapbe = np.loadtxt('CO.ralda_rapbe_CO.dat')
    C_rapbe = np.loadtxt('CO.ralda_rapbe_C.dat')
    O_rapbe = np.loadtxt('CO.ralda_rapbe_O.dat')

    a = CO_rapbe
    a[:, 1] -= (C_rapbe[:, 1] + O_rapbe[:, 1])
    ext, A, B, sigma = extrapolate(a[:, 0], a[:, 1], reg=3, plot=False)
    assert abs(A - -3.51) < 0.01
    assert abs(E_hf + A - -10.88) < 0.01
