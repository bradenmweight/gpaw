import os

import numpy as np

from ase import Atoms
from gpaw import GPAW, PW
from gpaw.mpi import world
from gpaw.test import equal, findpeak
from gpaw.utilities import compiled_with_sl
from gpaw.response.df import DielectricFunction
from ase.build import bulk

# Comparing the plasmon peaks found in bulk sodium for two different
# atomic structures. Testing for idential plasmon peaks. Not using
# physical sodium cell.

a1 = bulk('Na')
a2 = bulk('Na')
a2.set_initial_magnetic_moments([[0.1, ]])
a1.calc = GPAW(mode=PW(300),
               kpts={'size': (8, 8, 8), 'gamma': True},
               parallel={'band': 1},
               txt='na_spinpaired.txt')
a2.calc = GPAW(mode=PW(300),
               kpts={'size': (8, 8, 8), 'gamma': True},
               parallel={'band': 1},
               txt='na_spinpol.txt')
a1.get_potential_energy()
a2.get_potential_energy()

# Use twice as many bands for expanded structure
a1.calc.diagonalize_full_hamiltonian(nbands=20)
a2.calc.diagonalize_full_hamiltonian(nbands=20)

a1.calc.write('intraband_spinpaired.gpw', 'all')
a2.calc.write('intraband_spinpolarized.gpw', 'all')

# Calculate the dielectric functions
try:
    os.remove('intraband_spinpaired+0+0+0.pckl')
except OSError:
    pass
    
df1 = DielectricFunction('intraband_spinpaired.gpw',
                         domega0=0.03,
                         ecut=10,
                         rate=0.1,
                         integrationmode='tetrahedron integration',
                         name='intraband_spinpaired')

df1NLFCx, df1LFCx = df1.get_dielectric_function(direction='x')
df1NLFCy, df1LFCy = df1.get_dielectric_function(direction='y')
df1NLFCz, df1LFCz = df1.get_dielectric_function(direction='z')
wp1_vv = df1.chi0.plasmafreq_vv**0.5
wp1 = wp1_vv[0, 0]
try:
    os.remove('intraband_spinpolarized+0+0+0.pckl')
except OSError:
    pass
    
df2 = DielectricFunction('intraband_spinpolarized.gpw',
                         domega0=0.03,
                         ecut=10,
                         rate=0.1,
                         integrationmode='tetrahedron integration',
                         name='intraband_spinpolarized')

df2NLFCx, df2LFCx = df2.get_dielectric_function(direction='x')
df2NLFCy, df2LFCy = df2.get_dielectric_function(direction='y')
df2NLFCz, df2LFCz = df2.get_dielectric_function(direction='z')
wp2_vv = df2.chi0.plasmafreq_vv**0.5
wp2 = wp2_vv[0, 0]

# Compare plasmon frequencies and intensities
w_w = df1.chi0.omega_w

from matplotlib import pyplot as plt
from ase.units import Hartree, Bohr
# Analytical Drude result
n = 1 / (a1.get_volume() * Bohr**-3)

wp = np.sqrt(4 * np.pi * n)

plt.plot(w_w * Hartree, -(1. / df1LFCx).imag)
plt.plot(w_w * Hartree, -(1. / df2LFCx).imag)
plt.plot([wp1 * Hartree, wp1 * Hartree], [-1, 3])
plt.plot([wp2 * Hartree, wp2 * Hartree], [-1, 3])
plt.plot([wp * Hartree, wp * Hartree], [-1, 3], '--')
print(wp, wp1, wp2)

plt.show()

w1, I1 = findpeak(w_w, -(1. / df1LFCx).imag)
w2, I2 = findpeak(w_w, -(1. / df2LFCx).imag)
equal(w1, w2, 1e-2)
equal(I1, I2, 1e-3)

w1, I1 = findpeak(w_w, -(1. / df1LFCy).imag)
w2, I2 = findpeak(w_w, -(1. / df2LFCy).imag)
equal(w1, w2, 1e-2)
equal(I1, I2, 1e-3)

w1, I1 = findpeak(w_w, -(1. / df1LFCz).imag)
w2, I2 = findpeak(w_w, -(1. / df2LFCz).imag)
equal(w1, w2, 1e-2)
equal(I1, I2, 1e-3)
