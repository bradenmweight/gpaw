# web-page: ag-ddos.png
import numpy as np
from scipy.integrate import trapezoid
import matplotlib.pyplot as plt
from gpaw import GPAW

calc = GPAW('au.gpw', txt=None)
energy, pdos = calc.get_orbital_ldos(a=0, angular='d')
energy -= calc.get_fermi_level()
I = trapezoid(pdos, energy)
center = trapezoid(pdos * energy, energy) / I
width = np.sqrt(trapezoid(pdos * (energy - center)**2, energy) / I)
plt.plot(energy, pdos)
plt.xlabel('Energy (eV)')
plt.ylabel('d-projected DOS on atom 0')
plt.title(f'd-band center = {center:.1f} eV, d-band width = {width:.1f} eV')
# plt.show()
plt.savefig('ag-ddos.png')
