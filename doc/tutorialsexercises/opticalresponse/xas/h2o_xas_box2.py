import numpy as np
import matplotlib.pyplot as plt

from gpaw import GPAW, setup_paths
from gpaw.xas import XAS

setup_paths.insert(0, '.')

h = 0.2

offset = 0.0
for L in np.arange(4, 14, 2) * 8 * h:
    calc = GPAW(f'h2o_hch_{L:.1f}.gpw')
    xas = XAS(calc)
    x, y = xas.get_spectra(fwhm=0.4)
    plt.plot(x, sum(y) + offset, label=f'{L:.1f}')
    offset += 0.005

plt.legend()
plt.xlim(-6, 4)
plt.ylim(-0.002, 0.03)
plt.savefig('h2o_xas_box.png')
