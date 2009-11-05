import numpy as np
import matplotlib
matplotlib.use('Agg')
import pylab as plt
from ase.utils.eos import EquationOfState
from ase.io import read
def f(width, k, g):
    filename = 'Fe-B-%.2f-%02d-%2d.traj' % (width, k, g)
    configs = read(filename + '@::2')
    # Extract volumes and energies:
    volumes = [a.get_volume() for a in configs]
    energies = [a.get_potential_energy() for a in configs]
    eos = EquationOfState(volumes, energies)
    v0, e0, B = eos.fit()
    return v0, e0, B

kk = [2, 4, 6, 8, 10, 12]

plt.figure(figsize=(6, 4))
for width in [0.05, 0.1, 0.15, 0.2]:
    a = [f(width, k, 16)[0]**(1.0 / 3.0) for k in kk]
    print ('%7.3f ' * 7) % ((width,) + tuple(a))
    plt.plot(kk, a, label='width = %.2f eV' % width)
plt.legend(loc='lower right')
plt.axis(ymin=2.83, ymax=2.85)
plt.xlabel('number of k-points')
plt.ylabel('lattice constant [Ang]')
plt.savefig('Fe_conv_k.png')

plt.figure(figsize=(6, 4))
gg = np.arange(12, 40, 4)
a = [f(0.1, 6, g)[0]**(1.0 / 3.0) for g in gg]
plt.plot(2.84 / gg, a, 'o-')
plt.axis(ymin=2.83, ymax=2.85)
plt.xlabel('grid-spacing [Ang]')
plt.ylabel('lattice constant [Ang]')
plt.savefig('Fe_conv_h.png')
