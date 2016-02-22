# creates: HLi.png
import matplotlib.pyplot as plt
from ase import Atoms
from ase.units import Bohr
from gpaw.utilities.wavefunctions import WaveFunctionInterpolator
from gpaw import GPAW

hli = Atoms('HLi', positions=[[0, 0, 0], [0, 0, 1.6]])
hli.center(vacuum=2.5)
hli.set_calculator(GPAW(txt='HLi.txt'))
hli.get_potential_energy()

wfs = WaveFunctionInterpolator(hli.calc, h=0.05)

for n, color in enumerate(['green', 'red']):
    ps = wfs.get_wave_function(n, ae=False)
    ae = wfs.get_wave_function(n)
    norm = wfs.gd.integrate(ae**2)
    print('Norm:', norm)
    assert abs(norm - 1) < 1e-2
    i = ps.shape[0] // 2
    x = wfs.gd.coords(2) * Bohr
    plt.plot(x, ps[i, i], '--', color=color, label=r'$\tilde\psi_{}$'.format(n))
    plt.plot(x, ae[i, i], '-', color=color, label=r'$\psi_{}$'.format(n))

    ps0 = hli.calc.get_pseudo_wave_function(n, pad=True) * Bohr**1.5
    gd0 = hli.calc.wfs.gd
    i = ps0.shape[0] // 2
    x = gd0.coords(2) * Bohr
    plt.plot(x, ps0[i, i], 'o', color=color)
    
plt.plot(x, 0 * x, 'k')
plt.xlabel('z [Ang]')
plt.legend()
plt.savefig('HLi.png')
