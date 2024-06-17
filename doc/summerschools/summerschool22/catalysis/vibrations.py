# %%
"""
Extra exercise - vibrational energy
===================================

The energy calculated with DFT is the electronic energy at 0K. However, to
model catalytic reactions we are usually intereseted in the energy landscape
at finite temperature. In this exercise we will calculate the energy
contributions from vibrations and see how they affect the splitting of
N<sub>2</sub> on Ru.

We calculate the vibrational energy in the harmonic approximation using the
finite displacement method. For further reading see for example:

* [Stoffel et al, Angewandte Chemie Int. Ed., 49, 5242-5266 (2010)][1]
* [Togo and Tanaka, Scripta Materialia 108, 1-5 (2015)][2]

[1]: https://onlinelibrary.wiley.com/doi/abs/10.1002/anie.200906780
[2]: https://www.sciencedirect.com/science/article/pii/S1359646215003127


Vibrational energy of the initial and final states
--------------------------------------------------

a) Calculating vibrations requires tighter convergence than normal energy
   calculations. Therefore you should first take your already optimised initial
   and final state geometries from the NEB calculations and relax them further
   to `fmax=0.01eV/Å` with the QuasiNewton optimiser and an energy cutoff of
   500eV. Converge the eigenstates to 1e-8. (Note that for other systems you
   might need even tighter convergence!)

Submit the structures to the queue. The optimisation should take 10-15 mins
for each structure on 8 cores.
"""

# %%
# teacher:
# Script for reading previous structure and converging
from gpaw import GPAW, PW
from ase.constraints import FixAtoms
from ase.optimize import QuasiNewton
import numpy as np
from ase.io import read, write

for name in ['N2Ru-top.traj', '2Nads.traj']:
    slab = read(name)

    z = slab.positions[:, 2]
    constraint = FixAtoms(mask=(z < z.min() + 1.0))
    slab.set_constraint(constraint)

    calc = GPAW(xc='PBE',
                mode=PW(500),
                txt='tight-' + name[:-4] + 'txt',
                kpts={'size': (4, 4, 1), 'gamma': True},
                convergence={'eigenstates': 1e-8})
    slab.calc = calc

    dyn = QuasiNewton(slab,
                      trajectory='tight-' + name,
                      maxstep=0.02,
                      logfile='tight-' + name[:-4] + 'log')
    dyn.run(fmax=0.01)

# %%
"""
b) Once you have done this you can calculate the vibrations using the
   [vibrations module in ASE][3] following
   the template below. We only calculate the vibrations of the adsorbate and
   assume that the frequencies of the substrate are unchanged - this is a
   common assumption. Use 4 displacements to fit the frequencies and the same
   calculator parameters as in a).

[3]: https://wiki.fysik.dtu.dk/ase/ase/vibrations/vibrations.html

Submit the calculations for the initial and the final state to the queue. It
will take a while to run, but you can start preparing your analysis (part c
and d) or read some of the references in the meantime.
"""

# %%
# You can't run this cell - use it as a starting point for your python script!
from ase.io import read
from ase.vibrations import Vibrations
from gpaw import GPAW, PW

slab = read('tight-N2Ru-top.traj')  # student: slab = read('your_structure')
calc = GPAW(xc='PBE',
            mode=PW(500),
            kpts={'size': (4, 4, 1), 'gamma': True},  # student: kpts=?,
            symmetry={'point_group': False},
            txt='vib.txt')
slab.calc = calc
Uini = slab.get_potential_energy()

vib = Vibrations(slab,
                 name='vib',
                 indices=[8, 9],  # student: indices=[?, ?],
                 nfree=4)
vib.run()
vib.summary(log='vib_summary.log')
for i in range(6):
    vib.write_mode(i)

# %%
from ase.io import read
from ase.vibrations import Vibrations
from gpaw import GPAW, PW

slab = read('tight-2Nads.traj')
calc = GPAW(xc='PBE',
            mode=PW(500),
            kpts={'size': (4, 4, 1), 'gamma': True},  # student: kpts=?,
            symmetry={'point_group': False},
            txt='vib2.txt')
slab.calc = calc
Ufin = slab.get_potential_energy()

vib2 = Vibrations(slab,
                  name='vib2',
                  indices=[8, 9],  # student: indices=[?, ?],
                  nfree=4)
vib2.run()
vib2.summary(log='vib2_summary.log')
for i in range(6):
    vib2.write_mode(i)

# %%
"""

The final three lines write out the frequencies and `.traj` files with
animations of the different modes in order of their energy. Take a look at
the vibrational modes in the files. Do you understand why the first mode has
a low energy, while the last one has a high energy?
"""

# %%
"""
c) Analyse the frequencies in the harmonic approximation:
"""

# %%
from ase.thermochemistry import HarmonicThermo

T = 300  # Kelvin
# An example only - insert your calculated energy levels here - in eV:
energies = [0.01, 0.05, 0.10]
vibs = HarmonicThermo(energies)
Efree = vibs.get_helmholtz_energy(T, verbose=True)
print('Free energy at 300 K: ', Efree)

# %%
"""
The `verbose` keyword gives a detailed description of the different
contributions to the free energy. For more information on what the different
contributions are see the [ASE background webpage][4]
(go to the **Harmonic limit** sub-heading).

[4]: https://wiki.fysik.dtu.dk/ase/ase/thermochemistry/thermochemistry.html

Now try to calculate how the different contributions change with temperature.
You can for example make a `for` loop and use the `get_entropy()` and
`get_internal_energy()` methods [(see description here)][5].

[5]: https://wiki.fysik.dtu.dk/ase/ase/thermochemistry/thermochemistry.html#ase.thermochemistry.IdealGasThermo.get_enthalpy
"""

# %%
"""
d) Calculate how the vibrational energy affects the overall reaction energy.
   How does it change with temperature? Which contribution is important for the
   change in reaction energy?
"""

# %%
# teacher:
# Script for analysing energy contributions from vibrations
from ase.io import read
from ase.thermochemistry import HarmonicThermo
from matplotlib import pyplot as plt

energies = vib.get_energies()  # initial vib. energy levels
energiesfin = vib2.get_energies()  # final vib. energy levels

# The actual analysis
ads = HarmonicThermo(energies)
adsfin = HarmonicThermo(energiesfin)
print('DFT reaction energy: ', Ufin - Uini)
Sini = []
Sfin = []
U = []
U2 = []
temp = [10, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500]
for T in temp:
    Sini.append(-T * ads.get_entropy(T, verbose=False))
    U.append(ads.get_internal_energy(T, verbose=False))
    Sfin.append(-T * adsfin.get_entropy(T, verbose=False))
    U2.append(adsfin.get_internal_energy(T, verbose=False))
    Efree = ads.get_helmholtz_energy(T, verbose=False)
    Efreefin = adsfin.get_helmholtz_energy(T, verbose=False)
    print(f'Reaction free energy at {T} K: ',
          Ufin + Efreefin - (Uini + Efree))
plt.plot(temp, Sini, 'r', label='-T*S ini')
plt.plot(temp, U, 'b', label='U ini')
plt.plot(temp, Sfin, 'r--', label='-T*S fin')
plt.plot(temp, U2, 'b--', label='U fin')
plt.plot(temp, np.array(U) + np.array(Sini), 'm', label='Efree ini')
plt.plot(temp, np.array(U2) + np.array(Sfin), 'm--', label='Efree fin')
plt.legend()
plt.savefig('vib.png')

# %%
"""
e) To make sure that your NEB is converged you should also calculate the
   vibrational energy of the transition state. Again, this requires tighter
   convergence than we have used in the NEB exercise. This takes a while to run
   so to save time, we provide the transition state geometry from a reasonably
   converged NEB (i.e. `fmax=0.01`, a cutoff energy of 800eV and 6x6
   k-points) in the file `TS.xyz`. Calculate the vibrations with these
   parameters. How many imaginary modes do you get and how do they look? What
   does this mean?
"""

# %%
# teacher:
from ase.optimize import BFGS
from ase.mep import NEB
from ase.io import Trajectory

initial = read('N2Ru.traj')
final = read('2Nads.traj')

images = read('neb.traj@-5:-2')
constraint = FixAtoms(indices=list(range(4)))
ts = images[1]
ts.set_constraint(constraint)
calc = GPAW(xc='PBE',
            mode=PW(800),
            txt='ts.txt',
            kpts={'size': (6, 6, 1), 'gamma': True})
ts.calc = calc

neb = NEB(images, k=1.0, climb=True)
qn = BFGS(neb, logfile='neb.log')
traj = Trajectory('ts.traj', 'w', ts)
qn.attach(traj)
qn.run(fmax=0.01)
del ts.info['adsorbate_info']  # xyz-format doesn't support that
write('TS.xyz', ts)


# %%
# teacher:
ts.calc = GPAW(xc='PBE',
               mode=PW(800),
               txt='vibts.txt',
               kpts={'size': (6, 6, 1), 'gamma': True},
               symmetry={'point_group': False})
vib = Vibrations(ts, name='vibts', indices=(8, 9), nfree=4, delta=0.02)
vib.run()
vib.summary(log='vibts_summary.log')
for i in range(6):
    vib.write_mode(i)

# ---------------------
#   #    meV     cm^-1
# ---------------------
#   0   72.3i    583.1i
#   1    5.4      43.9
#   2   41.7     336.0
#   3   47.4     382.3
#   4   70.1     565.7
#   5   71.3     575.0
# ---------------------
# Zero-point energy: 0.118 eV
# The imaginary mode is beautifully along the reaction coordinate!
