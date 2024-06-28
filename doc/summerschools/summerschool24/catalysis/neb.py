# %%
# teacher
import ase.visualize as viz
viz.view = lambda atoms, repeat=None: None

# %%
"""
## Introduction to Nudged Elastic Band (NEB) calculations

This tutorial describes how to use the NEB method to calculate the diffusion
barrier for an Au atom on Al(001). If you are not familiar with the NEB
method some relevant references are listed
[here.](https://wiki.fysik.dtu.dk/ase/ase/neb.html)

The tutorial uses the EMT potential in stead of DFT, as this is a lot faster.
It is based on a [tutorial found on the ASE
webpage](https://wiki.fysik.dtu.dk/ase/tutorials/neb/diffusion.html#diffusion-tutorial).

"""

# %%
# magic: %matplotlib notebook
import matplotlib.pyplot as plt
from ase.build import fcc100, add_adsorbate
from ase.constraints import FixAtoms
from ase.calculators.emt import EMT
from ase.optimize import BFGS
from ase.visualize import view
from ase.mep import NEB
from ase.io import read

# %%
"""
First we set up the initial state and check that it looks ok:
"""

# %%
# 2x2-Al(001) surface with 3 layers and an
# Au atom adsorbed in a hollow site:
slab = fcc100('Al', size=(2, 2, 3))
add_adsorbate(slab, 'Au', 1.7, 'hollow')
slab.center(axis=2, vacuum=4.0)
view(slab)

# %%
"""
Then we optimise the structure and save it
"""

# %%
# Fix second and third layers:
mask = [atom.tag > 1 for atom in slab]
slab.set_constraint(FixAtoms(mask=mask))

# Use EMT potential:
slab.calc = EMT()

# Optimise initial state:
qn = BFGS(slab, trajectory='initial.traj')
qn.run(fmax=0.05)

# %%
"""
We make the final state by moving the Au atom one lattice constant and
optimise again.
"""

# %%
# Optimise final state:
slab[-1].x += slab.get_cell()[0, 0] / 2
qn = BFGS(slab, trajectory='final.traj')
qn.run(fmax=0.05)

# %%
"""
Now we make a NEB calculation with 3 images between the inital and final
states. The images are initially made as copies of the initial state and the
command `interpolate()` makes a linear interpolation between the initial and
final state. As always, we check that everything looks ok before we run the
calculation.

NOTE: The linear interpolation works well in this case but not for e.g.
rotations. In this case an improved starting guess can be made with the [IDPP
method.](https://wiki.fysik.dtu.dk/ase/tutorials/neb/idpp.html#idpp-tutorial)
"""

# %%
initial = read('initial.traj')
final = read('final.traj')

constraint = FixAtoms(mask=[atom.tag > 1 for atom in initial])

n_im = 3  # number of images
images = [initial]
for i in range(n_im):
    image = initial.copy()
    image.calc = EMT()
    image.set_constraint(constraint)
    images.append(image)

images.append(final)

neb = NEB(images)
neb.interpolate()
view(images)

# %%
qn = BFGS(neb, trajectory='neb.traj')
qn.run(fmax=0.05)

# %%
"""
We visualize the final path with:
"""

# %%
view(images)

# %%
"""
You can find the barrier by selecting Tools->NEB in the gui (unfortunately,
the gui cannot show graphs when started from a notebook), or you can make a
script using
[NEBTools](https://wiki.fysik.dtu.dk/ase/ase/neb.html#ase.neb.NEBTools),
e.g.:
"""

# %%
from ase.mep import NEBTools

images = read('neb.traj@-5:')

nebtools = NEBTools(images)

# Get the calculated barrier and the energy change of the reaction.
Ef, dE = nebtools.get_barrier()
print('Barrier:', Ef)
print('Reaction energy:', dE)

# Generate new matplotlib axis - otherwise nebtools plot double.
fig, ax = plt.subplots()
nebtools.plot_band(ax)

# %%
"""
## Exercise

Now you should make your own NEB using the configuration with N<sub>2</sub>
lying down as the initial state and the configuration with two N atoms
adsorbed on the surface as the final state. The NEB needs to run in parallel
so you should make it as a python script, however you can use the Notebook to
test your configurations (but not the parallelisation) if you like and export
it as a script in the end.

### Parallelisation

The NEB should be parallelised over images. An example can be found in [this
GPAW tutorial](https://wiki.fysik.dtu.dk/gpaw/tutorialsexercises/moleculardynamics/neb/neb.html). The
script enumerates the cpu's and uses this number (``rank``) along with the
total number of cpu's (``size``) to distribute the tasks.
"""

# %%
# This code is just for illustration
from ase.parallel import world
n_im = 4              # Number of images
n = world.size // n_im      # number of cpu's per image
j = 1 + world.rank // n     # image number on this cpu
assert n_im * n == world.size

# %%
"""
For each image we assign a set of cpu's identified by their rank. The rank
numbers are given to the calculator associated with this image.
"""

# %%
# This code is just for illustration
from gpaw import GPAW, PW
images = [initial]
for i in range(n_im):
    ranks = range(i * n, (i + 1) * n)
    image = initial.copy()
    image.set_constraint(constraint)
    if world.rank in ranks:
        calc = GPAW(mode=PW(350),
                    nbands='130%',
                    xc='PBE',  # student: ...,
                    txt=f'{i}.txt',
                    communicator=ranks)
        image.calc = calc
    images.append(image)
images.append(final)

# %%
"""
When running the parallel NEB, you should choose the number of CPU cores
properly.  Let Ncore = N_im * Nk where N_im is the number of images, and Nk
is a divisor of the number of k-points; i.e. if there are 6 irreducible
k-point, Nk should be 1, 2, 3 or 6.  Keep the total number of cores to 24 or
less, or your job will wait too long in the queue.
"""

# %%
"""
### Input parameters

Some suitable parameters for the NEB are given below:

* Use the same calculator and constraints as for the initial and final images, but remember to set the `communicator` as described above
* Use 6 images. This gives a reasonable description of the energy landscape and can be run e.g. on 12 cores.
* Use a spring constant of 1.0 between the images. A lower value will slow the convergence
* Relax the initial NEB until `fmax=0.1eV/Å`, then switch on the climbing image and relax until `fmax=0.05eV/Å`.
"""

# %%
# teacher:
from gpaw import GPAW, PW
from ase.visualize import view
from ase.optimize import BFGS
from ase.parallel import world

initial = read('N2Ru.traj')
final = read('2Nads.traj')

images = [initial]
N = 4  # No. of images

z = initial.positions[:, 2]
constraint = FixAtoms(mask=(z < z.min() + 1.0))

j = world.rank * N // world.size
n = world.size // N  # number of cpu's per image

for i in range(N):
    ranks = range(i * n, (i + 1) * n)
    image = initial.copy()
    image.set_constraint(constraint)
    if world.rank in ranks:
        calc = GPAW(xc='PBE',
                    mode=PW(350),
                    nbands='130%',
                    communicator=ranks,
                    txt=f'{i}.txt',
                    kpts={'size': (4, 4, 1), 'gamma': True},
                    convergence={'eigenstates': 1e-7})
        image.calc = calc
    images.append(image)
images.append(final)

neb = NEB(images, k=1.0, parallel=True)
neb.interpolate()
qn = BFGS(neb, logfile='neb.log', trajectory='neb.traj')
qn.run(fmax=0.1)
neb.climb = True
qn.run(fmax=0.05)

# %%
"""
Once the calculation is done you should check that the final path looks
reasonable. What is the N-N distance at the saddle point? Use NEBTools to
calculate the barrier. Is N<sub>2</sub> likely to dissociate on the surface
at room temperature?
"""

# %%
# teacher:
# Check results:
images[1:5] = read('neb.traj@-5:-1')
energies = [image.get_potential_energy() for image in images]
emax = max(energies)
assert energies[2] == emax
assert abs(emax - energies[0] - 1.32) < 0.02
d = images[2].get_distance(-1, -2)
assert abs(d - 1.777) < 0.004
