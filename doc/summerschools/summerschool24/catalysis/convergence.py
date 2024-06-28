# %%
# teacher
import ase.visualize as viz
viz.view = lambda atoms, repeat=None: None

# %%
"""
Convergence checks
==================

In this notebook we look at the adsorption energy and height of a nitrogen
atom on a Ru(0001) surface in the hcp site.  We check for convergence with
respect to:

* number of layers
* number of k-points in the BZ
* plane-wave cutoff energy
"""

# %%
"""
Nitrogen atom
-------------

First step is an isolated nitrogen atom which has a magnetic moment of 3.
More information:
[Atoms][1] and [GPAW parameters][2].

[1]: https://wiki.fysik.dtu.dk/ase/ase/atoms.html#ase.Atoms
[2]: https://wiki.fysik.dtu.dk/gpaw/documentation/basic.html#parameters
"""

# %%
ecut = 400.0
vacuum = 4.0
from ase import Atoms
from gpaw import GPAW, PW, Davidson
nitrogen = Atoms('N', magmoms=[3])
nitrogen.center(vacuum=4.0)
nitrogen.calc = GPAW(txt='N.txt',
                     mode=PW(ecut),
                     eigensolver=Davidson(niter=2))
en = nitrogen.get_potential_energy()

# %%
print(en, 'eV')

# %%
"""
Clean slab
----------

We use the [ase.build.hcp0001()][3] function to build the Ru(0001) surface.

[3]: https://wiki.fysik.dtu.dk/ase/ase/build/surface.html#ase.build.hcp0001
"""

# %%
nlayers = 2
a = 2.72
c = 1.58 * a
from ase.build import hcp0001
slab = hcp0001('Ru', a=a, c=c, size=(1, 1, nlayers), vacuum=vacuum)

# %%
from ase.visualize import view
view(slab, repeat=(3, 3, 2))

# %%
nkpts = 7
slab.calc = GPAW(txt='Ru.txt',
                 mode=PW(ecut),
                 eigensolver=Davidson(niter=2),
                 kpts={'size': (nkpts, nkpts, 1), 'gamma': True},
                 xc='PBE')
eru = slab.get_potential_energy()

# %%
"""
# N/Ru(0001):
"""

# %%
import numpy as np
height = 1.1
nslab = hcp0001('Ru', a=a, c=c, size=(1, 1, nlayers))
# Calculate the coordianates of the N-atoms:
z = slab.positions[:, 2].max() + height
x, y = np.dot([2 / 3, 2 / 3], slab.cell[:2, :2])
nslab.append('N')
nslab.positions[-1] = [x, y, z]
nslab.center(vacuum=vacuum, axis=2)  # 2: z-axis

# %%
"""
Alternatively, you can just use the [add_adsorbate()][4] function:

[4]: https://wiki.fysik.dtu.dk/ase/ase/build/surface.html#ase.build.add_adsorbate
"""

# %%
height = 1.1
nslab = hcp0001('Ru', a=a, c=c, size=(1, 1, nlayers))
from ase.build import add_adsorbate
add_adsorbate(nslab, 'N', position='hcp', height=height)
nslab.center(vacuum=vacuum, axis=2)

# %%
view(nslab, repeat=(3, 3, 2))

# %%
from ase.io import write
write('nru2.png', nslab.repeat((3, 3, 1)))

# %%
"""
![rnu](nru2.png)
"""

# %%
nslab.calc = GPAW(txt='NRu.txt',
                  mode=PW(ecut),
                  eigensolver=Davidson(niter=2),
                  poissonsolver={'dipolelayer': 'xy'},
                  kpts={'size': (nkpts, nkpts, 1), 'gamma': True},
                  xc='PBE')
enru0 = nslab.get_potential_energy()

# %%
print('Unrelaxed adsoption energy:', enru0 - eru - en, 'eV')

# %%
nslab.get_forces()

# %%
"""
The force on the N-atom is quite big.  Let's freeze the surface and relax the
adsorbate.  We use
[ase.optimize.BFGSLineSearch][5] and
[ase.constraints.FixAtoms][6]
for this task.

[5]: https://wiki.fysik.dtu.dk/ase/ase/optimize.html#module-ase.optimize
[6]: https://wiki.fysik.dtu.dk/ase/ase/constraints.html#ase.constraints.FixAtoms
"""

# %%
# This cell will take a few minutes to finish ...
from ase.constraints import FixAtoms
from ase.optimize import BFGSLineSearch
nslab.constraints = FixAtoms(indices=list(range(nlayers)))
optimizer = BFGSLineSearch(nslab, trajectory='NRu.traj')
optimizer.run(fmax=0.01)
height = nslab.positions[-1, 2] - nslab.positions[:-1, 2].max()
print('Height:', height, 'Ang')

# %%
enru = nslab.get_potential_energy()
print('Relaxed adsorption energy:', enru - eru - en, 'eV')

# %%
"""
In order to make it easy to check for convergence of the adsorption energy
and height we write a little function that does all of the stuff above taking
`nlayers`, `nkpts` and `ecut` as input parameters.

The `adsorb()` function is shown below for completenes, but you should not
use it inside this notebook.  Instead, please take a look at the
[check_convergence.py][7]
script that also contains the definition of the `adsorb()` function.  The
script will do a bunch of calculations with different parameters and store
the results in a database file (`convergence.db`) that we analyse below ...

[7]: https://gitlab.com/gpaw/gpaw/blob/master/doc/summerschools/summerschool18/catalysis/check_convergence.py
"""


# %%
def adsorb(db, height=1.2, nlayers=3, nkpts=7, ecut=400):
    """Adsorb nitrogen in hcp-site on Ru(0001) surface.

    Do calculations for N/Ru(0001), Ru(0001) and a nitrogen atom
    if they have not already been done.

    db: Database
        Database for collecting results.
    height: float
        Height of N-atom above top Ru-layer.
    nlayers: int
        Number of Ru-layers.
    nkpts: int
        Use a (nkpts * nkpts) Monkhorst-Pack grid that includes the
        Gamma point.
    ecut: float
        Cutoff energy for plane waves.

    Returns height.
    """

    name = f'Ru{nlayers}-{nkpts}x{nkpts}-{ecut:.0f}'

    parameters = dict(mode=PW(ecut),
                      eigensolver=Davidson(niter=2),
                      poissonsolver={'dipolelayer': 'xy'},
                      kpts={'size': (nkpts, nkpts, 1), 'gamma': True},
                      xc='PBE')

    # N/Ru(0001):
    slab = hcp0001('Ru', a=a, c=c, size=(1, 1, nlayers))
    z = slab.positions[:, 2].max() + height
    x, y = np.dot([2 / 3, 2 / 3], slab.cell[:2, :2])
    slab.append('N')
    slab.positions[-1] = [x, y, z]
    slab.center(vacuum=vacuum, axis=2)  # 2: z-axis

    # Fix first nlayer atoms:
    slab.constraints = FixAtoms(indices=list(range(nlayers)))

    id = db.reserve(name=f'N/{nlayers}Ru(0001)', nkpts=nkpts, ecut=ecut)
    if id is not None:  # skip calculation if already done
        slab.calc = GPAW(txt='N' + name + '.txt',
                         **parameters)
        optimizer = BFGSLineSearch(slab, logfile='N' + name + '.opt')
        optimizer.run(fmax=0.01)
        height = slab.positions[-1, 2] - slab.positions[:-1, 2].max()
        db.write(slab, id=id,
                 name=f'N/{nlayers}Ru(0001)', nkpts=nkpts, ecut=ecut,
                 height=height)

    # Clean surface (single point calculation):
    id = db.reserve(name=f'{nlayers}Ru(0001)', nkpts=nkpts, ecut=ecut)
    if id is not None:
        del slab[-1]  # remove nitrogen atom
        slab.calc = GPAW(txt=name + '.txt',
                         **parameters)
        slab.get_forces()
        db.write(slab, id=id,
                 name=f'{nlayers}Ru(0001)', nkpts=nkpts, ecut=ecut)

    # Nitrogen atom:
    id = db.reserve(name='N-atom', ecut=ecut)
    if id is not None:
        # Create spin-polarized nitrogen atom:
        molecule = Atoms('N', magmoms=[3])
        molecule.center(vacuum=4.0)
        # Remove parameters that make no sense for an isolated atom:
        del parameters['kpts']
        del parameters['poissonsolver']
        # Calculate energy:
        molecule.calc = GPAW(txt=name + '.txt', **parameters)
        molecule.get_potential_energy()
        db.write(molecule, id=id, name='N-atom', ecut=ecut)

    return height


# %%
"""
Read more about ASE databases
[here](https://wiki.fysik.dtu.dk/ase/ase/db/db.html#module-ase.db).
"""

# %%
from ase.db import connect
db = connect('convergence.db')

# %%
assert len(db) == 72
# We assume that the calculations have already been done for you

# %%
# This cell should not take any time because the adsorb() function
# is clever enough to skip calculations already in the database.
h = 1.2
for n in range(1, 10):  # layer
    h = adsorb(db, h, n, 7, 400)
for k in range(4, 18):  # k-points
    h = adsorb(db, h, 2, k, 400)
for ecut in range(350, 801, 50):  # plane-wave cutoff
    h = adsorb(db, h, 2, 7, ecut)

# %%
"""
You can inspect database file with the
[command line tool](https://wiki.fysik.dtu.dk/ase/ase/db/db.html#ase-db)
`ase db` like this:
"""

# %%
# magic: !ase db convergence.db -c ++ -L 0  # show all columns (-c ++); show all rows (-L 0)

# %%
# magic: !ase db convergence.db formula=Ru2N,nkpts=7 -c ecut,height -s ecut

# %%
"""
Now we can analyse the results of the convergence tests.  We extract the
result from the database with a little helper function `select()`:
"""


# %%
def select(nlayers, nkpts, ecut):
    """Extract adsorption energy and height from database."""
    en = db.get(N=1, Ru=0, ecut=ecut).energy
    eru = db.get(N=0, Ru=nlayers, nkpts=nkpts, ecut=ecut).energy
    row = db.get(N=1, Ru=nlayers, nkpts=nkpts, ecut=ecut)
    enru = row.energy
    return row.height, enru - eru - en


# %%
# magic: %matplotlib notebook
from matplotlib import pyplot as plt

# %%
n = np.arange(1, 10)
h, e = np.array([select(nlayers, 7, 400) for nlayers in n]).T
fig, axs = plt.subplots(2, 1, sharex=True)
fig.subplots_adjust(hspace=0)
axs[0].plot(n, h)
axs[1].plot(n, e)
axs[0].set_ylabel('height [Å]')
axs[1].set_ylabel('ads. energy [eV]')
axs[1].set_xlabel('number of layers')

# %%
k = np.arange(4, 18)
h, e = np.array([select(2, nkpts, 400) for nkpts in k]).T
fig, axs = plt.subplots(2, 1, sharex=True)
fig.subplots_adjust(hspace=0)
axs[0].plot(k, h)
axs[1].plot(k, e)
axs[0].set_ylabel('height [Å]')
axs[1].set_ylabel('ads. energy [eV]')
axs[1].set_xlabel('number of k-points')

# %%
x = np.arange(350, 801, 50)
h, e = np.array([select(2, 7, ecut) for ecut in x]).T
fig, axs = plt.subplots(2, 1, sharex=True)
fig.subplots_adjust(hspace=0)
axs[0].plot(x, h)
axs[1].plot(x, e)
axs[0].set_ylabel('height [Å]')
axs[1].set_ylabel('ads. energy [eV]')
axs[1].set_xlabel('plane-wave cutoff energy [eV]')

# %%
"""
Conclusion
----------

For accurate calculations you would need:

* a plane-wave cutoff of 600 eV
* 5 layers of Ru
* 9x9 Monkhorst-Pack grid for BZ sampling (for a 1x1 unit cell)

For our quick'n'dirty calculations we will use 350 eV, 2 layers and a
4x4 $\\Gamma$-centered Monkhorst-Pack grid (for a 2x2 unit cell).
"""

# %%
