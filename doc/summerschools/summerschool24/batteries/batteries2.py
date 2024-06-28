# %%
# teacher
import ase.visualize as viz
viz.view = lambda atoms, repeat=None: None

# %%
"""
# Battery project
"""

# %%
"""
## Day 3 - Equilibrium potential
"""

# %%
"""
Today you will study the LiFePO$_4$ cathode. You will calculate the equilibrium potential and use Bayesian error estimation to quantify how sensitive the calculated equilibrium potential is towards choice of functional. After today you should be able to discuss:

-  The volume change during charge/discharge.

-  The maximum gravimetric and volumetric energy density of a FePO$_4$/C battery assuming the majority of weight and volume will be given by the electrodes.

-  Uncertainty in the calculations.

Some of calculations you will perform today will be tedious to be run in this notebook. You will automatically submit some calculations to the HPC cluster directly from this notebook. When you have to wait for calculations to finish you can get started on addressing the bullet points above.
"""

# %%
"""
## Initialize
"""

# %%
# magic: %matplotlib notebook
import numpy as np
from ase.visualize import view
import matplotlib.pyplot as plt
from ase.io import read, write, Trajectory
from ase.parallel import paropen
from gpaw import GPAW, FermiDirac, Mixer, PW
from ase.dft.bee import BEEFEnsemble
from ase import Atoms

# %%
"""
## FePO$_4$
"""

# %%
"""
First we will construct an atoms object for FePO$_4$. ASE can read files from in a large number of different [formats](https://wiki.fysik.dtu.dk/ase/ase/io/io.html?highlight=read%20formats#file-input-and-output). However, in this case you will build it from scratch using the below information:
"""

# %%
# Positions:
# Fe      2.73015081       1.46880951       4.56541172
# Fe      2.23941067       4.40642872       2.14957739
# Fe      7.20997230       4.40642925       0.26615813
# Fe      7.70070740       1.46880983       2.68199421
# O       1.16033403       1.46881052       3.40240205
# O       3.80867172       4.40642951       0.98654342
# O       8.77981469       4.40642875       1.42923946
# O       6.13142032       1.46881092       3.84509827
# O       4.37288562       1.46880982       0.81812712
# O       0.59764596       4.40643021       3.23442747
# O       5.56702590       4.40642886       4.01346264
# O       9.34268360       1.46880929       1.59716233
# O       1.64001691       0.26061277       1.17298291
# O       3.32931769       5.61463705       3.58882629
# O       8.30013707       3.19826250       3.65857000
# O       6.61076951       2.67698811       1.24272700
# O       8.30013642       5.61459688       3.65856912
# O       6.61076982       0.26063178       1.24272567
# O       1.64001666       2.67700652       1.17298270
# O       3.32931675       3.19822249       3.58882660
# P       0.90585688       1.46880966       1.89272372
# P       4.06363530       4.40642949       4.30853266
# P       9.03398503       4.40642957       2.93877879
# P       5.87676435       1.46881009       0.52297232


# Unit cell:
#            periodic     x          y          z
#   1. axis:    yes    9.94012    0.00000    0.00000
#   2. axis:    yes    0.00000    5.87524    0.00000
#   3. axis:    yes    0.00000    0.00000    4.83157


# %%
"""
You *can* use the cell below as a starting point.
"""

# %%
# fepo4 = Atoms('Fe4O...',
#               positions=[[x0, y0, z0],[x1, y1, z1]...],
#               cell=[x, y, z],
#               pbc=[True, True, True])

# Teacher:
fepo4 = Atoms('Fe4O16P4',
             positions=[[2.73015081, 1.46880951, 4.56541172],
                [2.23941067, 4.40642872, 2.14957739],
                [7.20997230, 4.40642925, 0.26615813],
                [7.70070740, 1.46880983, 2.68199421],
                [1.16033403, 1.46881052, 3.40240205],
                [3.80867172, 4.40642951, 0.98654342],
                [8.77981469, 4.40642875, 1.42923946],
                [6.13142032, 1.46881092, 3.84509827],
                [4.37288562, 1.46880982, 0.81812712],
                [0.59764596, 4.40643021, 3.23442747],
                [5.56702590, 4.40642886, 4.01346264],
                [9.34268360, 1.46880929, 1.59716233],
                [1.64001691, 0.26061277, 1.17298291],
                [3.32931769, 5.61463705, 3.58882629],
                [8.30013707, 3.19826250, 3.65857000],
                [6.61076951, 2.67698811, 1.24272700],
                [8.30013642, 5.61459688, 3.65856912],
                [6.61076982, 0.26063178, 1.24272567],
                [1.64001666, 2.67700652, 1.17298270],
                [3.32931675, 3.19822249, 3.58882660],
                [0.90585688, 1.46880966, 1.89272372],
                [4.06363530, 4.40642949, 4.30853266],
                [9.03398503, 4.40642957, 2.93877879],
                [5.87676435, 1.46881009, 0.52297232]
                        ],
             cell=[9.94012, 5.87524, 4.83157],
             pbc=[1, 1, 1])

# %%
"""
Visualize the structure you have made. Explore the different functions in the visualizer and determine the volume of the cell (`View -> Quick Info`).
"""

# %%
view(fepo4)

# %%
"""
For better convergence of calculations you should specify initial magnetic moments to iron. The iron will in this structure be Fe$^{3+}$ as it donates two *4s* electrons and one *3d* electron to PO$_4$$^{3-}$. What is the magnetic moment of iron? For simplicity you should assume that FePO$_4$ is ferromagnetic.
"""

# %%
# Teacher:
for atom in fepo4:
    if atom.symbol == 'Fe':
        atom.magmom = 5.0  # student: atom.magmom = ?

# %%
"""
Now examine the initial magnetic moments of the system using an [appropriate method](https://wiki.fysik.dtu.dk/ase/ase/atoms.html?highlight=get_initial#list-of-all-methods).
"""

# %%
magmoms = fepo4.get_initial_magnetic_moments()  # student: magmoms = fepo4.xxx()
print(magmoms)

# %%
"""
Write your atoms object to file.
"""

# %%
write('fepo4.traj', fepo4)

# %%
"""
For this calculation you will use the BEEF-vdW functional developed by [Wellendorff et al.](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.85.235149) Although there are better alternatives for calculating the energy of bulk systems, the BEEF-vdW has a build-in ensemble for error estimation of calculated energies. In the set-up of this calculator you will append relevant keyword values into a dictionary, which is inputted the calculator object.
"""

# %%
params_GPAW = {}

# %%
"""
To save computational time while keeping the calculations physically sound, the following should be used:
"""

# %%
params_GPAW['mode']        = PW(500)                      # The used plane wave energy cutoff
params_GPAW['nbands']      = -40                          # The number on empty bands had the system been spin-paired
params_GPAW['kpts']        = {'size': (2, 4, 5),          # The k-point mesh
                              'gamma': True}
params_GPAW['spinpol']     = True                         # Performing spin polarized calculations
params_GPAW['xc']          = 'BEEF-vdW'                   # The used exchange-correlation functional
params_GPAW['occupations'] = FermiDirac(width=0.1,        # The smearing
                                        fixmagmom=True)   # Total magnetic moment fixed to the initial value
params_GPAW['convergence'] = {'eigenstates': 1.0e-4,      # eV^2 / electron
                              'energy':      2.0e-4,      # eV / electron
                              'density':     1.0e-3}
params_GPAW['mixer']       = Mixer(0.1, 5, weight=100.0)  # The mixer used during SCF optimization


# %%
"""
DFT suffers from a so-called self-interaction error. An electron interacts with the system electron density, to which it contributes itself. The error is most pronounced for highly localized orbitals. [Hubbard U correction](https://wiki.fysik.dtu.dk/gpaw/tutorialsexercises/energetics/hubbardu/hubbardu.html) is used to mitigate the self-interaction error of the highly localized *3d*-electrons of Fe. This is done in GPAW using the `setups` keyword.
"""

# %%
params_GPAW['setups']      = {'Fe': ':d,4.3'}             # U=4.3 applied to d orbitals


# %%
"""
Make a GPAW calculator and attach it to the atoms object. Here you will use [get_potential_energy](https://wiki.fysik.dtu.dk/ase/ase/atoms.html#ase.Atoms.get_potential_energy) to start the calculation.
"""

# %%
calc = GPAW(**params_GPAW)
fepo4.calc = calc
epot_fepo4_cell = fepo4.get_potential_energy()
print(epot_fepo4_cell)
write('fepo4_out.traj', fepo4)


# %%
"""
You will use the ensemble capability of the BEEF-vdW functional. You will need this later so you should write it to file so you do not have to start all over again later. Start by obtaining the required data from the calculator, i.e., the individual energy of each term in the BEEF-vdW functional expansion. Get the energy difference compared to BEEF-vdW for 2000 ensemble functionals.
"""

# %%
from ase.dft.bee import BEEFEnsemble

ens = BEEFEnsemble(calc)
dE = ens.get_ensemble_energies(2000)


# %%
"""
Print the energy differences to file. This is not the most efficient way of printing to file but can allow easier subsequent data treatment.
"""

# %%
with paropen('ensemble_fepo4.dat', 'a') as result:
    for e in dE:
        print(e, file=result)


# %%
"""
You now have what you need to make a full script. Make it in the cell below and execute it to make sure the script runs. Once you have made sure the calculation is able to run, stop it by `interupt the kernel`.
"""

# %%
# %%writefile 'fepo4.py'
from ase.parallel import paropen
from ase.io import read, write
from ase.dft.bee import BEEFEnsemble
from gpaw import GPAW, FermiDirac, Mixer, PW

# Read in the structure you made and wrote to file above
fepo4 = read('fepo4.traj')

params_GPAW = {...}

# do calculation ...
# BEEF ...
# write ensemble_fepo4.dat file ...

write('fepo4_out.traj', fepo4)


# Teacher:
from ase.parallel import paropen
from ase.io import read, write
from ase.io.trajectory import Trajectory
from ase.dft.bee import BEEFEnsemble
from gpaw import GPAW, FermiDirac, Mixer, PW
fepo4 = read('fepo4.traj')

params_GPAW = {}
params_GPAW['mode']        = PW(500)                     #The used plane wave energy cutoff
params_GPAW['nbands']      = -40                           #The number on empty bands had the system been spin-paired
params_GPAW['kpts']        = {'size':  (2,4,5),            #The k-point mesh
                              'gamma': True}
params_GPAW['spinpol']     = True                          #Performing spin polarized calculations
params_GPAW['xc']          = 'BEEF-vdW'                    #The used exchange-correlation functional
params_GPAW['occupations'] = FermiDirac(width = 0.1,      #The smearing
                                        fixmagmom = True)  #Total magnetic moment fixed to the initial value
params_GPAW['convergence'] = {'eigenstates': 1.0e-4,       #eV^2 / electron
                              'energy':      2.0e-4,       #eV / electron
                              'density':     1.0e-3,}
params_GPAW['mixer']       = Mixer(0.1, 5, weight=100.0)   #The mixer used during SCF optimization
params_GPAW['setups']      = {'Fe': ':d,4.3'}              #U=4.3 applied to d orbitals

calc = GPAW(**params_GPAW)
fepo4.calc = calc
epot_fepo4_cell=fepo4.get_potential_energy()
print('E_Pot=', epot_fepo4_cell)

write('fepo4_out.traj', fepo4)

ens = BEEFEnsemble(calc)
dE = ens.get_ensemble_energies(2000)

with paropen('ensemble_fepo4.dat', 'a') as result:
    for e in dE:
        print(e, file=result)

# %%
"""
Uncomment the `%%writefile` line and execute the cell again and submit the calculation to the HPC cluster. The calculation should take around 10 minutes.
"""

# %%
# magic: !mq submit fepo4.py -R 8:1h  # submits the calculation to 8 cores, 1 hour

# %%
"""
Run the below cell to examine the status of your calculation. If no output is returned, the calculation has either finished or failed.
"""

# %%
# magic: !mq ls

# %%
"""
Once the calculation begins, you can run the cells below to open the error log and output of the calculation in a new window. This can be done while the calculation is running.
"""

# %%
# Error log
# magic: !cat "$(ls -t fepo4.py.*err | head -1)"

# %%
# Output
# magic: !cat "$(ls -t fepo4.py.*out | head -1)"

# %%
"""
Once the calculation has finished, load in the result. You can skip past this cell and return later.
"""

# %%
try:
    fepo4 = read('fepo4_out.traj')
    print('Calculation finished')
except FileNotFoundError:
    print('Calculation has not yet finished')

# %%
"""
## LiFePO$_4$


"""

# %%
"""
You will now do similar for LiFePO$_4$. In this case you will load in a template structure called `lifepo4_wo_li.traj` missing only the Li atoms. It is located in the resources folder.
"""

# %%
lifepo4_wo_li = read('lifepo4_wo_li.traj')

# %%
"""
Visualize the structure.
"""

# %%
view(lifepo4_wo_li)

# %%
"""
You should now add Li into the structure using the fractional coordinates below:
"""

# %%
# Li  0    0    0
# Li  0    0.5  0
# Li  0.5  0.5  0.5
# Li  0.5  0    0.5


# %%
"""
Add Li atoms into the structure, e.g., by following the example in [this ASE tutorial](https://wiki.fysik.dtu.dk/ase/gettingstarted/manipulating_atoms/manipulating_atoms.html?highlight=set_cell#manipulating-atoms).
"""

# %%
from numpy import identity
from ase import Atom

cell = lifepo4_wo_li.get_cell()

# ...

# lifepo4 = lifepo4_wo_li.copy()

# Teacher:
from numpy import identity
from ase import Atom

lifepo4 = lifepo4_wo_li.copy()
cell = lifepo4.get_cell()
xyzcell = identity(3)
lifepo4.set_cell(xyzcell, scale_atoms=True)  # Set the unit cell and rescale
lifepo4.append(Atom('Li', (0, 0, 0)))
lifepo4.append(Atom('Li', (0, 0.5, 0)))
lifepo4.append(Atom('Li', (0.5, 0.5, 0.5)))
lifepo4.append(Atom('Li', (0.5, 0, 0.5)))
lifepo4.set_cell(cell, scale_atoms=True)

# %%
"""
Visualize the structure with added Li.
"""

# %%
view(lifepo4)

# %%
"""
Ensure that the magnetic moments are as they should be, once again assuming ferromagnetism for simplicity.
"""

# %%
# ...

# teacher
print(lifepo4.get_initial_magnetic_moments())

# %%
"""
At this point you should save your structure by writing it to a trajectory file.
"""

# %%
write('lifepo4.traj', lifepo4)

# %%
"""
You should now calculate the potential energy of this sytem using the method and same calculational parameters as for FePO$_4$ above. Make a full script in the cell below similar to what you did above for FePO$_4$ and make sure that it runs.
"""

# %%
# %% writefile 'lifepo4.py'
from ase.parallel import paropen
from ase.io import read, write
from ase.dft.bee import BEEFEnsemble
from gpaw import GPAW, FermiDirac, Mixer, PW

# Read in the structure you made and wrote to file above
lifepo4 = read('lifepo4.traj')

params_GPAW = {...}

# ...
# ...
# ...

# write('lifepo4_out.traj', lifepo4)

# teacher
from ase.parallel import paropen
from ase.io import read
from ase.io.trajectory import Trajectory
from ase.dft.bee import BEEFEnsemble
from gpaw import GPAW, FermiDirac, Mixer, PW

#Read in the structure you made and wrote to file above
lifepo4 = read('lifepo4.traj')

params_GPAW = {}
params_GPAW['mode']        = PW(500)                     #The used plane wave energy cutoff
params_GPAW['nbands']      = -40                           #The number on empty bands had the system been spin-paired
params_GPAW['kpts']        = {'size':  (2,4,5),            #The k-point mesh
                              'gamma': True}
params_GPAW['spinpol']     = True                          #Performing spin polarized calculations
params_GPAW['xc']          = 'BEEF-vdW'                    #The used exchange-correlation functional
params_GPAW['occupations'] = FermiDirac(width = 0.1,      #The smearing
                                        fixmagmom = True)  #Total magnetic moment fixed to the initial value
params_GPAW['convergence'] = {'eigenstates': 1.0e-4,       #eV^2 / electron
                              'energy':      2.0e-4,       #eV / electron
                              'density':     1.0e-3,}
params_GPAW['mixer']       = Mixer(0.1, 5, weight=100.0)   #The mixer used during SCF optimization
params_GPAW['setups']      = {'Fe': ':d,4.3'}              #U=4.3 applied to d orbitals

calc = GPAW(**params_GPAW)
lifepo4.calc = calc
epot_lifepo4_cell=lifepo4.get_potential_energy()
print('E_Pot=', epot_lifepo4_cell)

traj=Trajectory('lifepo4_out.traj', mode='w', atoms=lifepo4)
traj.write()

ens = BEEFEnsemble(calc)
dE = ens.get_ensemble_energies(2000)
result = paropen('ensemble_lifepo4.dat','a')
for e in dE:
    print(e, file=result)
result.close()

# %%
"""
If the code runs, submit to the HPC cluster as you did above. The calculation takes approximately 10 minutes.
"""

# %%
# magic: !mq submit lifepo4.py -R 8:1h  # submits the calculation to 8 cores, 1 hour

# %%
"""
Run the below cell to examine the status of your calculation. If no output is returned, the calculation has either finished or failed.
"""

# %%
# magic: !mq ls

# %%
"""
Once the calculation begins, you can run the cells below to open the error log and output of the calculation in a new window.
"""

# %%
# Error log
# magic: !cat "$(ls -t lifepo4.py.*err | head -1)"

# %%
# Output
# magic: !cat "$(ls -t lifepo4.py.*out | head -1)"

# %%
"""
When calculation has finished, load in the result. You can skip past this cell and return later.
"""

# %%
try:
    lifepo4=read('lifepo4_out.traj')
    print('Calculation finished')
except FileNotFoundError:
    print('Calculation has not yet finished')

# %%
"""
### Li metal

We use a Li metal reference to calculate the equilibrium potential. On exercise day 2 you used a Li metal reference to calculate the intercalation energy in the graphite anode. The approach is similar here. Just read in the result of the calculation with DFTD3 and attach a new calulator to it.
"""

# %%
from ase import Atoms
from gpaw import GPAW, FermiDirac, PW
from ase.io import read

li_metal = read('Li-metal-DFTD3.traj')  # Change file name accordingly

calc = GPAW(mode=PW(500),
            kpts=(8, 8, 8),
            occupations=FermiDirac(0.15),
            nbands=-10,
            txt=None,
            xc='BEEF-vdW')

li_metal.calc = calc
li_metal.get_potential_energy()
li_metal.write('li_metal.traj')

# %%
"""
Now calculate the ensemble in the same way as for FePO$_4$ and LiFePO$_4$.
"""

# %%
ens = BEEFEnsemble(calc)
li_metal_ens_cell= ens.get_ensemble_energies(2000)
with paropen('ensemble_li_metal.dat', 'a') as result:
    for e in li_metal_ens_cell:
        print(e, file=result)

# %%
"""
## Calculate equilibrium potential and uncertainty
"""

# %%
"""
You can now calculate the equilibrium potential for the case of a FePO$_4$/Li metal battery from the intercallation energy of Li in FePO$_4$. For simplicity, use that assumption that all vibrational energies and entropic terms cancel each other. You should now have completed all submitted calculations before you proceed.
"""

# %%
"""
The calculated energies are for the full cells. Convert them to the energy per formula unit. The [len(...)](https://docs.python.org/3.6/library/functions.html#len) function can be quite helpful for this.
"""

# %%
epot_fepo4_cell=fepo4.get_potential_energy()
epot_lifepo4_cell=lifepo4.get_potential_energy()
epot_li_metal_cell=li_metal.get_potential_energy()
print('epot_fepo4_cell =', epot_fepo4_cell)
print('epot_lifepo4_cell =', epot_lifepo4_cell)
print('epot_li_metal_cell =', epot_li_metal_cell)

# %%
epot_fepo4 = epot_fepo4_cell / len(fepo4) * 6  # student: epot_fepo4 = ...
epot_lifepo4 = epot_lifepo4_cell / len(lifepo4) * 7  # student: epot_lifepo4 = ...
epot_li_metal = epot_li_metal_cell / len(li_metal)  # student: epot_li_metal = ...
# print(epot_fepo4, ...)

# %%
"""
No calculate the equilibrium potential under the assumption that it is given by $V_{eq} = \Delta U /e $, where $U$ is the electronic potential energy of the system and $e$ is the number of electrons transfered.
"""

# %%
# V_eq = ...

# teacher
V_eq = epot_lifepo4 - epot_fepo4 - epot_li_metal
print(V_eq)

# %%
"""
You will now calculate the error estimate for the Li intercallation energy in FePO$_4$ using the BEEF ensemble results. Start by loading in the files. Wait a few minutes and rerun the cell if the number is not 2000 for all of them.
"""

# %%
fepo4_ens_cell = np.genfromtxt('ensemble_fepo4.dat', max_rows=2000)
lifepo4_ens_cell = np.genfromtxt('ensemble_lifepo4.dat', max_rows=2000)

print('number of functionals in ensemble=', len(fepo4_ens_cell))
print('number of functionals in ensemble=', len(lifepo4_ens_cell))
print('number of functionals in ensemble=', len(li_metal_ens_cell))

# %%
"""
Note that these are energies per cell and not per formula unit. Convert them as you did the potential energies above. Note that you are now performing the operation on a list of length 2000 and not a single float value as before.
"""

# %%
# fepo4_ens = fepo4_ens_cell / ...
# ...
# ...

# teacher
fepo4_ens = fepo4_ens_cell / len(fepo4) * 6
lifepo4_ens = lifepo4_ens_cell / len(lifepo4) * 7
li_metal_ens = li_metal_ens_cell / len(li_metal)

# %%
"""
Make a list of equilibrium potentials.
"""

# %%
# V_eq_ens = lifepo4_ens - ...

# teacher
V_eq_ens = lifepo4_ens - fepo4_ens - li_metal_ens

# %%
"""
Use the plot command below to visualize the distribution.
"""

# %%
plt.hist(V_eq_ens, 50)
plt.grid(True)

# %%
"""
Use the [NumPy function standard deviation function](https://docs.scipy.org/doc/numpy/reference/generated/numpy.std.html#numpy-std) to obtain the standard deviation of the ensemble.
"""

# %%
# error = ...
# print(error)

# teacher
error = np.std(V_eq_ens)
print(error)

# %%
"""
The equilibrium potential for a FePO$_4$/Li battery is thus as a good estimate:
"""

# %%
print(f'{V_eq:.2f} V +- {error:.2f} V')

# %%
"""
You can get the equilibrium potential for the FePO$_4$/C battery using the intercallation energy of Li in graphite, that you calculated on Day 2. What equilibrium potential do you find? How does that compare to the cell voltage you can obtain from FePO$_4$/C batteries?
"""

# %%
# You can use this cell for FePO4/C potential calculation

# %%
"""
Make sure you are able to discuss the bullet points at the top of this notebook. You can use the cell below for calculations.
"""

# %%


# %%
"""
## Bonus
"""

# %%
"""
How does the predicted error estimate change if you consider the full reaction from Li in graphite + FePO4  to empty graphite + LiFePO4.
"""
