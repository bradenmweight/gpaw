# %%
# teacher
import ase.visualize as viz
viz.view = lambda atoms, repeat=None: None
# %%
"""
# Battery Project

## Day 4 - Transport barriers and Voltage profile
"""

# %%
"""
Today you will calculate the energy barriers for transport of Li intercalated in the graphite anode. You will examine how sensitive this barrier is to the interlayer distance in graphite.  You will also examine the energy of intermediate states during the charge/discharge process. This will allow some basic discussion of the voltage profile of the battery.

You will in general be provided less code than yesterday, especially towards the end of this notebook. You will have to use what you have already seen and learned so far.

There will be some natural pauses while you wait for calculations to finish. If you do not finish this entire notebook today, do not despair.
"""

# %%
"""
## Initialize
"""

# %%
# magic: %matplotlib notebook
from ase import Atom
from ase.visualize import view
import matplotlib.pyplot as plt
from ase.io import read, write
from ase.mep import NEB
from ase.optimize import BFGS
from ase.parallel import paropen
from gpaw import GPAW, FermiDirac, Mixer, PW
from ase.constraints import FixAtoms

# %%
"""
## Transport barrier of Li in graphite
"""
# %%
"""
You will now calculate the energy barrier for Li diffusion in the graphite anode. You will do this using the [Nudged Elastic Band (NEB) method](https://wiki.fysik.dtu.dk/ase/ase/neb.html#module-ase.neb)

You can use your work from Day 2, but for simplicity you are advised to load in the initial atomic configuration from file.
"""

# %%
initial = read('NEB_init.traj')

# %%
"""
Visualize the structure.
"""

# %%
view(initial)

# %%
"""
You will now make a final structure, where the Li atom has been moved to a neighbouring equivalent site. The [`get_positions`](https://wiki.fysik.dtu.dk/ase/ase/atoms.html?highlight=get_positions#ase.Atoms.get_positions), [`set_positions`](https://wiki.fysik.dtu.dk/ase/ase/atoms.html?highlight=get_positions#ase.Atoms.set_positions) and [`get_cell`](https://wiki.fysik.dtu.dk/ase/ase/atoms.html?highlight=get_positions#ase.Atoms.get_cell) functions are highly useful for such a task. HINT: Displace the Li atom $\frac{1}{n} (\vec{a}+\vec{b})$
"""

# %%
final = initial.copy()

# %%
# ...
# ...

# teacher
cell = final.get_cell()
pos = final.get_positions()
pos[6] = pos[6] + cell[1] / 3 + cell[0] / 3
final.set_positions(pos)

# %%
"""
Visualize that you have made the final strcuture correctly.
"""

# %%
view(final)

# %%
"""
Make a band consisting of 7 images including the initial and final.
"""

# %%
images = [initial]
images += [initial.copy() for i in range(5)]  # These will become the minimum energy path images.
images += [final]

# %%
"""
It this point `images` consist of 6 copies of `initial` and one entry of `final`. Use the `NEB` method to create an initial guess for the minimum energy path (MEP). In the cell below a simple interpolation between the `initial` and `final` image is used as initial guess.
"""

# %%
neb = NEB(images)
neb.interpolate()

# %%
"""
Visualize the NEB images.
"""

# %%
view(images)

# %%
"""
It turns out, that while running the NEB calculation, the largest amount of resources will be spend translating the carbon layer without any noticeable buckling. You will thus [constrain](https://wiki.fysik.dtu.dk/ase/ase/constraints.html#constraints) the positions of the carbon atoms to save computational time.

Each image in the NEB requires a unique calculator.

This very simple case is highly symmetric. To better illustrate how the NEB method works, the symmetry is broken using the [rattle](https://wiki.fysik.dtu.dk/ase/ase/atoms.html#ase.Atoms.rattle) function.
"""

# %%
for image in images[0:7]:
    calc = GPAW(mode=PW(500), kpts=(5, 5, 6), xc='LDA', txt=None, symmetry={'point_group': False})
    image.calc = calc
    image.set_constraint(FixAtoms(mask=[atom.symbol == 'C' for atom in image]))

images[3].rattle(stdev=0.05, seed=42)


# %%
"""
Start by calculating the energy and forces of the first (`initial`) and last (`final`) images as this is not done during the actual NEB calculation.

Note, that this can take a while if you opt to do it inside the notebook.
"""

# %%
images[0].get_potential_energy()
images[0].get_forces()
images[6].get_potential_energy()
images[6].get_forces()


# %%
"""
You can run the NEB calculation by running an optimization on the NEB object the same way you would on an atoms object. Note the `fmax` is larger for this tutorial example than you would normally use.
"""

# %%
optimizer = BFGS(neb, trajectory='neb.traj', logfile='neb.log' )
optimizer.run(fmax=0.10)


# %%
"""
Submit the calculation to the HPC cluster. Do this by first building a complete script in the cell below using the cells above (minus the `view()` commands). Make sure the cell runs and then interrupt the kernel.
"""

# %%
#from ase.io import read, write
#from ase.mep import NEB
#from ase.optimize import BFGS
#from ase.parallel import paropen
#from gpaw import GPAW, FermiDirac, Mixer, PW
#from ase.constraints import FixAtoms

# initial = read('NEB_init.traj')

# final = ...

# ...
# ...

# optimizer.run(fmax=0.10)

# teacher
from ase.io import read, write
from ase.mep import NEB
from ase.optimize import BFGS
from ase.parallel import paropen
from gpaw import GPAW, FermiDirac, Mixer, PW
from ase.constraints import FixAtoms

initial=read('NEB_init.traj')

final=initial.copy()
cell=final.get_cell()
pos=final.get_positions()
pos[6]=pos[6]+cell[1]/3.+cell[0]/3.
final.set_positions(pos)

images = [initial]
images += [initial.copy() for i in range(5)]  #These will become the minimum energy path images.
images += [final]

neb = NEB(images)
neb.interpolate()

for image in images[0:7]:
    calc = GPAW(mode=PW(500), kpts=(5, 5, 6), xc='LDA', symmetry={'point_group': False})
    image.calc = calc
    image.set_constraint(FixAtoms(mask=[atom.symbol == 'C' for atom in image]))

images[3].rattle(stdev=0.05, seed=42)

images[0].get_potential_energy()
images[0].get_forces()
images[6].get_potential_energy()
images[6].get_forces()

optimizer = BFGS(neb, trajectory = 'neb.traj', logfile = 'neb.log' )
optimizer.run(fmax = 0.10)

# %%
"""
You can use the cell below to submit the calculation in the same way as on earlier days.
"""

# %%
# magic: !mq submit NEB.py -R 8:1h  # submits the calculation to 8 cores, 1 hour

# %%
"""
Run the below cell to examine the status of your calculation.
"""

# %%
# magic: !mq ls

# %%
"""
You can run the cells below to open the error log and output of the calculation in a new window. This can be done while the calculation is running.
"""

# %%
# magic: !cat "$(ls -t NEB.py.*err | head -1)"

# %%
# magic: !cat "$(ls -t NEB.py.*out | head -1)"

# %%
"""
The optimiziation progress can be seen by running the below cell.
"""

# %%
# magic: !cat neb.log

# %%
"""
You can move on while you wait for the calculation to finish.

Once the maximum force (`fmax`) in the log is below 0.1, the calculation is finished.
Load in the full trajectory.
"""

# %%
full_neb = read('neb.traj@:')

# %%
"""
You will use the `ase gui` to inspect the result. The below line reads in the last 7 images in the file. In this case the MEP images.
"""

# %%
# magic: !ase gui neb.traj@-7:

# %%
"""
In the GUI use `Tools` $\rightarrow$ `NEB`.

Now inspect how the TS image has developed.
"""

# %%
# magic: !ase gui neb.traj@3::7

# %%
"""
For more complicated MEP's, use the [climbing image method](https://wiki.fysik.dtu.dk/ase/ase/neb.html?highlight=neb#climbing-image) to determine the transition state. Why is it not required here?
"""

# %%
"""
## Bonus

You will now study the influence of changing the interlayer graphite distance on the energy barrier. Due to the high degree of symmetry, this can be done easily in this case. Load in the initial state (IS) and transition state (TS) images from the converged MEP.
"""

# %%
IS_image = images[0]
TS_image = images[3]

# %%
"""
Now calculate the energy of the initial state (IS) image and the transition state (TS) image using [`get_potential_energy()`](https://wiki.fysik.dtu.dk/ase/ase/atoms.html?highlight=get_potential_energy#ase.Atoms.get_potential_energy)
"""

# %%
epot_IS = IS_image.get_potential_energy()
#epot_TS= ...

# teacher
epot_TS = TS_image.get_potential_energy()

# %%
barrier = epot_TS - epot_IS
print('Energy barrier:', barrier)

# %%
"""
Why does this not fully align with what you found before? New reduce the graphite layer distance by change the the size of the unit cell in the *z* direction by 3 %.
"""

# %%
cell = IS_image.get_cell()
IS_image97 = IS_image.copy()
IS_image97.set_cell([cell[0], cell[1], cell[2] * 0.97], scale_atoms=True)
TS_image97 = TS_image.copy()
TS_image97.set_cell([cell[0], cell[1], cell[2] * 0.97], scale_atoms=True)

# %%
"""
Use the same calculator object as you did above and calculate the potential energy of the compressed initial and final state.
"""

# %%
# calc = ...
# ...

# teacher
calc = GPAW(mode=PW(500), kpts=(5, 5, 6), xc='LDA', symmetry={'point_group': False})
TS_image97.calc = calc
calc = GPAW(mode=PW(500), kpts=(5, 5, 6), xc='LDA', symmetry={'point_group': False})
IS_image97.calc = calc

# %%
"""
Now calculate the energy of the compressed IS and TS.
"""

# %%
# epot_TS97 = ...

# teacher
epot_TS97 = TS_image97.get_potential_energy()
epot_IS97 = IS_image97.get_potential_energy()

# %%
"""
What is the energy barrier now?
"""

# %%
# barrier97 = ...
# print('Energy barrier:', barrier97)

# teacher
barrier97=epot_TS97-epot_IS97
print("Energy barrier:", barrier97)

# %%
"""
Now repeat the procedure but expanding the intergraphite distance by 3 %.
"""

# %%
# IS_image103 = IS_image.copy()
# IS_image103.set_cell(...

# calc ...


# epot_TS103 = ...
# ...

# teacher
IS_image103=IS_image.copy()
IS_image103.set_cell([cell[0],cell[1],cell[2]*1.03], scale_atoms=True)
TS_image103=TS_image.copy()
TS_image103.set_cell([cell[0],cell[1],cell[2]*1.03], scale_atoms=True)

calc = GPAW(mode=PW(500), kpts=(5, 5, 6), xc='LDA', symmetry={'point_group': False})
TS_image103.calc = calc
calc = GPAW(mode=PW(500), kpts=(5, 5, 6), xc='LDA', symmetry={'point_group': False})
IS_image103.calc = calc

epot_TS103=TS_image103.get_potential_energy()
epot_IS103=IS_image103.get_potential_energy()

# %%
"""
What is the energy barrier now?
"""

# %%
# barrier103 = ...
# print('Energy barrier:', barrier103)

# teacher
barrier103 = epot_TS103 - epot_IS103
print('Energy barrier:', barrier103)

# %%
"""
## FePO$_4$ with one Li
"""

# %%
"""
You will now calculate the energy gain of adding a single Li atom into the FePO$_4$ cell you made on Day 3. This corresponds to a charge of 25 %. You can compare this energy to the equilibrium potential.

Load in the FePO$_4$ structure you wrote to file on in a previous exercise and add Li. Assume that the cell dimension remain unchanged.
"""

# %%
#fepo4=read('fepo4.traj')
#fepo4_1li=fepo4.copy()

# teacher
from ase import Atoms
fepo4=Atoms('FeFeFeFeOOOOOOOOOOOOOOOOPPPP',
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

for atom in fepo4:
    if atom.symbol == 'Fe':
        atom.magmom = 5.0

fepo4_1li = fepo4.copy()

# %%
# fepo4_1li.append(...)

# teacher
fepo4_1li.append('Li')

# %%
"""
Visualize the structure you made.
"""

# %%
view(fepo4_1li)

# %%
"""
Adjust the total magnetic moment of the cell such that it is 19.
"""

# %%
for atom in fepo4_1li:
    if atom.symbol == 'Fe':
        atom.magmom = 4.75

print(sum(fepo4_1li.get_initial_magnetic_moments()))

# %%
"""
Write your atoms object to file giving it the name `fepo4_1li.traj`.
"""

# %%
write('fepo4_1li.traj', fepo4_1li)

# %%
"""
Make a full script in the cell below similar to those you made yesterday. Make sure the cell runs before interupting the notebook kernel.
"""

# %%
# %%writefile 'fepo4_1li.py'
#from ase.parallel import paropen
#from ase.io import read, write
#from ase.dft.bee import BEEFEnsemble
#from gpaw import GPAW, FermiDirac, Mixer, PW

# Read in the structure you made and wrote to file above
fepo4_1li = read('fepo4_1li.traj')

#...
#...

# write('fepo4_1li_out.traj', fepo4_1li)

# ens = BEEFEnsemble(calc)
# with paropen('ensemble_fepo4_1li.dat', 'a') as result:
#     for e in dE:
#         print(e, file=result)

# teacher
from ase.parallel import paropen
from ase.io import read
from ase.dft.bee import BEEFEnsemble
from gpaw import GPAW, FermiDirac, Mixer, PW

#Read in the structure you made and wrote to file above
fepo4_1li=read('fepo4_1li.traj')

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
fepo4_1li.calc = calc
epot_fepo4_1li_cell=fepo4_1li.get_potential_energy()
print('E_Pot=', epot_fepo4_1li_cell)

write('fepo4_1li_out.traj', fepo4_1li)

ens = BEEFEnsemble(calc)
dE = ens.get_ensemble_energies(2000)
result = paropen('ensemble_fepo4_1li.dat','a')
for i in range(0,len(dE)):
    print(dE[i], file=result)
result.close()

# %%
"""
Submit this calculation to the HPC cluster as you did on exercise day 3.
"""

# %%
# magic: !mq submit fepo4_1li.py -R 8:1h # submits the calculation to 8 cores, 1 hour

# %%
"""
Run the below cell to examine the status of your calculation.
"""

# %%
# magic: !mq ls

# %%
"""
You can run the cells below to open the error log and output of the calculation in a new window. This can be done while the calculation is running.
"""

# %%
# magic: !cat "$(ls -t fepo4_1li.py.*err | head -1)"

# %%
# magic: !cat "$(ls -t fepo4_1li.py.*out | head -1)"

# %%
"""
You can move on while you wait for the calculation to finish. Once the calculation is finished load in the structure by running the cell below.
"""

# %%
try:
    fepo4_1li=read('fepo4_1li_out.traj')
    print('Calculation finished')
except FileNotFoundError:
    print('Calculation has not yet finished')

# %%
"""
You are now ready to calculate the energy gained by intercalating a single Li ion into the cathode. Start by loading in the relevant reference structures and obtain the potential energies. This should not require any DFT calculations.
"""

# %%
# Loading in files from exercise day 3.
li_metal = read('li_metal.traj')
fepo4 = read('fepo4_out.traj')

epot_li_metal = li_metal.get_potential_energy() / len(li_metal)

# %%
# epot_fepo4 = ...
# ...

# teacher
epot_fepo4=fepo4.get_potential_energy()
epot_fepo4_1li=fepo4_1li.get_potential_energy()

# %%
"""
Calculate the energy of intercalting a single Li in the FePO$_4$ cell. How does this energy compare with the equilibirum potential? What can it tell you about the charge/discharge potential curves?
"""

# %%
# ...
# print(...)

# teacher
li_cost=epot_fepo4_1li-epot_fepo4-epot_li_metal
print(li_cost)

# %%
"""
## Bonus: LiFePO$_4$ with one vacancy
"""

# %%
"""
If time permits, you will now do a similar calculation but this time with LiFePO$_4$ contraining one vacancy. Once again you should assume that the cell dimension remain unchanged compaired to LiFePO$_4$.

There are numerous ways to obtain this structure. You can get inspiration from the way LiFePO$_4$ was made on Exercise day 3, use the [`del` or `pop()` methods](https://wiki.fysik.dtu.dk/ase/ase/atoms.html?highlight=pop#list-methods), or even use the GUI to delete an atom and save the structure afterwards.
"""

# %%
# In this cell you create the vacancy in LiFePO4

# lifepo4_vac = ...

# ...

# teacher
lifepo4_wo_li=read('lifepo4_wo_li.traj')
from numpy import identity
cell=lifepo4_wo_li.get_cell()
xyzcell = identity(3)
lifepo4_wo_li.set_cell(xyzcell, scale_atoms=True)  # Set the unit cell and rescale
#lifepo4_wo_li.append(Atom('Li', (0, 0, 0)))
lifepo4_wo_li.append(Atom('Li', (0, 0.5, 0)))
lifepo4_wo_li.append(Atom('Li', (0.5, 0.5, 0.5)))
lifepo4_wo_li.append(Atom('Li', (0.5, 0, 0.5)))
lifepo4_wo_li.set_cell(cell, scale_atoms=True)
lifepo4_vac=lifepo4_wo_li.copy()

# %%
"""
Visualize the structure
"""

# %%
view(lifepo4_vac)

# %%
"""
Now ensure that the total magnetic moment is equal to 17.
"""

# %%
for atom in fepo4_1li:
    if atom.symbol == 'Fe':
        atom.magmom = 4.25

print(sum(fepo4_1li.get_initial_magnetic_moments()))

# %%
"""
Write your atoms object to file giving it the name `lifepo4_vac.traj`.
"""

# %%
# ...

# teacher
write('lifepo4_vac.traj', lifepo4_vac)

# %%
"""
Make a full script in the cell below similar to that you made above. Make sure the cell runs before interupting the notebook kernel.
"""

# %%
# %%writefile 'lifepo4_vac.py'
# from ase.parallel import paropen
# from ase.io import read, write
# from ase.dft.bee import BEEFEnsemble
# from gpaw import GPAW, FermiDirac, Mixer, PW

# Read in the structure you made and wrote to file above
# lifepo4_vac = read('lifepo4_vac.traj')


# ...

# write('lifepo4_vac_out.traj', lifepo4_vac)

# ens = BEEFEnsemble(calc)
# dE = ens.get_ensemble_energies(2000)
# with paropen('ensemble_lifepo4_vac.dat','a') as results:
#     for e in dE:
#         print(e, file=result)

# teacher
from ase.parallel import paropen
from ase.io import read, write
from ase.dft.bee import BEEFEnsemble
from gpaw import GPAW, FermiDirac, Mixer, PW

# Read in the structure you made and wrote to file above
lifepo4_vac = read('lifepo4_vac.traj')

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
lifepo4_vac.calc = calc
epot_lifepo4_vac_cell=lifepo4_vac.get_potential_energy()
print('E_Pot=', epot_lifepo4_vac_cell)

write('lifepo4_vac_out.traj', lifepo4_vac)

ens = BEEFEnsemble(calc)
dE = ens.get_ensemble_energies(2000)
result = paropen('ensemble_lifepo4_vac.dat','a')
for i in range(0,len(dE)):
    print(dE[i], file=result)
result.close()

# %%
"""
Once you have made sure the cell runs, submit it to the HPC cluster.
"""

# %%
# magic: !mq submit lifepo4_vac.py -R 8:1h  # submits the calculation to 8 cores, 1 hour

# %%
"""
Once the calculation has finished, load in the trajectory.
"""

# %%
try:
    lifepo4_vac=read('lifepo4_vac_out.traj')
    print('Calculation finished')
except FileNotFoundError:
    print('Calculation has not yet finished')

# %%
"""
Once the calculation has finished you are ready to calculate the energy cost of creating a li vacancy in the fully lithiated LiFePO$_4$. Start by loading in the relevant reference structures and obtain the potential energies. This should not require any calculations.
"""

# %%
# Loading in files from exercise day 3.
li_metal = read('li_metal.traj')   # you should have already read this in above
lifepo4 = read('lifepo4_out.traj')

epot_li_metal = li_metal.get_potential_energy() / len(li_metal)

# %%
# epot_lifepo4 = ...
# ...

# teacher
epot_lifepo4=lifepo4.get_potential_energy()
epot_lifepo4_vac=lifepo4_vac.get_potential_energy()

# %%
# vac_cost = ...
# print(vac_cost)

# teacher
vac_cost=epot_lifepo4_vac-epot_lifepo4+epot_li_metal
print(vac_cost)

# %%
"""
How does this energy compare with the equilibirum potential? What can it tell you about the charge/discharge potential curves?
"""

# %%
"""
## Bonus
Calculate the error estimates of the energy for the added Li atom and vacancy formation using the ensembles.
"""

# %%
# Cell for bonus question

# %%
# Cell for bonus question

# %%
# Cell for bonus question
