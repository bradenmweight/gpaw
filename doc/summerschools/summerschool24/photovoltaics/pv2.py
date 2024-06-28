# %%
"""
# Calculating absorption spectra
"""

# %%
"""
In this exercise we are going to calculate the absorption spectra of the material you have considered so far.
In general, the absorption spectrum is given by the imaginary part of the macroscopic dielectric function. This means that our computational task is to calculate the macroscopic dielectric function for our material, and then plot the imaginary part. For more information about the dielectric function please consult:

https://wiki.fysik.dtu.dk/gpaw/tutorialsexercises/opticalresponse/dielectric_response/dielectric_response.html

To calculate the dielectric function, we will use the Random Phase Approximation (RPA) correlation energy (so we say we calculate the absorption spectrum within the random phase approximation).
(Details about RPA to calculate the total energy can be found here: https://wiki.fysik.dtu.dk/gpaw/documentation/xc/rpa.html#rpa)

As discussed earlier, it is of greatest importance to make sure our calculations are converged. Since the convergence parameters is not necessarily the same for band gaps and RPA absorption spectra (and any other material property) with respect to for instance plane-wave cut-off or k-points, we need to do do a new ground state calculation with parameters. For RPA absorption spectra we will here look at the number of bands included in the calculation and the k-point mesh. We will therefore restart from the previous ground state file, and update some of the parameters. First we will look at a too rough k-point mesh and then do more calculations with a finer k-point mesh. In this way we can see the importance of converging the calculations.

The new ground state is calculated below. For BN use (24,24,1) k-points while for the bulk materials use (12,12,4) k-points.
"""

# %%
from gpaw import GPAW
from gpaw.mpi import world
from gpaw.occupations import FermiDirac


#Load and initialize ground state gpw file from previous exercise
calc_old = GPAW('CdTe_gs_LDA.gpw', txt=None) # student: calc_old = GPAW('.gpw', txt=None)

#Extract number of valence bands:
nval = calc_old.wfs.nvalence

# Do new ground state calculations with more k-points.
# This is because in general RPA calculations requires more k-poins to be converged.

calc = GPAW('CdTe_gs_LDA.gpw').fixed_density(  # student: calc = GPAW('???.gpw').fixed_density(
    kpts=(12, 12, 4),  # student: kpts=???,
    nbands=8 * nval,  # number of bands to include in calculation
    convergence={'bands': 6 * nval},  # number of bands to convergence
    txt='es.txt',  # student: txt = '???',
    occupations=FermiDirac(width=1e-4))


calc.get_potential_energy()
#Now save the .gpw file. 'all' means we save all the wave functions to the gpw file. This is required for the rpa calculations
calc.write('CdTe_12x12x4.gpw', 'all') # student: calc.write('', 'all')

# %%
"""
Now that we have obtained the ground state, it is time to calculate the dielectric function. This is done below where we first initialize the parameters needed for an RPA calculation, then calculate the dielectric function, and finally obtain the polarizability. Note that one also would have to converge the parameters initialized below, for a fully converged study, however in the interest of time and computational power we will not consider this here. The principle is however exactly the same, as you will see with the different k-point meshes we will employ here.

These calculations are both time-wise and memory-wise heavier than what you previously encountered. Therefore we will need to submit these calculations to the databar so they can run over night. Open a new SSH terminal, edit and copy the below code into a script format (.py file), and submit the calculations from the terminal using the following command:

mq submit -R 8:15h script.py

This will submit the script with the name "script.py" to 8 cores with a maximum time of 15 hours.
"""

# %%
#Define parameters for rpa calculations. You should change the name of the output file so it corresponds with your calculation:
from gpaw import GPAW
from gpaw.mpi import world
from gpaw.occupations import FermiDirac
from gpaw.response.df import DielectricFunction

#Insert the name of your structure and the k-point grid you are using instead of "???".

calc_old = GPAW('CdTe_gs_LDA.gpw', txt=None) # student: calc_old = GPAW('.gpw', txt=None)
nval = calc_old.wfs.nvalence

#Note: For a 2D material (here BN) we use an additional keyword in the parameters below: 'truncation': '2D'
#This truncates the Coulomb interaction in the lateral direction to avoid non-physical periodic interactions.

kwargs = {
    'frequencies': {
        'type': 'nonlinear',
        'domega0': 0.01},  # Define spacing of frequency grid
    'eta': 0.05,           # Broadening parameter
    'intraband': False,    # Here we do not include intraband transitions for calculating the absorption spectrum
    'nblocks': 8,          # Number of blocks used for parallelization
    'ecut': 50,            # Plane wave cutoff in eV
    'nbands': 3 * nval}    # Number of bands included in rpa calculation

#Calculate dielectric function. Takes ground state calculation and defined parameters in "kwargs" as input:
df = DielectricFunction('CdTe_12x12x4.gpw', **kwargs) # student: df = DielectricFunction('.gpw', **kwargs)


#Finally we calculate he polarizability in the x, y, and z direction. The output is a .csv file (one for each direction) which can be plotted.
#Consider the symmetries of your material and figure out if calculating the absorption spectrum in all 3 directions
#is needed or two or more of them will be the same. In that case only do the directions you will need (this is for saving
#time and memory). Also remember to change the filename so it fits with your calculation.

df.get_polarizability(xc='RPA',                         #We want to calculate the absorption spectrum within RPA
                      q_c = [0, 0, 0],                  #We consider the zero momentum wave vector
                      direction = 'x',                  #Define real space direction
                      filename='CdTe_rpa_x.csv')  # student: filename='=???_rpa_x.csv'       #Name of output file

df.get_polarizability(xc='RPA',
                      q_c = [0, 0, 0],
                      direction = 'y',
                      filename='CdTe_rpa_y.csv')  # student: filename='=???_rpa_y.csv'
df.get_polarizability(xc='RPA',
                      q_c = [0, 0, 0],
                      direction = 'z',
                      filename='CdTe_rpa_z.csv') # student: filename='=???_rpa_z.csv'

# %%
"""
Now it is time to do more calculations with finer k-point mesh to see the importance of converging our calculations. Repeat the above calculations with the following parameters:

For BN try: (24,24,1), (40,40,1), and (60,60,1) k-points
For bulk materials try e.g.: (12,12,4), (18,18,6), and (24,24,8) k-points
(Remember to use a meshes according to the length of lattice vectors.)

For the calculations with more k-points you probably want to combine the (new) ground state calculation and the calculation of the absorption spectrum into one script which you will submit over night.

"""
