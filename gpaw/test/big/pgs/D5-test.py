from gpaw.pgs import GPAWULMSymmetryCalculator
from gpaw.pgs import tools

import gpaw.mpi
from gpaw import GPAW
from ase import Atom, Atoms

from gpaw.test import equal
import math
import numpy as np

# Build molecule:
name = 'ferrocene-chiral'
system = Atoms()
system.append(Atom('Fe', position=[0,0,0]))

height = 1.68 #distance of cyclopentadienyl from Fe
CC = 1.425 #C-C distance
CH = 1.080 #C-H distance
FeC_xy = CC / (2. * math.sin(2 * math.pi / 5. / 2.)) #distance of C from z axis
chirality_angle = 14. * (2 * math.pi / 360.)

# Top ring:
for i in range(5):
    x = FeC_xy * math.cos(i * 2 * math.pi / 5. + chirality_angle / 2.)
    y = FeC_xy * math.sin(i * 2 * math.pi / 5. + chirality_angle / 2.)
    z = height
    system.append(Atom('C', position=[x, y, z]))

    x = (FeC_xy + CH) * math.cos(i * 2 * math.pi / 5. + chirality_angle / 2.)
    y = (FeC_xy + CH) * math.sin(i * 2 * math.pi / 5. + chirality_angle / 2.)
    z = height
    system.append(Atom('H', position=[x, y, z]))

# Bottom ring:
for i in range(5):
    x = FeC_xy * math.cos(i * 2 * math.pi / 5. - chirality_angle / 2.)
    y = FeC_xy * math.sin(i * 2 * math.pi / 5. - chirality_angle / 2.)
    z = -height
    system.append(Atom('C', position=[x, y, z]))

    x = (FeC_xy + CH) * math.cos(i * 2 * math.pi / 5. - chirality_angle / 2.)
    y = (FeC_xy + CH) * math.sin(i * 2 * math.pi / 5. - chirality_angle / 2.)
    z = -height
    system.append(Atom('H', position=[x, y, z]))


system.center(vacuum=5.0)

h = 0.2
calc = GPAW(h=h,
            nbands=50,
            txt='%s-gpaw.txt' % name
            )
system.set_calculator(calc)
e = system.get_potential_energy()

calc.write('%s.gpw' % name, mode='all')


# Symmetry analysis:

symcalc = GPAWULMSymmetryCalculator(filename='%s.gpw' % name,
                                    statelist=range(28),
                                    pointgroup='D5',
                                    mpi=gpaw.mpi,
                                    overlapfile='overlaps_%s.txt' % name,
                                    symmetryfile='symmetries_%s.txt' % name)

symcalc.initialize()

# Define atom indices around which the analysis is run:
coreatoms = range(len(system))

# Deliver the required rotations:
Rx, Ry, Rz =  (0.,0.,0.)
symcalc.set_initialrotations(Rx, Ry, Rz)

# Determine some parameters from the data:
wfshape = tools.get_wfshape(symcalc)

# Give the grid spacing to the symmetry calculator for shifting
# the atoms to center:
h = tools.get_h(symcalc)
symcalc.set_gridspacing(h)

# Set up the volume where the analysis is restricted:
symcalc.set_cutarea(tools.calculate_cutarea(atoms=symcalc.atoms,
                                            coreatoms=coreatoms,
                                            wfshape=wfshape,
                                            gridspacing=h,
                                            cutlimit=3.00))

# Set up the shift vector based on the center-of-mass of `coreatoms`:
symcalc.set_shiftvector(tools.calculate_shiftvector(atoms=symcalc.atoms,
                                                coreatoms=coreatoms,
                                                gridspacing=h))

# Calculate the symmetry representation weights of the wave functions:
symcalc.calculate(analyze=True)

if gpaw.mpi.rank == 0:
    f = open('symmetries_%s.txt' % name, 'r')
    results = []
    for line in f:
        if line.startswith('#'):
            continue
        results.append(line.split()[:-1])
    f.close()

    results = np.array(results).astype(float)
    for i in range(len(results)):
        norm = results[i, 2]
        bestweight = (results[i, 3:]).max()
        equal(bestweight / norm, 1.0, 0.1)

