import numpy as np
import sys
import os
import time
from ase.units import Bohr
from ase.lattice import bulk
from gpaw import GPAW, PW
from gpaw.test import findpeak
from gpaw.atom.basis import BasisMaker
from gpaw.response.df import DF
from gpaw.mpi import serial_comm, rank, size
from gpaw.utilities import devnull


if rank != 0:
  sys.stdout = devnull 

assert size <= 4**3

# Ground state calculation

t1 = time.time()

a = 4.043
atoms = bulk('Al', 'fcc', a=a)
atoms.center()
calc = GPAW(mode=PW(200),
            kpts=(4,4,4),
            parallel={'band':1},
            idiotproof=False,  # allow uneven distribution of k-points
            xc='LDA')

atoms.set_calculator(calc)
atoms.get_potential_energy()
t2 = time.time()

# Excited state calculation
q = np.array([1/4.,0.,0.])
w = np.linspace(0, 24, 241)

df = DF(calc=calc, q=q, w=w, eta=0.2, ecut=(50,50,50))
df.get_EELS_spectrum(filename='EELS_Al')
df.check_sum_rule()
df.write('Al.pckl')

t3 = time.time()

print 'For ground  state calc, it took', (t2 - t1) / 60, 'minutes'
print 'For excited state calc, it took', (t3 - t2) / 60, 'minutes'

d = np.loadtxt('EELS_Al')

# New results are compared with test values

wpeak1,Ipeak1 = findpeak(d[:,0],d[:,1])
wpeak2,Ipeak2 = findpeak(d[:,0],d[:,2])


test_wpeak1 = 15.70 # eV
test_Ipeak1 = 28.90 # eV
test_wpeak2 = 15.725 # eV
test_Ipeak2 = 26.24 # eV


if np.abs(test_wpeak1-wpeak1)<1e-2 and np.abs(test_wpeak2-wpeak2)<1e-2:
    pass
else:
    print test_wpeak1-wpeak1,test_wpeak2-wpeak2
    raise ValueError('Plasmon peak not correct ! ')

if np.abs(test_Ipeak1-Ipeak1)>1e-2 or np.abs(test_Ipeak2-Ipeak2)>1e-2:
    print Ipeak1, Ipeak2
    raise ValueError('Please check spectrum strength ! ')






