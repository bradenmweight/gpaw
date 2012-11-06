from gpaw.tddft import *
from ase import Atoms
from gpaw import GPAW
from sys import argv
from gpaw.tddft.lcao_tddft import LCAOTDDFT
import gpaw.lrtddft as lrtddft
from sys import exit
from gpaw.mpi import world
N = 3
xc = 'oldLDA'
c = +1
h = 0.4
b = 'sz'
sy = 'Na'+str(N)
positions = []
for i in range(N):
    positions.append([0.00,0.00,i*2.50])
atoms = Atoms(symbols=sy, positions = positions)
atoms.center(vacuum=3)

# LCAO-RT-TDDFT
calc = LCAOTDDFT(mode='lcao', xc=xc, h=h, basis=b, 
                dtype=complex, charge=c, width=0, convergence={'density':1e-6}, 
                propagator_debug=True)#, propagator='taylor')
                #parallel={'sl_default':(1,1,64)})
atoms.set_calculator(calc)
atoms.get_potential_energy()
dmfile = sy+'_lcao_'+b+'_rt_z.dm'+str(world.size)
specfile = sy+'_lcao_'+b+'_rt_z.spectrum'+str(world.size)
calc.absorption_kick([0.0,0,0.001])
calc.propagate(10, 20, dmfile)
if world.rank == 0:
    photoabsorption_spectrum(dmfile, specfile)

if 0:
    # Reference RS-LR-TDDFT
    calc = GPAW(xc=xc, h=h, charge=c, width=0, nbands=4)
    atoms.set_calculator(calc)
    atoms.get_potential_energy()
    lr = lrtddft.LrTDDFT(calc, finegrid=0)
    lr.diagonalize()
    if world.rank == 0:
        lrtddft.photoabsorption_spectrum(lr, sy+'_rs_lr.spectrum', e_min=0.0, e_max=40)

    # Reference LCAO-LR-TDDFT
    calc = GPAW(mode='lcao', xc=xc, h=h, basis=b, charge=c, width=0)
    atoms.set_calculator(calc)
    atoms.get_potential_energy()
    lr = lrtddft.LrTDDFT(calc, finegrid=0)
    lr.diagonalize()
    if world.rank == 0:
        lrtddft.photoabsorption_spectrum(lr, sy+'_lcao_'+b+'_lr.spectrum', e_min=0.0, e_max=400)

    # XXX Actually check that the spectums match
