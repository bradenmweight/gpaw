import sys
import numpy as np

from ase.lattice.hexagonal import Graphene
#import ase.io as io
#from ase import Atoms, Atom
from ase.visualize import view
from ase.parallel import parprint as pp

from gpaw import GPAW, restart
from gpaw.response.df import DF
from gpaw.mpi import world
from gpaw.version import version
from ase.structure import bulk

system = Graphene(symbol='C',
                  latticeconstant={'a': 2.45,'c': 1.0},
                  size=(1,1,1))
system.pbc = (1, 1, 0)
system.center(axis=2, vacuum=4.0)

nkpts = 5


communicator = world.new_communicator(np.arange(1))
gpwname = 'dump.graphene.gpw'

if world.rank == 0:
    calc = GPAW(kpts=(nkpts, nkpts, 1),
                communicator=communicator,
                h=0.24,
                xc='oldLDA',
                nbands=len(system) * 6,
                txt='gpaw.graphene.txt')
    system.set_calculator(calc)
    system.get_potential_energy()
    calc.write(gpwname, mode='all')

world.barrier()

parallel = dict(domain=(1, 1, 1), band=1)
if world.size == 8:
    #parallel['domain'] = (1, 1, 2)
    parallel['band'] = 2
calc = GPAW(gpwname,
            txt=None,
            parallel=parallel,
            idiotproof=False)
pp('after restart')

q = np.array([1.0 / nkpts, 0., 0.])
w = np.linspace(0, 31.9, 320)
dw = w[1] - w[0]

def check(name, energies, loss, ref_energy, ref_peakloss):
    arg = loss.argmax()
    energy = energies[arg]
    peakloss = loss[arg]
    pp('check %s :: energy = %5.2f [%5.2f], peakloss = %.12f [%.12f]' 
       % (name, energy, ref_energy, peakloss, ref_peakloss))

    data = dict(name=name, peakloss=peakloss, energy=energy,
                ref_peakloss=ref_peakloss, ref_energy=ref_energy)
    return data


def getpeak(energies, loss):
    arg = loss.argmax()
    energy = energies[arg]
    peakloss = loss[arg]
    return energy, peakloss

scriptlines = []

def check(name, energy, peakloss, ref_energy, ref_loss):
    pp('check %s :: energy = %5.2f [%5.2f], peakloss = %.12f [%.12f]' 
       % (name, energy, ref_energy, peakloss, ref_loss))
    energy_errs.append(energy - ref_energy)
    loss_errs.append(peakloss - ref_loss)
    

array = np.array
template = """\
check_df('%s', %s, %s, %s, %s,
         **%s)"""

loss_errs = []
energy_errs = []

def check_df(name, ref_energy, ref_loss, ref_energy_lfe, ref_loss_lfe,
             **kwargs_override):
    kwargs = dict(calc=calc, q=q, w=w.copy(), eta=0.5, ecut=(30, 30, 30),
                  txt='df.%s.txt' % name)
    kwargs.update(kwargs_override)
    df = DF(**kwargs)
    fname = 'dfdump.%s.dat' % name
    df.get_EELS_spectrum(filename=fname)
    world.barrier()
    d = np.loadtxt(fname)

    loss = d[:, 1]
    loss_lfe = d[:, 2]
    energies = d[:, 0]

    #import pylab as pl
    #fig = pl.figure()
    #ax1 = fig.add_subplot(111)
    #ax1.plot(d[:, 0], d[:, 1]/np.max(d[:, 1]))
    #ax1.plot(d[:, 0], d[:, 2]/np.max(d[:, 2]))
    #ax1.axis(ymin=0, ymax=1)
    #fig.savefig('fig.%s.pdf' % name)

    energy, peakloss = getpeak(energies, loss)
    energy_lfe , peakloss_lfe = getpeak(energies, loss_lfe)

    check(name, energy, peakloss, ref_energy, ref_loss)
    check('%s-lfe' % name, energy_lfe, peakloss_lfe, ref_energy_lfe,
          ref_loss_lfe)

    line = template % (name, energy, peakloss, energy_lfe, peakloss_lfe,
                       repr(kwargs_override))
    scriptlines.append(line)


# These lines can be generated by the loop over 'scriptlines' below,
# in case new reference values are wanted
check_df('3d', 20.8, 2.11246190374, 27.2, 1.9156598619,
         **{'rpad': array([1, 1, 1])})
check_df('2d', 20.1, 1.98584827779, 27.1, 1.67243546061,
         **{'rpad': array([1, 1, 1]), 'vcut': '2D'})
check_df('3d-rpad2', 17.8, 1.18466868855, 7.4, 0.638418137197,
         **{'rpad': array([1, 1, 2])})
check_df('2d-rpad2', 17.8, 1.18129266083, 7.4, 0.634712407924,
         **{'rpad': array([1, 1, 2]), 'vcut': '2D'})

pp()
pp('Insert lines into script to set new reference values:')
for line in scriptlines:
    pp(line)
pp()

for err in energy_errs:
    # with the current grid this error just means 
    assert err < dw / 4.0, err
for err in loss_errs:
    assert err < 1e-6, err
