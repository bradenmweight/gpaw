import numpy as np
import pytest

from ase.lattice.hexagonal import Graphene
from ase.parallel import parprint as pp

from gpaw import GPAW
from gpaw.response.df import DielectricFunction
from gpaw.mpi import world


@pytest.mark.response
@pytest.mark.skip(reason='TODO')
def test_graphene_EELS():
    system = Graphene(symbol='C',
                      latticeconstant={'a': 2.45, 'c': 1.0},
                      size=(1, 1, 1))
    system.pbc = (1, 1, 0)
    system.center(axis=2, vacuum=4.0)

    nkpts = 5

    communicator = world.new_communicator(np.array([world.rank]))
    gpwname = 'dump.graphene.gpw'

    if world.rank == 0:
        calc = GPAW(mode='pw',
                    kpts=(nkpts, nkpts, 1),
                    communicator=communicator,
                    xc='oldLDA',
                    nbands=len(system) * 6,
                    txt='gpaw.graphene.txt')
        system.calc = calc
        system.get_potential_energy()
        calc.write(gpwname, mode='all')

    world.barrier()

    parallel = dict(domain=(1, 1, 1), band=1)
    if world.size == 8:
        # parallel['domain'] = (1, 1, 2)
        parallel['band'] = 2
    calc = GPAW(gpwname,
                txt=None,
                parallel=parallel,
                idiotproof=False)
    pp('after restart')

    q = np.array([1.0 / nkpts, 0., 0.])
    w = np.linspace(0, 31.9, 320)
    dw = w[1] - w[0]

    def getpeak(energies, loss):
        arg = loss.argmax()
        energy = energies[arg]
        peakloss = loss[arg]
        return energy, peakloss

    scriptlines = []

    loss_errs = []
    energy_errs = []

    def check(name, energy, peakloss, ref_energy, ref_loss):
        pp('check %s :: energy = %5.2f [%5.2f], peakloss = %.12f [%.12f]'
           % (name, energy, ref_energy, peakloss, ref_loss))
        energy_errs.append(abs(energy - ref_energy))
        loss_errs.append(abs(peakloss - ref_loss))

    template = """\
    check_df('%s', %s, %s, %s, %s,
             **%s)"""

    def check_df(name, ref_energy, ref_loss, ref_energy_lfe, ref_loss_lfe,
                 **kwargs_override):
        kwargs = dict(calc=calc, frequencies=w.copy(), eta=0.5, ecut=30,
                      txt='df.%s.txt' % name)
        kwargs.update(kwargs_override)
        df = DielectricFunction(**kwargs)
        fname = 'dfdump.%s.dat' % name
        df.get_eels_spectrum('RPA', q_c=q, filename=fname)
        world.barrier()
        d = np.loadtxt(fname, delimiter=',')

        loss = d[:, 1]
        loss_lfe = d[:, 2]
        energies = d[:, 0]

        # import pylab as pl
        # fig = pl.figure()
        # ax1 = fig.add_subplot(111)
        # ax1.plot(d[:, 0], d[:, 1]/np.max(d[:, 1]))
        # ax1.plot(d[:, 0], d[:, 2]/np.max(d[:, 2]))
        # ax1.axis(ymin=0, ymax=1)
        # fig.savefig('fig.%s.pdf' % name)

        energy, peakloss = getpeak(energies, loss)
        energy_lfe, peakloss_lfe = getpeak(energies, loss_lfe)

        check(name, energy, peakloss, ref_energy, ref_loss)
        check('%s-lfe' % name, energy_lfe, peakloss_lfe, ref_energy_lfe,
              ref_loss_lfe)

        line = template % (name, energy, peakloss, energy_lfe, peakloss_lfe,
                           repr(kwargs_override))
        scriptlines.append(line)

    # These lines can be generated by the loop over 'scriptlines' below,
    # in case new reference values are wanted

    # The implementation for vcut choices has still to be done

    check_df('3d', 20.20, 2.505295593820, 26.90, 1.748517033160)
    #         **{'rpad': array([1, 1, 1])}) #, 'vcut': '3D'})
    # check_df('2d', 20.10, 2.449662058530, 26.80, 1.538080502420,
    #          **{'rpad': array([1, 1, 1]), 'vcut': '2D'})

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
