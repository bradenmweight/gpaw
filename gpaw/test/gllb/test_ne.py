from ase import Atoms, Atom
from gpaw import GPAW
from gpaw.atom.generator import Generator
from gpaw.atom.configurations import parameters
from gpaw import setup_paths
from gpaw.test import equal
from gpaw.mpi import world


def test_gllb_ne(in_tmp_dir):
    atom = 'Ne'
    setup_paths.insert(0, '.')

    for xcname in ['GLLBSC','GLLB']:
        if world.rank == 0:
            g = Generator(atom, xcname =xcname, scalarrel=False,nofiles=True)
            g.run(**parameters[atom])
            eps = g.e_j[-1]
        else:
            eps = 0.0
        eps = world.sum(eps)
        world.barrier()

        a = 5
        Ne = Atoms([Atom(atom, (0, 0, 0))],
                   cell=(a, a, a), pbc=False)
        Ne.center()
        calc = GPAW(nbands=7, h=0.25, xc=xcname)
        Ne.set_calculator(calc)
        e = Ne.get_potential_energy()
        # Calculate the discontinuity
        response = calc.hamiltonian.xc.xcs['RESPONSE']
        response.calculate_delta_xc()
        response.calculate_delta_xc_perturbation()

        eps3d = calc.wfs.kpt_u[0].eps_n[3]
        #if world.rank == 0:
        equal(eps, eps3d, 1e-3)
        # Correct for small cell +0.14eV (since the test needs to be fast
        # in test suite)
        equal(e+0.147106041, 0, 5e-2)