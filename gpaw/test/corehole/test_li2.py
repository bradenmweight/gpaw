"""Test all electron density for right interpretation of coreholes"""
from ase.build import molecule
from ase.units import Bohr
from gpaw import GPAW, PoissonSolver
from gpaw.test import equal, gen
from gpaw.mixer import Mixer


def test_aed_with_corehole_li():
    """Compare number of electrons for different channels with corehole"""
    li_setup = gen('Li', name='fch1s', corehole=(1, 0, 1), xcname='PBE')
    grf = 1
    atoms = molecule('Li2')
    atoms.center(vacuum=2.5)

    calc = GPAW(xc='PBE',
                mixer=Mixer(),
                spinpol=False,
                setups={0: li_setup}, charge=-1,
                poissonsolver=PoissonSolver('fd'))
    atoms.calc = calc
    atoms.get_potential_energy()

    n_sg = calc.get_all_electron_density(gridrefinement=grf)

    ne_sz = calc.density.gd.integrate(
        n_sg, global_integral=False) * (Bohr / grf)**3
    equal(ne_sz, 6.0, 1e-8)

    atoms.set_initial_magnetic_moments([0.66, .34])
    calc.set(spinpol=True)
    atoms.get_potential_energy()

    for sz in range(2):
        n_sg = calc.get_all_electron_density(spin=sz, gridrefinement=grf)
        ne_sz = calc.density.gd.integrate(
            n_sg, global_integral=False) * (Bohr / grf)**3
        equal(ne_sz, 3.0, 1e-5)

    if not False:  # Did I break non-corehole in sp case?
        atoms.set_initial_magnetic_moments([-0.5, 0.5])

        calc = GPAW(xc='PBE',
                    mixer=Mixer(),
                    spinpol=True,
                    poissonsolver=PoissonSolver('fd'))
        atoms.calc = calc
        atoms.get_potential_energy()

        for sz in range(2):
            n_sg = calc.get_all_electron_density(spin=sz, gridrefinement=grf)
            ne_sz = calc.density.gd.integrate(
                n_sg, global_integral=False) * (Bohr / grf)**3
            equal(ne_sz, 3.0, 1e-8)


if __name__ == '__main__':
    test_aed_with_corehole_li()
