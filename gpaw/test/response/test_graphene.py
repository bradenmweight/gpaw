import pytest
import numpy as np
from ase import Atoms

from gpaw import GPAW, PW, FermiDirac, Mixer
from gpaw.utilities import compiled_with_sl
from gpaw.response.df import DielectricFunction
from gpaw.mpi import world

# This test assures that some things that
# should be equal, are.


@pytest.mark.response
@pytest.mark.slow
def test_response_graphene(in_tmp_dir):
    a = 2.5
    c = 3.22

    GR = Atoms(symbols='C2',
               positions=[(0.5 * a, 0.2 - np.sqrt(3) / 6 * a, 0.0),
                          (0.5 * a, 0.2 + np.sqrt(3) / 6 * a, 0.0)],
               cell=[(0.5 * a, -0.5 * 3**0.5 * a, 0),
                     (0.5 * a, 0.5 * 3**0.5 * a, 0),
                     (0.0, 0.0, c * 2.0)])
    GR.set_pbc((True, True, True))
    atoms = GR
    # These are not physically-motivated kpoint grid sizes.
    # They are the (nearly) the minimum number we need to test
    # the symmetry off-on and gamma off-on difference.
    GSsettings = [
        {'symmetry': 'off', 'kpts': {'size': [3, 2, 1], 'gamma': False}},
        {'symmetry': {}, 'kpts': {'size': [3, 2, 1], 'gamma': False}},
        {'symmetry': 'off', 'kpts': {'size': [3, 2, 1], 'gamma': True}},
        {'symmetry': {}, 'kpts': {'size': [3, 2, 1], 'gamma': True}}]

    DFsettings = [{'disable_point_group': True,
                   'disable_time_reversal': True},
                  {'disable_point_group': False,
                   'disable_time_reversal': True},
                  {'disable_point_group': True,
                   'disable_time_reversal': False},
                  {'disable_point_group': False,
                   'disable_time_reversal': False}]

    if world.size > 1 and compiled_with_sl():
        DFsettings.append({'disable_point_group': False,
                           'disable_time_reversal': False,
                           'nblocks': 2})

    for GSkwargs in GSsettings:
        calc = GPAW(mode=PW(200),
                    occupations=FermiDirac(0.2),
                    mixer=Mixer(0.4),
                    convergence={'eigenstates': 1e-4, 'density': 1e-3},
                    **GSkwargs)

        atoms.calc = calc
        atoms.get_potential_energy()
        calc.write('gr.gpw', 'all')

        dfs = []
        for kwargs in DFsettings:
            DF = DielectricFunction(calc='gr.gpw',
                                    frequencies={'type': 'nonlinear',
                                                 'domega0': 0.2},
                                    eta=0.2,
                                    ecut=15.0,
                                    **kwargs)
            df1, df2 = DF.get_dielectric_function()
            if world.rank == 0:
                dfs.append(df1)

        # Check the calculated dielectric functions against
        # each other.
        while len(dfs):
            df = dfs.pop()
            for DFkwargs, df2 in zip(DFsettings[-len(dfs):], dfs):
                try:
                    assert np.allclose(df, df2)
                    print('Ground state settings:', GSkwargs)
                    print('DFkwargs1:', DFsettings[-len(dfs) - 1])
                    print('DFkwargs2:', DFkwargs)
                    print(np.max(np.abs((df - df2) / df)))
                except AssertionError:
                    print('Some symmetry or block-par. related problems')
                    print('for calculation with following ground state')
                    print('settings')
                    print('Ground state settings:', GSkwargs)
                    print('The following DF settings do not return the ' +
                          'same results')
                    print('DFkwargs1:', DFsettings[-len(dfs) - 1])
                    print('DFkwargs2:', DFkwargs)
                    print(np.max(np.abs((df - df2) / df)))
                    raise AssertionError
