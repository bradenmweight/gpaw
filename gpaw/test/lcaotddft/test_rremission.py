import numpy as np

from ase.build import molecule
from gpaw import GPAW
from gpaw.lcaotddft import LCAOTDDFT
from gpaw.lcaotddft.dipolemomentwriter import DipoleMomentWriter
from gpaw.lcaotddft.qed import RRemission
from gpaw.mpi import world
from .test_molecule import calculate_error


def check_mm(ref_fpath, data_fpath, atol):
    world.barrier()
    ref = np.loadtxt(ref_fpath)
    data = np.loadtxt(data_fpath)
    err = calculate_error(data, ref)
    assert err < atol


def test_lcaotddft_simple(in_tmp_dir):
    atoms = molecule('Na2')
    atoms.center(vacuum=4.0)
    calc = GPAW(mode='lcao', h=0.4, basis='dzp',
                setups={'Na': '1'},
                convergence={'density': 1e-12})
    atoms.calc = calc
    atoms.get_potential_energy()
    calc.write('gs.gpw', mode='all')
    td_calc = LCAOTDDFT('gs.gpw', rremission=RRemission(0.5, [0, 0, 1]))
    DipoleMomentWriter(td_calc, 'dm.dat')
    td_calc.absorption_kick([0.0, 0.0, 1e-5])
    td_calc.propagate(40, 20)
    world.barrier()

    with open('dm_ref.dat', 'w') as f:
        f.write('''
# DipoleMomentWriter[version=1](center=False, density='comp')
#            time            norm                    dmx                    dmy                    dmz
# Start; Time = 0.00000000
          0.00000000      -1.06086854e-15     8.084419232390e-15     1.010552404049e-14    -2.070038560300e-14
# Kick = [    0.000000000000e+00,     0.000000000000e+00,     1.000000000000e-05]; Time = 0.00000000
          0.00000000      -5.31739841e-16     2.309834066397e-15     1.154917033199e-15    -2.472546058136e-14
          1.65365493      -7.01836907e-16     6.929502199192e-15     7.795689974091e-15     3.419488133320e-05
          3.30730987      -1.23152514e-15     1.010552404049e-14     7.795689974091e-15     6.321572810574e-05
          4.96096480      -4.03794022e-16     5.774585165993e-15     3.464751099596e-15     8.589211243093e-05
          6.61461974      -3.96893154e-16     4.619668132795e-15     3.176021841296e-15     1.030156529177e-04
          8.26827467      -7.47531843e-16     5.774585165993e-15     4.619668132795e-15     1.153269761059e-04
          9.92192960      -1.53870701e-15     1.299281662348e-14     1.270408736518e-14     1.235126215245e-04
         11.57558454      -7.51075532e-16     1.010552404049e-14     1.039425329879e-14     1.282024146875e-04
         13.22923947      -6.98199963e-16     2.887292582997e-15    -0.000000000000e+00     1.299679651227e-04
         14.88289440      -7.46039763e-16     2.598563324697e-15     6.640772940892e-15     1.293224201649e-04
         16.53654934      -1.74107030e-16     1.154917033199e-15     1.443646291498e-15     1.267215154536e-04
         18.19020427      -2.12621333e-16     5.774585165993e-16     2.309834066397e-15     1.225658544645e-04
         19.84385921      -9.03827173e-16     6.929502199192e-15     6.929502199192e-15     1.172042561434e-04
         21.49751414      -7.79145278e-16     2.309834066397e-15     3.753480357896e-15     1.109379443511e-04
         23.15116907      -2.12434823e-16     4.619668132795e-15     5.774585165993e-15     1.040253093614e-04
         24.80482401      -9.08489922e-16     3.464751099596e-15     4.330938874495e-15     9.668693730311e-05
         26.45847894      -9.38145002e-16     6.352043682592e-15     8.084419232390e-15     8.911057664562e-05
         28.11213387      -3.02425869e-16     3.464751099596e-15     6.929502199192e-15     8.145570238167e-05
         29.76578881       3.49706139e-17    -0.000000000000e+00    -8.661877748990e-16     7.385737279135e-05
         31.41944374       1.28878369e-16     2.021104808098e-15    -2.887292582997e-15     6.642917009752e-05
         33.07309868      -9.49335599e-16     7.218231457491e-15     4.908397391094e-15     5.926517291299e-05
'''.strip())  # noqa: E501

    check_mm('dm.dat', 'dm_ref.dat', atol=5e-14)
