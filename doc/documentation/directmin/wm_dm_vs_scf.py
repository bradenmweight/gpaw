# Generate the data visualized in web-page: water.png
import numpy as np
from doc.documentation.directmin import tools_and_data
from ase import Atoms
from gpaw import LCAO
from ase.parallel import paropen
from gpaw.mpi import world
from gpaw.atom.basis import BasisMaker
from gpaw import setup_paths
setup_paths.insert(0, '.')

for symbol in ['H', 'O']:
    bm = BasisMaker(symbol, xc='PBE')
    basis = bm.generate(zetacount=3, polarizationcount=2)
    basis.write_xml()

positions = tools_and_data.positions

L = 9.8553729
r = [[1, 1, 1],
     [2, 1, 1],
     [2, 2, 1]]
calc_args = {'xc': 'PBE', 'h': 0.2,
             'convergence': {'density': 1.0e-6,
                             'eigenstates': 100},
             'maxiter': 333, 'basis': 'tzdp',
             'mode': LCAO(), 'symmetry': 'off',
             'parallel': {'domain': world.size}}
# Results (total energy, number of iterations) obtained
# in a previous calculation. Used to compare with the
# current results.
saved_results = tools_and_data.wm_saved_results

t = np.zeros(2)
iters = np.zeros(2)
with paropen('water-results.txt', 'w') as fd:
    for i, x in enumerate(r):
        atoms = Atoms('32(OH2)', positions=positions)
        atoms.set_cell((L, L, L))
        atoms.set_pbc(1)
        atoms = atoms.repeat(x)
        for dm in [0, 1]:
            if dm:
                txt = str(len(atoms) // 3) + '_H2Omlcls_dm.txt'
            else:
                txt = str(len(atoms) // 3) + '_H2Omlcls_scf.txt',
            tools_and_data.set_calc(atoms, calc_args, txt, dm)

            e, iters[dm], t[dm] = \
                tools_and_data.get_energy_and_iters(atoms, dm)
            # Compare with saved data from previous calculation
            assert abs(saved_results[dm][i, 0] - e) < 1.0e-2
            assert abs(saved_results[dm][i, 1] - iters[dm]) < 3

        # Ratio of elapsed times per iteration
        # 2 is added to account for the diagonalization
        # performed at the beginning and at the end of etdm
        ratio_per_iter = (t[0] / iters[0]) / (t[1] / (iters[1] + 2))

        print("{}\t{}\t{}".format(
              len(atoms), t[0] / t[1], ratio_per_iter),
              flush=True, file=fd)
