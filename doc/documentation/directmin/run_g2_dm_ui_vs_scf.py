import warnings
from ase.collections import g2
from doc.documentation.directmin import tools_and_data
from gpaw import LCAO, ConvergenceError
from ase.parallel import paropen

# Results (total energy, number of iterations) obtained
# in a previous calculation. Used to compare with the
# current results.
saved_data = \
    {0: tools_and_data.read_saved_data(tools_and_data.data_g2_scf),
     1: tools_and_data.read_saved_data(tools_and_data.data_g2_dm)}

calc_args = {'xc': 'PBE', 'h': 0.15,
             'convergence': {'density': 1.0e-6,
                             'eigenstates': 100},
             'maxiter': 333, 'basis': 'dzp',
             'mode': LCAO(), 'symmetry': 'off'}

eig_string = ['scf', 'dm']
with paropen('dm-g2-results.txt', 'w') as fdm, \
        paropen('scf-g2-results.txt', 'w') as fscf:
    fd = {0: fscf, 1: fdm}
    for name in saved_data[0].keys():
        atoms = g2[name]
        atoms.center(vacuum=7.0)
        for dm in [0, 1]:
            txt = name + eig_string[dm] + '.txt'
            tools_and_data.set_calc(atoms, calc_args, txt, dm)

            try:
                e, iters, t = tools_and_data.get_energy_and_iters(atoms, dm)

                # Compare with saved data from previous calculation
                e_diff_saved_calc = abs(saved_data[dm][name][1] - e)
                iters_diff_saved_calc = abs(saved_data[dm][name][0] - iters)
                if e_diff_saved_calc > 1.0e-2:
                    warnings.warn('Absolute difference in total energy '
                                  'for ' + eig_string[dm] + ' calculation of '
                                  + name + ' with respect to saved results '
                                  'is %f eV'
                                  % e_diff_saved_calc)
                if iters_diff_saved_calc > 3:
                    warnings.warn('Absolute difference in total number of '
                                  'iterations for ' + eig_string[dm] +
                                  ' calculation of ' + name + ' with respect '
                                  'to saved results is %d'
                                  % iters_diff_saved_calc)

                print(name + "\t{}".format(iters),
                      file=fd[dm], flush=True)

            except ConvergenceError:
                print(name + "\t{}".format(None),
                      file=fd[dm], flush=True)
