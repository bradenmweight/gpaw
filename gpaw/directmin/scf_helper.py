import warnings
import numpy as np

from ase.units import Ha
from gpaw.directmin.tools import (sort_orbitals_according_to_energies,
                                  get_n_occ)


def do_if_converged(eigensolver_name, wfs, ham, dens, log):
    if eigensolver_name == 'etdm':
        if hasattr(wfs.eigensolver, 'e_sic'):
            e_sic = wfs.eigensolver.e_sic
        else:
            e_sic = 0.0
        energy = ham.get_energy(
            0.0, wfs, kin_en_using_band=False, e_sic=e_sic)
        wfs.calculate_occupation_numbers(dens.fixed)
        wfs.eigensolver.get_canonical_representation(
            ham, wfs, dens, sort_eigenvalues=True)
        wfs.eigensolver.update_ks_energy(ham, wfs, dens)
        energy_converged = ham.get_energy(
            0.0, wfs, kin_en_using_band=False, e_sic=e_sic)
        energy_diff_after_scf = abs(energy - energy_converged) * Ha
        if energy_diff_after_scf > 1.0e-6:
            warnings.warn('Jump in energy of %f eV detected at the end of '
                          'SCF after getting canonical orbitals, SCF '
                          'might have converged to the wrong solution '
                          'or achieved energy convergence to the correct '
                          'solution above 1.0e-6 eV'
                          % (energy_diff_after_scf))

        log('\nOccupied states converged after'
            ' {:d} e/g evaluations'.format(wfs.eigensolver.eg_count))

    elif eigensolver_name == 'etdm-fdpw':
        solver = wfs.eigensolver
        solver.choose_optimal_orbitals(wfs)
        niter1 = solver.eg_count
        niter2 = 0
        niter3 = 0

        iloop1 = solver.iloop is not None
        iloop2 = solver.outer_iloop is not None
        if iloop1:
            niter2 = solver.total_eg_count_iloop
        if iloop2:
            niter3 = solver.total_eg_count_outer_iloop

        if iloop1 and iloop2:
            log(
                '\nOccupied states converged after'
                ' {:d} KS and {:d} SIC e/g '
                'evaluations'.format(niter3,
                                     niter2 + niter3))
        elif not iloop1 and iloop2:
            log(
                '\nOccupied states converged after'
                ' {:d} e/g evaluations'.format(niter3))
        elif iloop1 and not iloop2:
            log(
                '\nOccupied states converged after'
                ' {:d} KS and {:d} SIC e/g '
                'evaluations'.format(niter1, niter2))
        else:
            log(
                '\nOccupied states converged after'
                ' {:d} e/g evaluations'.format(niter1))
        if solver.converge_unocc:
            log('Converge unoccupied states:')
            max_er = wfs.eigensolver.error
            max_er *= Ha ** 2 / wfs.nvalence
            solver.run_unocc(ham, wfs, dens, max_er, log)
        else:
            solver.initialized = False
            log('Unoccupied states are not converged.')

        rewrite_psi = True
        sic_calc = 'SIC' in solver.func_settings['name']
        if sic_calc:
            rewrite_psi = False

        solver.get_canonical_representation(ham, wfs, rewrite_psi)

        occ_name = getattr(wfs.occupations, 'name', None)
        if occ_name == 'mom':
            f_sn = wfs.occupations.update_occupations()
            for kpt in wfs.kpt_u:
                k = wfs.kd.nibzkpts * kpt.s + kpt.q
                n_occ, occupied = get_n_occ(kpt)
                if n_occ != 0.0 and np.min(f_sn[k][:n_occ]) == 0:
                    warnings.warn('MOM has detected variational collapse '
                                  'after getting canonical orbitals. Check '
                                  'that the orbitals are consistent with the '
                                  'initial guess.')

        solver.get_energy_and_tangent_gradients(
            ham, wfs, dens)

        if occ_name == 'mom' and not sic_calc:
            # Sort orbitals according to eigenvalues
            sort_orbitals_according_to_energies(ham, wfs, use_eps=True)
            not_update = not wfs.occupations.update_numbers
            fixed_occ = wfs.occupations.use_fixed_occupations
            if not_update or fixed_occ:
                wfs.occupations.numbers = solver.initial_occupation_numbers


def check_eigensolver_state(eigensolver_name, wfs, ham, dens, log):

    solver = wfs.eigensolver
    name = eigensolver_name
    if name == 'etdm' or name == 'etdm-fdpw':
        solver.eg_count = 0
        solver.globaliters = 0

        if hasattr(solver, 'iloop'):
            if solver.iloop is not None:
                solver.iloop.total_eg_count = 0
        if hasattr(solver, 'outer_iloop'):
            if solver.outer_iloop is not None:
                solver.outer_iloop.total_eg_count = 0

        solver.check_assertions(wfs, dens)
        if (hasattr(solver, 'dm_helper') and solver.dm_helper is None) \
                or not solver.initialized:
            solver.initialize_dm_helper(wfs, ham, dens, log)
