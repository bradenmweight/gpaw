from ase.build import molecule
from gpaw import GPAW
import gpaw.dscf as dscf

# Ground state calculation
calc_params = dict(
    mode='fd',
    nbands=8,
    h=0.2,
    xc='PBE',
    spinpol=True,
    convergence={'energy': 100,
                 'density': 100,
                 'bands': -1})
calc_mol = GPAW(**calc_params)

CO = molecule('CO')
CO.center(vacuum=3)
CO.calc = calc_mol
E_gs = CO.get_potential_energy()

# Get the pseudowavefunctions and projector overlaps of the
# state which is to be occupied. n=5,6 is the 2pix and 2piy orbitals
n = 5
molecule = [0, 1]
wf_u = [kpt.psit_nG[n] for kpt in calc_mol.wfs.kpt_u]
p_uai = [dict([(molecule[a], P_ni[n]) for a, P_ni in kpt.P_ani.items()])
         for kpt in calc_mol.wfs.kpt_u]

# Excited state calculations
calc_1 = GPAW(**calc_params)
CO.calc = calc_1
weights = {0: [0., 0., 0., 1.], 1: [0., 0., 0., -1.]}
lumo = dscf.MolecularOrbital(calc_1, weights=weights)
dscf.dscf_calculation(calc_1, [[1.0, lumo, 1]], CO)
E_es1 = CO.get_potential_energy()
calc_1.write('dscf_CO_es1.gpw', mode='all')

calc_2 = GPAW(**calc_params)
CO.calc = calc_2
lumo = dscf.AEOrbital(calc_2, wf_u, p_uai)
dscf.dscf_calculation(calc_2, [[1.0, lumo, 1]], CO)
E_es2 = CO.get_potential_energy()
calc_2.write('dscf_CO_es2.gpw', mode='all')

assert abs(E_es1 - (E_gs + 5.8)) < 0.1
assert abs(E_es1 - E_es2) < 0.001
