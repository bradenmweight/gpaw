import numpy as np
from ase import Atoms
from gpaw import GPAW
from gpaw.tddft import TDDFT
from gpaw.inducedfield.inducedfield_tddft import TDDFTInducedField
from gpaw.poisson import PoissonSolver
from gpaw.test import equal

do_print_values = False  # Use this for printing the reference values
poisson_eps = 1e-20
density_eps = 1e-8

# PoissonSolver
poissonsolver = PoissonSolver(eps=poisson_eps)

# Na2 cluster
atoms = Atoms(symbols='Na2',
              positions=[(0, 0, 0), (3.0, 0, 0)],
              pbc=False)
atoms.center(vacuum=3.0)

# Standard ground state calculation
calc = GPAW(nbands=2, h=0.6, setups={'Na': '1'}, poissonsolver=poissonsolver,
        convergence={'density': density_eps})
atoms.set_calculator(calc)
energy = atoms.get_potential_energy()
calc.write('na2_gs.gpw', mode='all')

# Standard time-propagation initialization
time_step = 10.0
iterations = 60
kick_strength = [1.0e-3, 1.0e-3, 0.0]
td_calc = TDDFT('na2_gs.gpw', poissonsolver=poissonsolver)
td_calc.absorption_kick(kick_strength=kick_strength)

# Create and attach InducedField object
frequencies = [1.0, 2.08]     # Frequencies of interest in eV
folding = 'Gauss'             # Folding function
width = 0.1                   # Line width for folding in eV
ind = TDDFTInducedField(paw=td_calc,
                        frequencies=frequencies,
                        folding=folding,
                        width=width,
                        restart_file='na2_td.ind')

# Propagate as usual
td_calc.propagate(time_step, iterations // 2, 'na2_td_dm.dat', 'na2_td.gpw')

# Save TDDFT and InducedField objects
td_calc.write('na2_td.gpw', mode='all')
ind.write('na2_td.ind')

ind.paw = None

# Restart and continue
td_calc = TDDFT('na2_td.gpw', poissonsolver=poissonsolver)

# Load and attach InducedField object
ind = TDDFTInducedField(filename='na2_td.ind',
                        paw=td_calc,
                        restart_file='na2_td.ind')
    
# Continue propagation as usual
td_calc.propagate(time_step, iterations // 2, 'na2_td_dm.dat', 'na2_td.gpw')

# Calculate induced electric field
ind.calculate_induced_field(gridrefinement=2, from_density='comp',
                            poisson_eps=poisson_eps)

# Estimate tolerance (worst case error accumulation)
tol = iterations * ind.fieldgd.integrate(ind.fieldgd.zeros() + 1.0) * \
        max(density_eps, np.sqrt(poisson_eps))
#tol = 0.001167179811
if do_print_values:
    print('tol = %.12f' % tol)

# Test
val1 = ind.fieldgd.integrate(ind.Ffe_wg[0])
val2 = ind.fieldgd.integrate(np.abs(ind.Fef_wvg[0][0]))
val3 = ind.fieldgd.integrate(np.abs(ind.Fef_wvg[0][1]))
val4 = ind.fieldgd.integrate(np.abs(ind.Fef_wvg[0][2]))
val5 = ind.fieldgd.integrate(ind.Ffe_wg[1])
val6 = ind.fieldgd.integrate(np.abs(ind.Fef_wvg[1][0]))
val7 = ind.fieldgd.integrate(np.abs(ind.Fef_wvg[1][1]))
val8 = ind.fieldgd.integrate(np.abs(ind.Fef_wvg[1][2]))

if do_print_values:
    i = 1
    def equal(x, y, tol):
        global i
        print("equal(val%d, %20.12f, tol)" % (i, x))
        i += 1

equal(val1,    3314.463566436767, tol)
equal(val2,       2.254452178599, tol)
equal(val3,       2.208970141543, tol)
equal(val4,       1.697564042167, tol)
equal(val5,    3311.757129629823, tol)
equal(val6,       2.102499949779, tol)
equal(val7,       2.103871112604, tol)
equal(val8,       1.602251056016, tol)

ind.paw = None
