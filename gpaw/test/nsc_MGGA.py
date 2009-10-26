from ase import *
from gpaw import GPAW
from gpaw.test import equal, gen

# ??? g = Generator('H', 'TPSS', scalarrel=True, nofiles=True)

atoms = Atoms('H', magmoms=[1], pbc=True)
atoms.center(vacuum=3)
calc = GPAW(nbands=1, xc='PBE', txt='Hnsc.txt')
atoms.set_calculator(calc)
e1 = atoms.get_potential_energy()
niter1 = calc.get_number_of_iterations()
e1ref = calc.get_reference_energy()
de12t = calc.get_xc_difference('TPSS')
de12m = calc.get_xc_difference('M06L')


print '================'
print 'e1 = ', e1
print 'de12t = ', de12t
print 'de12m = ', de12m
print 'tpss = ', e1 + de12t
print 'm06l = ', e1 + de12m
print '================'

equal(e1 + de12t, -1.13140473115, 0.005)
equal(e1 + de12m, -1.19069632478, 0.005)

# ??? g = Generator('He', 'TPSS', scalarrel=True, nofiles=True)

atomsHe = Atoms('He', pbc=True)
atomsHe.center(vacuum=3)
calc = GPAW(nbands=1, xc='PBE', txt='Hensc.txt')
atomsHe.set_calculator(calc)
e1He = atomsHe.get_potential_energy()
niter_1He = calc.get_number_of_iterations()
e1refHe = calc.get_reference_energy()
de12tHe = calc.get_xc_difference('TPSS')
de12mHe = calc.get_xc_difference('M06L')

print '================'
print 'e1He = ', e1He
print 'de12tHe = ', de12tHe
print 'de12mHe = ', de12mHe
print 'tpss = ', e1He + de12tHe
print 'm06l = ', e1He + de12mHe
print '================'

equal(e1He+de12tHe, -0.448532905095, 0.005)
equal(e1He+de12mHe, -0.51400253951, 0.005)


energy_tolerance = 0.000001
niter_tolerance = 0
equal(e1, -1.12624857673, energy_tolerance) # svnversion 5252
equal(niter1, 25, niter_tolerance) # svnversion 5252
equal(e1He, 0.0104600878985, energy_tolerance) # svnversion 5252
equal(niter_1He, 10, niter_tolerance) # svnversion 5252
