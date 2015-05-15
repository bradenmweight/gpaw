"""
tests solvation parameters as in

O. Andreussi, I. Dabo, and N. Marzari,
The Journal of Chemical Physics, vol. 136, no. 6, p. 064102, 2012
"""

from gpaw import GPAW
from gpaw.cluster import Cluster
from gpaw.test import equal
from ase.structure import molecule
from ase.units import mol, kcal, Pascal, m, Bohr
from gpaw.solvation import (
    SolvationGPAW,
    ADM12SmoothStepCavity,
    LinearDielectric,
    GradientSurface,
    SurfaceInteraction,
    KB51Volume,
    VolumeInteraction,
    ElDensity
)
from gpaw.solvation.poisson import ADM12PoissonSolver
import warnings

SKIP_VAC_CALC = True

h = 0.24
vac = 4.0

epsinf = 78.36
rhomin = 0.0001 / Bohr ** 3
rhomax = 0.0050 / Bohr ** 3
st = 50. * 1e-3 * Pascal * m
p = -0.35 * 1e9 * Pascal

atoms = Cluster(molecule('H2O'))
atoms.minimal_box(vac, h)

if not SKIP_VAC_CALC:
    atoms.calc = GPAW(xc='PBE', h=h)
    Evac = atoms.get_potential_energy()
    print Evac
else:
    Evac = -14.8628983897  # h = 0.24, vac = 4.0

with warnings.catch_warnings():
    # ignore production code warning for ADM12PoissonSolver
    warnings.simplefilter("ignore")
    psolver = ADM12PoissonSolver()

atoms.calc = SolvationGPAW(
    xc='PBE', h=h, poissonsolver=psolver,
    cavity=ADM12SmoothStepCavity(
        rhomin=rhomin, rhomax=rhomax, epsinf=epsinf,
        density=ElDensity(),
        surface_calculator=GradientSurface(),
        volume_calculator=KB51Volume()
        ),
    dielectric=LinearDielectric(epsinf=epsinf),
    interactions=[
        SurfaceInteraction(surface_tension=st),
        VolumeInteraction(pressure=p)
        ]
    )
Ewater = atoms.get_potential_energy()
assert atoms.calc.get_number_of_iterations() < 40
atoms.get_forces()
DGSol = (Ewater - Evac) / (kcal / mol)
print 'Delta Gsol: %s kcal / mol' % (DGSol, )

equal(DGSol, -6.3, 2.)