from gpaw.cluster import Cluster
from ase.structure import molecule
from ase.data.vdw import vdw_radii
from gpaw.solvation import (
    SolvationGPAW,
    EffectivePotentialCavity,
    Power12Potential,
    LinearDielectric
)
from gpaw.solvation.poisson import ADM12PoissonSolver

h = 0.3
vac = 3.0
u0 = .180
epsinf = 80.
T = 298.15
vdw_radii = vdw_radii[:]
vdw_radii[1] = 1.09

atoms = Cluster(molecule('H2O'))
atoms.minimal_box(vac, h)
atomic_radii = [vdw_radii[n] for n in atoms.numbers]
atoms.pbc = True
atoms.calc = SolvationGPAW(
    xc='LDA', h=h,
    cavity=EffectivePotentialCavity(
        effective_potential=Power12Potential(atomic_radii=atomic_radii, u0=u0),
        temperature=T
        ),
    dielectric=LinearDielectric(epsinf=epsinf),
    poissonsolver=ADM12PoissonSolver(eps=1e-7)
    )
atoms.get_potential_energy()
atoms.get_forces()
