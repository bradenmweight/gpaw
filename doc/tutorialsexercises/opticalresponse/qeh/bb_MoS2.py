from pathlib import Path
from gpaw.mpi import world
from gpaw.response.df import DielectricFunction
from gpaw.response.qeh import BuildingBlock

df = DielectricFunction(calc='MoS2_gs_fulldiag.gpw',
                        eta=0.001,
                        frequencies={'type': 'nonlinear',
                                     'domega0': 0.05,
                                     'omega2': 10.0},
                        nblocks=8,
                        ecut=150,
                        truncation='2D')

buildingblock = BuildingBlock('MoS2', df, qmax=3.0)

buildingblock.calculate_building_block()

if world.rank == 0:
    Path('MoS2_gs_fulldiag.gpw').unlink()
