from __future__ import annotations
import numpy as np
from gpaw.hybrids import HybridXC
from ase.units import Bohr
from gpaw.mpi import MPIComm, Parallelization, world


def create_communicators(comm: MPIComm = None,
                         nibzkpts: int = 1,
                         domain: int | tuple[int, int, int] = None,
                         kpt: int = None,
                         band: int = None) -> dict[str, MPIComm]:
    parallelization = Parallelization(comm or world, nibzkpts)
    if domain is not None:
        domain = np.prod(domain)
    parallelization.set(kpt=kpt,
                        domain=domain,
                        band=band)
    comms = parallelization.build_communicators()
    return comms


def create_fourier_filter(grid):
    gamma = 1.6

    h = ((grid.icell**2).sum(1)**-0.5 / grid.size).max()

    def filter(rgd, rcut, f_r, l=0):
        gcut = np.pi / h - 2 / rcut / gamma
        ftmp = rgd.filter(f_r, rcut * gamma, gcut, l)
        f_r[:] = ftmp[:len(f_r)]

    return filter


class CalculationConfiguration:
    def __init__(self,
                 positions,
                 setups,
                 communicators,
                 grid,
                 xc,
                 ibz,
                 magmoms=None,
                 charge=0.0):
        self.positions = positions
        self.setups = setups
        self.magmoms = magmoms
        self.communicators = communicators
        self.ibz = ibz
        self.grid = grid
        self.xc = xc
        self.charge = charge

        self.grid2 = grid.new(size=grid.size * 2)
        # decomposition=[2 * d for d in grid.decomposition]

        self.band_comm = communicators['b']
        self.nelectrons = setups.nvalence - charge

    # Gather convergence criteria for SCF loop.
    custom = criteria.pop('custom', [])
    for name, criterion in criteria.items():
        if hasattr(criterion, 'todict'):
            # 'Copy' so no two calculators share an instance.
            criteria[name] = dict2criterion(criterion.todict())
        else:
            criteria[name] = dict2criterion({name: criterion})

    if not isinstance(custom, (list, tuple)):
        custom = [custom]
    for criterion in custom:
        if isinstance(criterion, dict):  # from .gpw file
            msg = ('Custom convergence criterion "{:s}" encountered, '
                   'which GPAW does not know how to load. This '
                   'criterion is NOT enabled; you may want to manually'
                   ' set it.'.format(criterion['name']))
            warnings.warn(msg)
            continue

        criteria[criterion.name] = criterion
        msg = ('Custom convergence criterion {:s} encountered. '
               'Please be sure that each calculator is fed a '
               'unique instance of this criterion. '
               'Note that if you save the calculator instance to '
               'a .gpw file you may not be able to re-open it. '
               .format(criterion.name))
        warnings.warn(msg)

    for criterion in criteria.values():
        criterion.reset()

    return criteria

    from gpaw.new.xc import XCFunctional
    from gpaw.xc import XC

    return XCFunctional(XC(value))


    from gpaw.new.modes import FDMode
    return FDMode()

    @staticmethod
    def from_parameters(atoms, params):
        parallel = params.parallel
        world = parallel['world']
        mode = params.mode
        xc = params.xc

        setups = params.setups(atoms.numbers,
                               params.basis,
                               xc.setup_name,
                               world)

        magmoms = params.magmoms(atoms)

        symmetry = params.symmetry(atoms, setups, magmoms)

        bz = params.kpts(atoms)
        ibz = symmetry.reduce(bz)

        d = parallel.pop('domain', None)
        k = parallel.pop('kpt', None)
        b = parallel.pop('band', None)

        if isinstance(xc, HybridXC):
            d = world.size

        communicators = create_communicators(world, len(ibz), d, k, b)
        communicators['w'] = world

        grid = mode.create_uniform_grid(params.h,
                                        params.gpts,
                                        atoms.cell / Bohr,
                                        atoms.pbc,
                                        symmetry,
                                        comm=communicators['d'])

        if mode.name == 'fd':
            pass  # filter = create_fourier_filter(grid)
            # setups = setups.filter(filter)

        return CalculationConfiguration(
            atoms.get_scaled_positions(),
            setups, communicators, grid, xc, ibz, magmoms,
            params.charge)
