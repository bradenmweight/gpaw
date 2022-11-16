from __future__ import annotations

import importlib
from types import SimpleNamespace
from typing import Any, Union

import numpy as np
from ase import Atoms
from ase.calculators.calculator import kpts2sizeandoffsets
from ase.units import Bohr

from gpaw.core import UniformGrid
from gpaw.core.domain import Domain
from gpaw.core.atom_arrays import (AtomArrays, AtomArraysLayout,
                                   AtomDistribution)
from gpaw.mixer import MixerWrapper, get_mixer_from_keywords
from gpaw.mpi import MPIComm, Parallelization, serial_comm, world
from gpaw.new import cached_property, prod
from gpaw.new.basis import create_basis
from gpaw.new.brillouin import BZPoints, MonkhorstPackKPoints
from gpaw.new.density import Density
from gpaw.new.input_parameters import InputParameters
from gpaw.new.scf import SCFLoop
from gpaw.new.smearing import OccupationNumberCalculator
from gpaw.new.symmetry import create_symmetries_object
from gpaw.new.xc import XCFunctional
from gpaw.setup import Setups
from gpaw.typing import DTypeLike, Array2D, ArrayLike1D, ArrayLike2D
from gpaw.utilities.gpts import get_number_of_grid_points


def builder(atoms: Atoms,
            params: dict[str, Any] | InputParameters) -> DFTComponentsBuilder:
    """Create DFT-components builder.

    * pw
    * lcao
    * fd
    * tb
    * atom
    """
    if isinstance(params, dict):
        params = InputParameters(params)

    mode = params.mode.copy()
    name = mode.pop('name')
    assert name in {'pw', 'lcao', 'fd', 'tb', 'atom'}
    mod = importlib.import_module(f'gpaw.new.{name}.builder')
    name = name.title() if name == 'atom' else name.upper()
    return getattr(mod, f'{name}DFTComponentsBuilder')(atoms, params, **mode)


class DFTComponentsBuilder:
    def __init__(self,
                 atoms: Atoms,
                 params: InputParameters):

        self.atoms = atoms.copy()
        self.mode = params.mode['name']
        self.params = params

        parallel = params.parallel
        world = parallel['world']

        self.check_cell(atoms.cell)

        self.initial_magmom_av, self.ncomponents = normalize_initial_magmoms(
            atoms, params.magmoms, params.spinpol or params.hund)

        self.soc = params.soc
        self.nspins = self.ncomponents % 3
        self.spin_degeneracy = self.ncomponents % 2 + 1

        self.xc = self.create_xc_functional()

        self.setups = Setups(atoms.numbers,
                             params.setups,
                             params.basis,
                             self.xc.setup_name,
                             world)

        if params.hund:
            c = params.charge / len(atoms)
            for a, setup in enumerate(self.setups):
                self.initial_magmom_av[a, 2] = setup.get_hunds_rule_moment(c)

        symmetries = create_symmetries_object(atoms,
                                              self.setups.id_a,
                                              self.initial_magmom_av,
                                              params.symmetry)
        assert not (self.ncomponents == 4 and len(symmetries) > 1)
        bz = create_kpts(params.kpts, atoms)
        self.ibz = symmetries.reduce(bz, strict=False)

        d = parallel.get('domain', None)
        k = parallel.get('kpt', None)
        b = parallel.get('band', None)
        self.communicators = create_communicators(world, len(self.ibz),
                                                  d, k, b)

        if self.mode == 'fd':
            pass  # filter = create_fourier_filter(grid)
            # setups = stups.filter(filter)

        self.nelectrons = self.setups.nvalence - params.charge

        self.nbands = calculate_number_of_bands(params.nbands,
                                                self.setups,
                                                params.charge,
                                                self.initial_magmom_av,
                                                self.mode == 'lcao')
        if self.ncomponents == 4:
            self.nbands *= 2

        self.dtype: DTypeLike
        if params.force_complex_dtype:
            self.dtype = complex
        elif self.ibz.bz.gamma_only:
            self.dtype = float
        else:
            self.dtype = complex

        self.grid, self.fine_grid = self.create_uniform_grids()

        self.fracpos_ac = self.atoms.get_scaled_positions()
        self.fracpos_ac %= 1
        self.fracpos_ac %= 1

    def create_uniform_grids(self):
        raise NotImplementedError

    def create_xc_functional(self):
        return XCFunctional(self.params.xc, self.ncomponents)

    def check_cell(self, cell):
        number_of_lattice_vectors = cell.rank
        if number_of_lattice_vectors < 3:
            raise ValueError(
                'GPAW requires 3 lattice vectors.  '
                f'Your system has {number_of_lattice_vectors}.')

    @cached_property
    def atomdist(self) -> AtomDistribution:
        return AtomDistribution(
            self.grid.ranks_from_fractional_positions(self.fracpos_ac),
            self.grid.comm)

    @cached_property
    def wf_desc(self) -> Domain:
        return self.create_wf_description()

    def __repr__(self):
        return f'{self.__class__.__name__}({self.atoms}, {self.params})'

    @cached_property
    def nct_R(self):
        out = self.grid.empty()
        nct_aX = self.get_pseudo_core_densities()
        nct_aX.to_uniform_grid(out=out,
                               scale=1.0 / (self.ncomponents % 3))
        return out

    def create_basis_set(self):
        return create_basis(self.ibz,
                            self.ncomponents % 3,
                            self.atoms.pbc,
                            self.grid,
                            self.setups,
                            self.dtype,
                            self.fracpos_ac,
                            self.communicators['w'],
                            self.communicators['k'])

    def density_from_superposition(self, basis_set):
        return Density.from_superposition(self.grid,
                                          self.nct_R,
                                          self.atomdist,
                                          self.setups,
                                          basis_set,
                                          self.initial_magmom_av,
                                          self.ncomponents,
                                          self.params.charge,
                                          self.params.hund)

    def create_occupation_number_calculator(self):
        return OccupationNumberCalculator(
            self.params.occupations,
            self.atoms.pbc,
            self.ibz,
            self.nbands,
            self.communicators,
            self.initial_magmom_av.sum(0),
            self.ncomponents,
            np.linalg.inv(self.atoms.cell.complete()).T)

    def create_scf_loop(self):
        hamiltonian = self.create_hamiltonian_operator()
        eigensolver = self.create_eigensolver(hamiltonian)

        mixer = MixerWrapper(
            get_mixer_from_keywords(self.atoms.pbc.any(),
                                    self.ncomponents, **self.params.mixer),
            self.ncomponents, self.grid._gd)

        occ_calc = self.create_occupation_number_calculator()
        return SCFLoop(hamiltonian, occ_calc,
                       eigensolver, mixer, self.communicators['w'],
                       {key: value
                        for key, value in self.params.convergence.items()
                        if key != 'bands'},
                       self.params.maxiter)

    def read_ibz_wave_functions(self, reader):
        raise NotImplementedError

    def create_potential_calculator(self):
        raise NotImplementedError

    def read_wavefunction_values(self, reader, ibzwfs):
        """ Read eigenvalues, occuptions and projections and fermi levels

        The values are read using reader and set as the appropriate properties
        of (the already instantiated) wavefunctions contained in ibzwfs
        """
        ha = reader.ha

        eig_skn = reader.wave_functions.eigenvalues
        occ_skn = reader.wave_functions.occupations
        P_sknI = reader.wave_functions.projections

        if self.params.force_complex_dtype:
            P_sknI = P_sknI.astype(complex)

        for wfs in ibzwfs:
            wfs._eig_n = eig_skn[wfs.spin, wfs.k] / ha
            wfs._occ_n = occ_skn[wfs.spin, wfs.k]
            layout = AtomArraysLayout([(setup.ni,) for setup in self.setups],
                                      dtype=self.dtype)
            if self.ncomponents < 4:
                wfs._P_ani = AtomArrays(layout,
                                        dims=(self.nbands,),
                                        data=P_sknI[wfs.spin, wfs.k])
            else:
                wfs._P_ani = AtomArrays(layout,
                                        dims=(self.nbands, 2),
                                        data=P_sknI[wfs.k])

        ibzwfs.fermi_levels = reader.wave_functions.fermi_levels / ha


def create_communicators(comm: MPIComm = None,
                         nibzkpts: int = 1,
                         domain: Union[int, tuple[int, int, int]] = None,
                         kpt: int = None,
                         band: int = None) -> dict[str, MPIComm]:
    parallelization = Parallelization(comm or world, nibzkpts)
    if domain is not None and not isinstance(domain, int):
        domain = prod(domain)
    parallelization.set(kpt=kpt,
                        domain=domain,
                        band=band)
    comms = parallelization.build_communicators()
    comms['w'] = comm
    return comms


def create_fourier_filter(grid):
    gamma = 1.6

    h = ((grid.icell**2).sum(1)**-0.5 / grid.size).max()

    def filter(rgd, rcut, f_r, l=0):
        gcut = np.pi / h - 2 / rcut / gamma
        ftmp = rgd.filter(f_r, rcut * gamma, gcut, l)
        f_r[:] = ftmp[:len(f_r)]

    return filter


def normalize_initial_magmoms(
        atoms: Atoms,
        magmoms: ArrayLike2D | ArrayLike1D | float | None = None,
        force_spinpol_calculation: bool = False) -> tuple[Array2D, int]:
    """Convert magnetic moments to (natoms, 3)-shaped array.

    Also return number of wave function components (1, 2 or 4).

    >>> h = Atoms('H', magmoms=[1])
    >>> normalize_initial_magmoms(h)
    (array([[0., 0., 1.]]), 2)
    >>> normalize_initial_magmoms(h, [[1, 0, 0]])
    (array([[1., 0., 0.]]), 4)
    """
    magmom_av = np.zeros((len(atoms), 3))
    ncomponents = 2

    if magmoms is None:
        magmom_av[:, 2] = atoms.get_initial_magnetic_moments()
    elif isinstance(magmoms, float):
        magmom_av[:, 2] = magmoms
    else:
        magmoms = np.asarray(magmoms)
        if magmoms.ndim == 1:
            magmom_av[:, 2] = magmoms
        else:
            magmom_av[:] = magmoms
            ncomponents = 4

    if (ncomponents == 2 and
        not force_spinpol_calculation and
        not magmom_av[:, 2].any()):
        ncomponents = 1

    return magmom_av, ncomponents


def create_kpts(kpts: dict[str, Any], atoms: Atoms) -> BZPoints:
    if 'points' in kpts:
        assert len(kpts) == 1, kpts
        return BZPoints(kpts['points'])
    size, offset = kpts2sizeandoffsets(**kpts, atoms=atoms)
    return MonkhorstPackKPoints(size, offset)


def calculate_number_of_bands(nbands: int | str | None,
                              setups: Setups,
                              charge: float,
                              initial_magmom_av: Array2D,
                              is_lcao: bool) -> int:
    nao = setups.nao
    nvalence = setups.nvalence - charge
    M = np.linalg.norm(initial_magmom_av.sum(0))

    orbital_free = any(setup.orbital_free for setup in setups)
    if orbital_free:
        return 1

    if isinstance(nbands, str):
        if nbands == 'nao':
            N = nao
        elif nbands[-1] == '%':
            cfgbands = (nvalence + M) / 2
            N = int(np.ceil(float(nbands[:-1]) / 100 * cfgbands))
        else:
            raise ValueError('Integer expected: Only use a string '
                             'if giving a percentage of occupied bands')
    elif nbands is None:
        # Number of bound partial waves:
        nbandsmax = sum(setup.get_default_nbands()
                        for setup in setups)
        N = int(np.ceil((1.2 * (nvalence + M) / 2))) + 4
        if N > nbandsmax:
            N = nbandsmax
        if is_lcao and N > nao:
            N = nao
    elif nbands <= 0:
        N = max(1, int(nvalence + M + 0.5) // 2 + (-nbands))
    else:
        N = nbands

    if N > nao and is_lcao:
        raise ValueError('Too many bands for LCAO calculation: '
                         f'{nbands}%d bands and only {nao} atomic orbitals!')

    if nvalence < 0:
        raise ValueError(
            f'Charge {charge} is not possible - not enough valence electrons')

    if nvalence > 2 * N:
        raise ValueError(
            f'Too few bands!  Electrons: {nvalence}, bands: {nbands}')

    return N


def create_uniform_grid(mode: str,
                        gpts,
                        cell,
                        pbc,
                        symmetry,
                        h: float = None,
                        interpolation: str = None,
                        ecut: float = None,
                        comm: MPIComm = serial_comm) -> UniformGrid:
    """Create grid in a backwards compatible way."""
    cell = cell / Bohr
    if h is not None:
        h /= Bohr

    realspace = (mode != 'pw' and interpolation != 'fft')
    if not realspace:
        pbc = (True, True, True)

    if gpts is not None:
        size = gpts
    else:
        modeobj = SimpleNamespace(name=mode, ecut=ecut)
        size = get_number_of_grid_points(cell, h, modeobj, realspace,
                                         symmetry.symmetry)
    return UniformGrid(cell=cell, pbc=pbc, size=size, comm=comm)
