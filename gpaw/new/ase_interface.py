from __future__ import annotations

import warnings
from functools import cached_property
from pathlib import Path
from types import SimpleNamespace
from typing import IO, Any, Union

import numpy as np
from ase import Atoms
from ase.units import Bohr, Ha
from gpaw import __version__
from gpaw.core import UGArray
from gpaw.dos import DOSCalculator
from gpaw.mpi import world, synchronize_atoms, broadcast as bcast
from gpaw.new import Timer, trace
from gpaw.new.builder import builder as create_builder
from gpaw.new.calculation import (DFTCalculation, DFTState,
                                  CalculationModeError,
                                  ReuseWaveFunctionsError, units)
from gpaw.new.gpw import read_gpw, write_gpw
from gpaw.new.input_parameters import (DeprecatedParameterWarning,
                                       InputParameters)
from gpaw.new.logger import Logger
from gpaw.new.pw.fulldiag import diagonalize
from gpaw.new.xc import create_functional
from gpaw.typing import Array1D, Array2D, Array3D
from gpaw.utilities import pack_density
from gpaw.utilities.memory import maxrss


def GPAW(filename: Union[str, Path, IO[str]] = None,
         txt: str | Path | IO[str] | None = '?',
         communicator=None,
         **kwargs) -> ASECalculator:
    """Create ASE-compatible GPAW calculator."""
    if txt == '?':
        txt = '-' if filename is None else None

    parallel = kwargs.get('parallel', {})
    comm = parallel.pop('world', None)
    if comm is None:
        comm = communicator or world
    else:
        warnings.warn(('Please use communicator=... '
                       'instead of parallel={''world'': ...}'),
                      DeprecatedParameterWarning)
    log = Logger(txt, comm)

    if filename is not None:
        if not {'parallel'}.issuperset(kwargs):
            illegal = set(kwargs) - {'parallel'}
            raise ValueError('Illegal arguments when reading from a file: '
                             f'{illegal}')
        atoms, dft, params, _ = read_gpw(filename,
                                         log=log,
                                         parallel=parallel)
        return ASECalculator(params,
                             log=log, dft=dft, atoms=atoms)

    params = InputParameters(kwargs)
    write_header(log, params)
    return ASECalculator(params, log=log)


LOGO = """\
  ___ ___ ___ _ _ _
 |   |   |_  | | | |
 | | | | | . | | | |
 |__ |  _|___|_____| - {version}
 |___|_|
"""


def write_header(log, params):
    from gpaw.io.logger import write_header as header
    log(LOGO.format(version=__version__))
    header(log, log.comm)
    log('---')
    with log.indent('input parameters:'):
        log(**dict(params.items()))


def compare_atoms(a1: Atoms, a2: Atoms) -> set[str]:
    if a1 is a2:
        return set()

    if len(a1.numbers) != len(a2.numbers) or (a1.numbers != a2.numbers).any():
        return {'numbers'}

    if (a1.pbc != a2.pbc).any():
        return {'pbc'}

    if abs(a1.cell - a2.cell).max() > 0.0:
        return {'cell'}

    if abs(a1.positions - a2.positions).max() > 0.0:
        return {'positions'}

    return set()


class ASECalculator:
    """This is the ASE-calculator frontend for doing a GPAW calculation."""

    name = 'gpaw'

    def __init__(self,
                 params: InputParameters,
                 *,
                 log: Logger,
                 dft: DFTCalculation | None = None,
                 atoms: Atoms | None = None):
        self.params = params
        self.log = log
        self.comm = log.comm
        self._dft = dft
        self._atoms = atoms
        self.timer = Timer()

    @property
    def dft(self) -> DFTCalculation:
        if self._dft is None:
            raise AttributeError
        return self._dft

    @property
    def atoms(self) -> Atoms:
        if self._atoms is None:
            raise AttributeError
        return self._atoms

    def __repr__(self):
        params = []
        for key, value in self.params.items():
            val = repr(value)
            if len(val) > 40:
                val = '...'
            params.append((key, val))
        p = ', '.join(f'{key}: {val}' for key, val in params)
        return f'ASECalculator({p})'

    def iconverge(self, atoms: Atoms | None):
        """Iterate to self-consistent solution.

        Will also calculate "cheap" properties: energy, magnetic moments
        and dipole moment.
        """
        if atoms is None:
            atoms = self.atoms
        else:
            synchronize_atoms(atoms, self.comm)

        converged = True

        if self._dft is not None:
            changes = compare_atoms(self.atoms, atoms)
            if changes & {'numbers', 'pbc', 'cell'}:
                if 'numbers' not in changes:
                    # Remember magmoms if there are any:
                    magmom_a = self.dft.results.get('magmoms')
                    if magmom_a is not None and magmom_a.any():
                        atoms = atoms.copy()
                        assert atoms is not None  # MYPY: why is this needed?
                        atoms.set_initial_magnetic_moments(magmom_a)

                if changes & {'numbers', 'pbc'}:
                    self._dft = None  # start from scratch
                else:
                    try:
                        self.create_new_calculation_from_old(atoms)
                    except ReuseWaveFunctionsError:
                        self._dft = None  # start from scratch
                    else:
                        converged = False
                        changes = set()

        if self._dft is None:
            self.create_new_calculation(atoms)
            converged = False
        elif changes:
            self.move_atoms(atoms)
            converged = False

        if converged:
            return

        with self.timer('SCF'):
            for ctx in self.dft.iconverge(
                    calculate_forces=self._calculate_forces):
                yield ctx

        self.log(f'Converged in {ctx.niter} steps')

        # Calculate all the cheap things:
        self.dft.energies()
        self.dft.dipole()
        self.dft.magmoms()

        self.dft.write_converged()

    def calculate_property(self,
                           atoms: Atoms | None,
                           prop: str) -> Any:
        """Calculate (if not already calculated) a property.

        The ``prop`` string must be one of

        * energy
        * forces
        * stress
        * magmom
        * magmoms
        * dipole
        """
        for _ in self.iconverge(atoms):
            pass

        if prop == 'forces':
            with self.timer('Forces'):
                self.dft.forces()
        elif prop == 'stress':
            with self.timer('Stress'):
                self.dft.stress()
        elif prop not in self.dft.results:
            raise KeyError('Unknown property:', prop)

        return self.dft.results[prop] * units[prop]

    def get_property(self,
                     name: str,
                     atoms: Atoms | None = None,
                     allow_calculation: bool = True) -> Any:
        if not allow_calculation and name not in self.dft.results:
            return None
        if atoms is None:
            atoms = self.atoms
        return self.calculate_property(atoms, name)

    @property
    def results(self):
        if self._dft is None:
            return {}
        return {name: value * units[name]
                for name, value in self.dft.results.items()}

    @trace
    def create_new_calculation(self, atoms: Atoms) -> None:
        with self.timer('Init'):
            self._dft = DFTCalculation.from_parameters(
                atoms, self.params, self.comm, self.log)
        self._atoms = atoms.copy()

    def create_new_calculation_from_old(self, atoms: Atoms) -> None:
        with self.timer('Morph'):
            self._dft = self.dft.new(
                atoms, self.params, self.log)
        self._atoms = atoms.copy()

    def move_atoms(self, atoms):
        with self.timer('Move'):
            self._dft = self.dft.move_atoms(atoms)
        self._atoms = atoms.copy()

    def _calculate_forces(self) -> Array2D:  # units: Ha/Bohr
        """Helper method for force-convergence criterium."""
        with self.timer('Forces'):
            self.dft.forces(silent=True)
        return self.dft.results['forces'].copy()

    def __del__(self):
        self.log('---')
        self.timer.write(self.log)
        try:
            mib = maxrss() / 1024**2
            self.log(f'\nMax RSS: {mib:.3f}  # MiB')
        except NameError:
            pass

    def get_potential_energy(self,
                             atoms: Atoms | None = None,
                             force_consistent: bool = False) -> float:
        return self.calculate_property(atoms,
                                       'free_energy' if force_consistent else
                                       'energy')

    @trace
    def get_forces(self, atoms: Atoms | None = None) -> Array2D:
        return self.calculate_property(atoms, 'forces')

    @trace
    def get_stress(self, atoms: Atoms | None = None) -> Array1D:
        return self.calculate_property(atoms, 'stress')

    def get_dipole_moment(self, atoms: Atoms | None = None) -> Array1D:
        return self.calculate_property(atoms, 'dipole')

    def get_magnetic_moment(self, atoms: Atoms | None = None) -> float:
        return self.calculate_property(atoms, 'magmom')

    def get_magnetic_moments(self, atoms: Atoms | None = None) -> Array1D:
        return self.calculate_property(atoms, 'magmoms')

    def get_non_collinear_magnetic_moment(self,
                                          atoms: Atoms | None = None
                                          ) -> Array1D:
        return self.calculate_property(atoms, 'non_collinear_magmom')

    def get_non_collinear_magnetic_moments(self,
                                           atoms: Atoms | None = None
                                           ) -> Array2D:
        return self.calculate_property(atoms, 'non_collinear_magmoms')

    def write(self, filename, mode=''):
        """Write calculator object to a file.

        Parameters
        ----------
        filename:
            File to be written
        mode:
            Write mode. Use ``mode='all'``
            to include wave functions in the file.
        """
        self.log(f'# Writing to {filename} (mode={mode!r})\n')

        write_gpw(filename, self.atoms, self.params,
                  self.dft, skip_wfs=mode != 'all')

    # Old API:

    implemented_properties = ['energy', 'free_energy',
                              'forces', 'stress',
                              'dipole', 'magmom', 'magmoms']

    def icalculate(self, atoms, system_changes=None):
        yield from self.iconverge(atoms)

    def new(self, **kwargs) -> ASECalculator:
        kwargs = {**dict(self.params.items()), **kwargs}
        return GPAW(**kwargs)

    def get_pseudo_wave_function(self, band, kpt=0, spin=None,
                                 periodic=False,
                                 broadcast=True) -> Array3D:
        state = self.dft.state
        collinear = state.ibzwfs.collinear
        if collinear:
            if spin is None:
                spin = 0
        else:
            assert spin is None
        wfs = state.ibzwfs.get_wfs(spin=spin if collinear else 0,
                                   kpt=kpt,
                                   n1=band, n2=band + 1)
        if wfs is not None:
            basis = getattr(self.dft.scf_loop.hamiltonian,
                            'basis', None)
            grid = state.density.nt_sR.desc.new(comm=None)
            if collinear:
                wfs = wfs.to_uniform_grid_wave_functions(grid, basis)
                psit_R = wfs.psit_nX[0]
            else:
                psit_sG = wfs.psit_nX[0]
                grid = grid.new(kpt=psit_sG.desc.kpt_c,
                                dtype=psit_sG.desc.dtype)
                psit_R = psit_sG.ifft(grid=grid)
            if not psit_R.desc.pbc.all():
                psit_R = psit_R.to_pbc_grid()
            if periodic:
                psit_R.multiply_by_eikr(-psit_R.desc.kpt_c)
            array_R = psit_R.data * Bohr**-1.5
        else:
            array_R = None
        if broadcast:
            array_R = bcast(array_R, 0, self.dft.comm)
        return array_R

    def get_atoms(self):
        atoms = self.atoms.copy()
        atoms.calc = self
        return atoms

    def get_fermi_level(self) -> float:
        return self.dft.state.ibzwfs.fermi_level * Ha

    def get_fermi_levels(self) -> Array1D:
        state = self.dft.state
        fl = state.ibzwfs.fermi_levels
        assert fl is not None and len(fl) == 2
        return fl * Ha

    def get_homo_lumo(self, spin: int = None) -> Array1D:
        state = self.dft.state
        return state.ibzwfs.get_homo_lumo(spin) * Ha

    def get_number_of_electrons(self):
        state = self.dft.state
        return state.ibzwfs.nelectrons

    def get_number_of_bands(self):
        state = self.dft.state
        return state.ibzwfs.nbands

    def get_number_of_grid_points(self):
        return self.dft.state.density.nt_sR.desc.size

    def get_effective_potential(self, spin=0):
        assert spin == 0
        vt_R = self.dft.state.potential.vt_sR[spin]
        return vt_R.to_pbc_grid().gather(broadcast=True).data * Ha

    def get_electrostatic_potential(self):
        density = self.dft.state.density
        potential, _ = self.dft.pot_calc.calculate(density)
        vHt_x = potential.vHt_x
        if isinstance(vHt_x, UGArray):
            return vHt_x.gather(broadcast=True).to_pbc_grid().data * Ha

        grid = self.dft.pot_calc.fine_grid
        return vHt_x.ifft(grid=grid).gather(broadcast=True).data * Ha

    def get_atomic_electrostatic_potentials(self):
        return self.dft.electrostatic_potential().atomic_potentials()

    def get_electrostatic_corrections(self):
        return self.dft.electrostatic_potential().atomic_corrections()

    def get_pseudo_density(self,
                           spin=None,
                           gridrefinement=1,
                           broadcast=True) -> Array3D:
        assert spin is None
        nt_sr = self.dft.densities().pseudo_densities(
            grid_refinement=gridrefinement)
        return nt_sr.gather(broadcast=broadcast).data.sum(0)

    def get_all_electron_density(self,
                                 spin=None,
                                 gridrefinement=1,
                                 broadcast=True,
                                 skip_core=False):
        n_sr = self.dft.densities().all_electron_densities(
            grid_refinement=gridrefinement,
            skip_core=skip_core)
        if spin is None:
            return n_sr.gather(broadcast=broadcast).data.sum(0)
        return n_sr[spin].gather(broadcast=broadcast).data

    def get_eigenvalues(self, kpt=0, spin=0, broadcast=True):
        state = self.dft.state
        eig_n = state.ibzwfs.get_eigs_and_occs(k=kpt, s=spin)[0] * Ha
        if broadcast:
            if self.comm.rank != 0:
                eig_n = np.empty(state.ibzwfs.nbands)
            self.comm.broadcast(eig_n, 0)
        return eig_n

    def get_occupation_numbers(self, kpt=0, spin=0, broadcast=True):
        state = self.dft.state
        weight = state.ibzwfs.ibz.weight_k[kpt] * state.ibzwfs.spin_degeneracy
        occ_n = state.ibzwfs.get_eigs_and_occs(k=kpt, s=spin)[1] * weight
        if broadcast:
            if self.comm.rank != 0:
                occ_n = np.empty(state.ibzwfs.nbands)
            self.comm.broadcast(occ_n, 0)
        return occ_n

    def get_reference_energy(self):
        return self.dft.setups.Eref * Ha

    def get_number_of_iterations(self):
        return self.dft.scf_loop.niter

    def get_bz_k_points(self):
        state = self.dft.state
        return state.ibzwfs.ibz.bz.kpt_Kc.copy()

    def get_ibz_k_points(self):
        state = self.dft.state
        return state.ibzwfs.ibz.kpt_kc.copy()

    def get_orbital_magnetic_moments(self):
        """Return the orbital magnetic moment vector for each atom."""
        state = self.dft.state
        if state.density.collinear:
            raise CalculationModeError(
                'Calculator is in collinear mode. '
                'Collinear calculations require spin–orbit '
                'coupling for nonzero orbital magnetic moments.')
        if not self.params.soc:
            warnings.warn('Non-collinear calculation was performed '
                          'without spin–orbit coupling. Orbital '
                          'magnetic moments may not be accurate.')
        return state.density.calculate_orbital_magnetic_moments()

    def calculate(self, atoms, properties=None, system_changes=None):
        if properties is None:
            properties = ['energy']

        for name in properties:
            self.calculate_property(atoms, name)
        # self.get_potential_energy(atoms)

    @cached_property
    def wfs(self):
        from gpaw.new.backwards_compatibility import FakeWFS
        return FakeWFS(self.dft, self.atoms)

    @property
    def density(self):
        from gpaw.new.backwards_compatibility import FakeDensity
        return FakeDensity(self.dft)

    @property
    def hamiltonian(self):
        from gpaw.new.backwards_compatibility import FakeHamiltonian
        return FakeHamiltonian(self.dft)

    @property
    def spos_ac(self):
        return self.atoms.get_scaled_positions()

    @property
    def world(self):
        return self.comm

    @property
    def setups(self):
        return self.dft.setups

    @property
    def initialized(self):
        return self._dft is not None

    def get_xc_functional(self):
        return self.dft.pot_calc.xc.name

    def get_xc_difference(self, xcparams):
        """Calculate non-selfconsistent XC-energy difference."""
        dft = self.dft
        pot_calc = dft.pot_calc
        state = dft.state
        density = dft.state.density
        xc = create_functional(xcparams, pot_calc.fine_grid)
        if xc.type == 'MGGA' and density.taut_sR is None:
            state.ibzwfs.make_sure_wfs_are_read_from_gpw_file()
            if isinstance(state.ibzwfs.wfs_qs[0][0].psit_nX, SimpleNamespace):
                params = InputParameters(dict(self.params.items()))
                builder = create_builder(self.atoms, params, self.comm)
                basis_set = builder.create_basis_set()
                ibzwfs = builder.create_ibz_wave_functions(
                    basis_set, state.potential, log=dft.log)
                ibzwfs.fermi_levels = state.ibzwfs.fermi_levels
                state.ibzwfs = ibzwfs
                dft.scf_loop.update_density_and_potential = False
                dft.converge()
            density.update_ked(state.ibzwfs)
        exct = pot_calc.calculate_non_selfconsistent_exc(
            xc, density.nt_sR, density.taut_sR)
        dexc = 0.0
        for a, D_sii in state.density.D_asii.items():
            setup = self.setups[a]
            dexc += xc.calculate_paw_correction(
                setup, np.array([pack_density(D_ii) for D_ii in D_sii.real]))
        dexc = state.ibzwfs.domain_comm.sum_scalar(dexc)
        return (exct + dexc - state.potential.energies['xc']) * Ha

    def diagonalize_full_hamiltonian(self,
                                     nbands: int | None = None,
                                     scalapack=None,
                                     expert: bool | None = None) -> None:
        if expert is not None:
            warnings.warn('Ignoring deprecated "expert" argument',
                          DeprecationWarning)
        state = self.dft.state

        if nbands is None:
            nbands = min(wfs.array_shape(global_shape=True)[0]
                         for wfs in self.dft.state.ibzwfs)
            nbands = self.dft.state.ibzwfs.kpt_comm.min_scalar(nbands)
            assert isinstance(nbands, int)

        self.dft.scf_loop.occ_calc._set_nbands(nbands)
        ibzwfs = diagonalize(state.potential,
                             state.ibzwfs,
                             self.dft.scf_loop.occ_calc,
                             nbands,
                             self.dft.pot_calc.xc)
        self.dft.state = DFTState(ibzwfs,
                                  state.density,
                                  state.potential)
        nbands = ibzwfs.nbands
        self.params.nbands = nbands
        self.params.keys.append('nbands')

    def gs_adapter(self):
        from gpaw.response.groundstate import ResponseGroundStateAdapter
        return ResponseGroundStateAdapter(self)

    def fixed_density(self, txt='-', **kwargs):
        kwargs = {**dict(self.params.items()), **kwargs}
        params = InputParameters(kwargs)
        log = Logger(txt, self.comm)
        builder = create_builder(self.atoms, params, self.comm)
        basis_set = builder.create_basis_set()
        state = self.dft.state
        comm1 = state.ibzwfs.kpt_band_comm
        comm2 = builder.communicators['D']
        potential = state.potential.redist(
            builder.grid,
            builder.electrostatic_potential_desc,
            builder.atomdist,
            comm1, comm2)
        density = state.density.redist(builder.grid,
                                       builder.interpolation_desc,
                                       builder.atomdist,
                                       comm1, comm2)
        ibzwfs = builder.create_ibz_wave_functions(basis_set, potential,
                                                   log=log)
        ibzwfs.fermi_levels = state.ibzwfs.fermi_levels
        state = DFTState(ibzwfs, density, potential)
        scf_loop = builder.create_scf_loop()
        scf_loop.update_density_and_potential = False

        dft = DFTCalculation(
            state,
            builder.setups,
            scf_loop,
            SimpleNamespace(fracpos_ac=self.dft.fracpos_ac,
                            poisson_solver=None),
            log)

        dft.converge()

        return ASECalculator(params,
                             log=log,
                             dft=dft,
                             atoms=self.atoms)

    def initialize(self, atoms):
        self.create_new_calculation(atoms)

    def converge_wave_functions(self):
        self.dft.state.ibzwfs.make_sure_wfs_are_read_from_gpw_file()

    def get_number_of_spins(self):
        return self.dft.state.density.ndensities

    @property
    def parameters(self):
        return self.params

    def dos(self,
            soc: bool = False,
            theta: float = 0.0,  # degrees
            phi: float = 0.0,  # degrees
            shift_fermi_level: bool = True) -> DOSCalculator:
        """Create DOS-calculator.

        Default is to ``shift_fermi_level`` to 0.0 eV.  For ``soc=True``,
        angles can be given in degrees.
        """
        return DOSCalculator.from_calculator(
            self, soc=soc,
            theta=theta, phi=phi,
            shift_fermi_level=shift_fermi_level)

    def band_structure(self):
        """Create band-structure object for plotting."""
        from ase.spectrum.band_structure import get_band_structure
        return get_band_structure(calc=self)

    @property
    def symmetry(self):
        return self.dft.state.ibzwfs.ibz.symmetries.symmetry
