from __future__ import annotations

from functools import cached_property
from typing import Any, Union

import numpy as np
from ase import Atoms
from ase.units import Bohr, Ha

from gpaw.core import UGDesc
from gpaw.core.atom_arrays import AtomDistribution
from gpaw.densities import Densities
from gpaw.electrostatic_potential import ElectrostaticPotential
from gpaw.gpu import as_np
from gpaw.mpi import broadcast_float, world
from gpaw.new import trace, zips
from gpaw.new.density import Density
from gpaw.new.ibzwfs import IBZWaveFunctions
from gpaw.new.input_parameters import InputParameters
from gpaw.new.logger import Logger
from gpaw.new.potential import Potential
from gpaw.new.scf import SCFLoop
from gpaw.setup import Setups
from gpaw.typing import Array1D, Array2D
from gpaw.utilities import (check_atoms_too_close,
                            check_atoms_too_close_to_boundary)
from gpaw.utilities.partition import AtomPartition


class ReuseWaveFunctionsError(Exception):
    """Reusing the old wave functions after cell-change failed.

    Most likekly, the number of k-points changed.
    """


class NonsenseError(Exception):
    """Operation doesn't make sense."""


class CalculationModeError(Exception):
    """Calculation mode does not match what is expected from a given method.

    For example, if a method only works in collinear mode and receives a
    calculator in non-collinear mode, this exception should be raised.
    """


units = {'energy': Ha,
         'free_energy': Ha,
         'forces': Ha / Bohr,
         'stress': Ha / Bohr**3,
         'dipole': Bohr,
         'magmom': 1.0,
         'magmoms': 1.0,
         'non_collinear_magmom': 1.0,
         'non_collinear_magmoms': 1.0}


class DFTState:
    def __init__(self,
                 ibzwfs: IBZWaveFunctions,
                 density: Density,
                 potential: Potential):
        """State of a Kohn-Sham calculation."""
        self.ibzwfs = ibzwfs
        self.density = density
        self.potential = potential

    def __repr__(self):
        return (f'DFTState({self.ibzwfs!r}, '
                f'{self.density!r}, {self.potential!r})')

    def __str__(self):
        return f'{self.ibzwfs}\n{self.density}\n{self.potential}'

    def move(self, fracpos_ac, atomdist):
        self.ibzwfs.move(fracpos_ac, atomdist)
        self.density.move(fracpos_ac, atomdist)
        self.potential.move(atomdist)


class DFTCalculation:
    def __init__(self,
                 state: DFTState,
                 setups: Setups,
                 scf_loop: SCFLoop,
                 pot_calc,
                 log: Logger):
        self.state = state
        self.setups = setups
        self.scf_loop = scf_loop
        self.pot_calc = pot_calc
        self.log = log
        self.comm = log.comm

        self.results: dict[str, Any] = {}
        self.fracpos_ac = self.pot_calc.fracpos_ac

    @classmethod
    def from_parameters(cls,
                        atoms: Atoms,
                        params: Union[dict, InputParameters],
                        comm=None,
                        log=None,
                        builder=None) -> DFTCalculation:
        """Create DFTCalculation object from parameters and atoms."""
        from gpaw.new.builder import builder as create_builder

        check_atoms_too_close(atoms)
        check_atoms_too_close_to_boundary(atoms)

        if params is None:
            params = {}
        if isinstance(params, dict):
            params = InputParameters(params)

        if not isinstance(log, Logger):
            log = Logger(log, comm or world)

        builder = builder or create_builder(atoms, params, log.comm)

        basis_set = builder.create_basis_set()

        density = builder.density_from_superposition(basis_set)
        density.normalize()

        # The SCF-loop has a Hamiltonian that has an fft-plan that is
        # cached for later use, so best to create the SCF-loop first
        # FIX this!
        scf_loop = builder.create_scf_loop()

        pot_calc = builder.create_potential_calculator()
        potential, _ = pot_calc.calculate_without_orbitals(
            density, kpt_band_comm=builder.communicators['D'])
        ibzwfs = builder.create_ibz_wave_functions(
            basis_set, potential, log=log)

        if ibzwfs.wfs_qs[0][0]._eig_n is not None:
            ibzwfs.calculate_occs(scf_loop.occ_calc)

        state = DFTState(ibzwfs, density, potential)

        write_atoms(atoms, builder.initial_magmom_av, builder.grid, log)
        log(state)
        log(builder.setups)
        log(scf_loop)
        log(pot_calc)

        return cls(state, builder.setups, scf_loop, pot_calc, log)

    def move_atoms(self, atoms) -> DFTCalculation:
        check_atoms_too_close(atoms)

        self.fracpos_ac = np.ascontiguousarray(atoms.get_scaled_positions())
        self.comm.broadcast(self.fracpos_ac, 0)

        atomdist = self.state.density.D_asii.layout.atomdist
        grid = self.state.density.nt_sR.desc
        rank_a = grid.ranks_from_fractional_positions(self.fracpos_ac)
        atomdist = AtomDistribution(rank_a, atomdist.comm)

        self.pot_calc.move(self.fracpos_ac, atomdist)
        self.state.move(self.fracpos_ac, atomdist)

        mm_av = self.results['non_collinear_magmoms']
        write_atoms(atoms, mm_av, self.state.density.nt_sR.desc, self.log)

        self.results = {}

        return self

    def iconverge(self, convergence=None, maxiter=None, calculate_forces=None):
        self.state.ibzwfs.make_sure_wfs_are_read_from_gpw_file()
        yield from self.scf_loop.iterate(self.state,
                                         self.pot_calc,
                                         convergence,
                                         maxiter,
                                         calculate_forces,
                                         log=self.log)

    @trace
    def converge(self,
                 convergence=None,
                 maxiter=None,
                 steps=99999999999999999,
                 calculate_forces=None):
        """Converge to self-consistent solution of Kohn-Sham equation."""
        for step, _ in enumerate(self.iconverge(convergence,
                                                maxiter,
                                                calculate_forces),
                                 start=1):
            if step == steps:
                break
        else:  # no break
            self.log(scf_steps=step)

    def energies(self):
        energies = combine_energies(self.state.potential, self.state.ibzwfs)

        self.log('Energy contributions relative to reference atoms:',
                 f'(reference = {self.setups.Eref * Ha:.6f})\n')

        for name, e in energies.items():
            if not name.startswith('total') and name != 'stress':
                self.log(f'{name + ":":10}   {e * Ha:14.6f}')
        total_free = energies['total_free']
        total_extrapolated = energies['total_extrapolated']
        self.log('----------------------------')
        self.log(f'Free energy: {total_free * Ha:14.6f}')
        self.log(f'Extrapolated:{total_extrapolated * Ha:14.6f}\n')

        total_free = broadcast_float(total_free, self.comm)
        total_extrapolated = broadcast_float(total_extrapolated, self.comm)

        self.results['free_energy'] = total_free
        self.results['energy'] = total_extrapolated

    def dipole(self):
        if 'dipole' in self.results:
            return
        dipole_v = self.state.density.calculate_dipole_moment(self.fracpos_ac)
        x, y, z = dipole_v * Bohr
        self.log(f'dipole moment: [{x:.6f}, {y:.6f}, {z:.6f}]  # |e|*Ang\n')
        self.results['dipole'] = dipole_v

    def magmoms(self) -> tuple[Array1D, Array2D]:
        mm_v, mm_av = self.state.density.calculate_magnetic_moments()
        self.results['magmom'] = mm_v[2]
        self.results['magmoms'] = mm_av[:, 2].copy()
        self.results['non_collinear_magmoms'] = mm_av
        self.results['non_collinear_magmom'] = mm_v

        if self.state.density.ncomponents > 1:
            x, y, z = mm_v
            self.log(f'total magnetic moment: [{x:.6f}, {y:.6f}, {z:.6f}]\n')
            self.log('local magnetic moments: [')
            for a, (setup, m_v) in enumerate(zips(self.setups, mm_av)):
                x, y, z = m_v
                c = ',' if a < len(mm_av) - 1 else ']'
                self.log(f'  [{x:9.6f}, {y:9.6f}, {z:9.6f}]{c}'
                         f'  # {setup.symbol:2} {a}')
            self.log()
        return mm_v, mm_av

    def forces(self, silent=False):
        """Calculate atomic forces."""
        if 'forces' not in self.results or silent:
            self._calculate_forces()

            if silent:
                return
            self.log('\nforces: [  # eV/Ang')
            F_av = self.results['forces'] * (Ha / Bohr)
            for a, setup in enumerate(self.setups):
                x, y, z = F_av[a]
                c = ',' if a < len(F_av) - 1 else ']'
                self.log(f'  [{x:9.3f}, {y:9.3f}, {z:9.3f}]{c}'
                         f'  # {setup.symbol:2} {a}')

    def _calculate_forces(self):
        xc = self.pot_calc.xc
        assert not hasattr(xc.xc, 'setup_force_corrections')

        # Force from projector functions (and basis set):
        F_av = self.state.ibzwfs.forces(self.state.potential)

        pot_calc = self.pot_calc
        Fcc_avL, Fnct_av, Ftauct_av, Fvbar_av = pot_calc.force_contributions(
            self.state)

        # Force from compensation charges:
        ccc_aL = \
            self.state.density.calculate_compensation_charge_coefficients()
        for a, dF_vL in Fcc_avL.items():
            F_av[a] += dF_vL @ ccc_aL[a]

        # Force from smooth core charge:
        for a, dF_v in Fnct_av.items():
            F_av[a] += dF_v[:, 0]

        if Ftauct_av is not None:
            # Force from smooth core ked:
            for a, dF_v in Ftauct_av.items():
                F_av[a] += dF_v[:, 0]

        # Force from zero potential:
        for a, dF_v in Fvbar_av.items():
            F_av[a] += dF_v[:, 0]

        F_av = as_np(F_av)

        domain_comm = ccc_aL.layout.atomdist.comm
        domain_comm.sum(F_av)

        F_av = self.state.ibzwfs.ibz.symmetries.symmetrize_forces(F_av)
        self.comm.broadcast(F_av, 0)
        self.results['forces'] = F_av

    def stress(self) -> None:
        if 'stress' in self.results:
            return
        stress_vv = self.pot_calc.stress(self.state)
        self.log('\nstress tensor: [  # eV/Ang^3')
        for (x, y, z), c in zips(stress_vv * (Ha / Bohr**3), ',,]'):
            self.log(f'  [{x:13.6f}, {y:13.6f}, {z:13.6f}]{c}')
        self.results['stress'] = stress_vv.flat[[0, 4, 8, 5, 2, 1]]

    def write_converged(self) -> None:
        self.state.ibzwfs.write_summary(self.log)
        vacuum_level = self.state.potential.get_vacuum_level()
        if not np.isnan(vacuum_level):
            self.log(f'vacuum-level: {vacuum_level:.3f}  # V')
            try:
                wf1, wf2 = self.workfunctions(vacuum_level=vacuum_level)
            except NonsenseError:
                pass
            else:
                self.log(f'Workfunctions: {wf1:.3f}, {wf2:.3f}  # eV')
        self.log.fd.flush()

    def workfunctions(self,
                      *,
                      vacuum_level: float | None = None
                      ) -> tuple[float, float]:
        if vacuum_level is None:
            vacuum_level = self.state.potential.get_vacuum_level()
        if np.isnan(vacuum_level):
            raise NonsenseError('No vacuum')
        try:
            correction = self.pot_calc.poisson_solver.dipole_layer_correction()
        except NotImplementedError:
            raise NonsenseError('No dipole layer')
        correction *= Ha
        fermi_level = self.state.ibzwfs.fermi_level * Ha
        wf = vacuum_level - fermi_level
        return wf - correction, wf + correction

    def electrostatic_potential(self) -> ElectrostaticPotential:
        return ElectrostaticPotential.from_calculation(self)

    def densities(self) -> Densities:
        return Densities.from_calculation(self)

    @cached_property
    def _atom_partition(self):
        # Backwards compatibility helper
        atomdist = self.state.density.D_asii.layout.atomdist
        return AtomPartition(atomdist.comm, atomdist.rank_a)

    def new(self,
            atoms: Atoms,
            params: InputParameters,
            log=None) -> DFTCalculation:
        """Create new DFTCalculation object."""
        from gpaw.new.builder import builder as create_builder

        if params.mode['name'] != 'pw':
            raise ReuseWaveFunctionsError

        ibzwfs = self.state.ibzwfs
        if ibzwfs.domain_comm.size != 1:
            raise ReuseWaveFunctionsError

        if not self.state.density.nt_sR.desc.pbc_c.all():
            raise ReuseWaveFunctionsError

        check_atoms_too_close(atoms)
        check_atoms_too_close_to_boundary(atoms)

        builder = create_builder(atoms, params, self.comm)

        kpt_kc = builder.ibz.kpt_kc
        old_kpt_kc = ibzwfs.ibz.kpt_kc
        if len(kpt_kc) != len(old_kpt_kc):
            raise ReuseWaveFunctionsError
        if abs(kpt_kc - old_kpt_kc).max() > 1e-9:
            raise ReuseWaveFunctionsError

        log('# Interpolating wave functions to new cell')

        density = self.state.density.new(builder.grid,
                                         builder.interpolation_desc,
                                         builder.fracpos_ac,
                                         builder.atomdist)
        density.normalize()

        # Make sure all have exactly the same density.
        # Not quite sure it is needed???
        # At the moment we skip it on GPU's because it doesn't
        # work!
        if density.nt_sR.xp is np:
            self.comm.broadcast(density.nt_sR.data, 0)

        scf_loop = builder.create_scf_loop()
        pot_calc = builder.create_potential_calculator()
        potential, _ = pot_calc.calculate(density)

        old_ibzwfs = ibzwfs

        def create_wfs(spin, q, k, kpt_c, weight):
            wfs = old_ibzwfs.wfs_qs[q][spin]
            return wfs.morph(
                builder.wf_desc,
                builder.fracpos_ac,
                builder.atomdist)

        ibzwfs = ibzwfs.create(
            ibz=builder.ibz,
            nelectrons=old_ibzwfs.nelectrons,
            ncomponents=old_ibzwfs.ncomponents,
            create_wfs_func=create_wfs,
            kpt_comm=old_ibzwfs.kpt_comm,
            kpt_band_comm=old_ibzwfs.kpt_band_comm,
            comm=self.comm)

        state = DFTState(ibzwfs, density, potential)

        write_atoms(atoms, builder.initial_magmom_av, builder.grid, log)
        log(state)
        log(builder.setups)
        log(scf_loop)
        log(pot_calc)

        return DFTCalculation(
            state, builder.setups, scf_loop, pot_calc, log)


def combine_energies(potential: Potential,
                     ibzwfs: IBZWaveFunctions) -> dict[str, float]:
    """Add up energy contributions."""
    energies = potential.energies.copy()
    energies.pop('stress', 0.0)
    energies['kinetic'] += ibzwfs.energies['band']
    energies['kinetic'] += ibzwfs.energies.get('exx_kinetic', 0.0)
    energies['xc'] += (ibzwfs.energies.get('exx_vv', 0.0) +
                       ibzwfs.energies.get('exx_vc', 0.0) +
                       ibzwfs.energies.get('exx_cc', 0.0))
    energies['entropy'] = ibzwfs.energies['entropy']
    energies['total_free'] = sum(energies.values())
    energies['total_extrapolated'] = (energies['total_free'] +
                                      ibzwfs.energies['extrapolation'])
    return energies


def write_atoms(atoms: Atoms,
                magmom_av: Array2D,
                grid: UGDesc,
                log) -> None:
    from gpaw.output import print_cell, print_positions
    print_positions(atoms, log, magmom_av)
    print_cell(grid._gd, grid.pbc, log)
