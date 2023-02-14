from __future__ import annotations

from math import pi
from types import SimpleNamespace

import numpy as np
from ase.data import atomic_numbers, covalent_radii
from ase.neighborlist import neighbor_list
from ase.units import Bohr, Ha

from gpaw.core.arrays import DistributedArrays
from gpaw.core.atom_arrays import AtomArraysLayout
from gpaw.core.domain import Domain
from gpaw.core.matrix import Matrix
from gpaw.lcao.tci import TCIExpansions
from gpaw.lfc import BasisFunctions
from gpaw.mpi import MPIComm, serial_comm
from gpaw.new import zip
from gpaw.new.calculation import DFTState
from gpaw.new.lcao.builder import LCAODFTComponentsBuilder, create_lcao_ibzwfs
from gpaw.new.lcao.hamiltonian import CollinearHamiltonianMatrixCalculator
from gpaw.new.lcao.wave_functions import LCAOWaveFunctions
from gpaw.new.pot_calc import PotentialCalculator
from gpaw.setup import Setup
from gpaw.spline import Spline
from gpaw.utilities.timing import NullTimer
from gpaw.typing import Array3D


class TBHamiltonianMatrixCalculator(CollinearHamiltonianMatrixCalculator):
    def _calculate_potential_matrix(self,
                                    wfs: LCAOWaveFunctions,
                                    V_xMM: Array3D = None) -> Matrix:
        return wfs.V_MM


class TBHamiltonian:
    def __init__(self,
                 basis: BasisFunctions):
        self.basis = basis

    def apply(self):
        raise NotImplementedError

    def create_hamiltonian_matrix_calculator(
            self,
            state: DFTState) -> TBHamiltonianMatrixCalculator:
        dH_saii = [{a: dH_sii[s]
                    for a, dH_sii in state.potential.dH_asii.items()}
                   for s in range(state.density.ncomponents)]

        V_sxMM = [np.zeros(0) for _ in range(state.density.ncomponents)]

        return TBHamiltonianMatrixCalculator(V_sxMM, dH_saii, self.basis)


class NoGrid(Domain):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._gd = SimpleNamespace(
            get_grid_spacings=lambda: [0, 0, 0],
            cell_cv=self.cell_cv,
            pbc_c=self.pbc_c,
            N_c=[0, 0, 0],
            dv=0.0)
        self.size = (0, 0, 0)

    def empty(self, shape=(), comm=serial_comm, xp=None):
        return DummyFunctions(self, shape, comm)

    def ranks_from_fractional_positions(self, fracpos_ac):
        return np.zeros(len(fracpos_ac), int)


class DummyFunctions(DistributedArrays[NoGrid]):
    def __init__(self,
                 grid: NoGrid,
                 dims: int | tuple[int, ...] = (),
                 comm: MPIComm = serial_comm):
        DistributedArrays. __init__(self, dims, (),
                                    comm, grid.comm, None, np.nan,
                                    grid.dtype)
        self.desc = grid

    def integrate(self):
        return np.ones(self.dims)

    def new(self):
        return self

    def __getitem__(self, index):
        return DummyFunctions(self.desc, comm=self.comm)

    def moment(self):
        return np.zeros(3)


class PSCoreDensities:
    def __init__(self, grid, fracpos_ac):
        self.layout = AtomArraysLayout([1] * len(fracpos_ac),
                                       grid.comm)

    def to_uniform_grid(self, out, scale):
        pass


class TBPotentialCalculator(PotentialCalculator):
    def __init__(self,
                 xc,
                 setups,
                 nct_R,
                 atoms):
        super().__init__(xc, None, setups, nct_R,
                         atoms.get_scaled_positions())
        self.atoms = atoms.copy()
        self.force_av = None
        self.stress_vv = None

    def calculate_charges(self, vHt_r):
        return {a: np.zeros(9) for a, setup in enumerate(self.setups)}

    def _calculate(self, density, ibzwfs, vHt_r):
        vt_sR = density.nt_sR

        atoms = self.atoms
        energy, force_av, stress_vv = pairpot(atoms)
        energy /= Ha
        self.force_av = force_av * Bohr / Ha

        vol = abs(np.linalg.det(atoms.cell[atoms.pbc][:, atoms.pbc]))
        self.stress_vv = stress_vv / vol * Bohr**atoms.pbc.sum() / Ha

        return {'kinetic': 0.0,
                'coulomb': 0.0,
                'zero': 0.0,
                'xc': energy,
                'external': 0.0}, vt_sR, vHt_r

    def _move(self, fracpos_ac, ndensities):
        self.atoms.set_scaled_positions(fracpos_ac)
        self.force_av = None
        self.stress_vv = None

    def force_contributions(self, state):
        return {}, {}, {a: self.force_av[a:a + 1]
                        for a in state.density.D_asii.keys()}

    def stress_contribution(self, state):
        return self.stress_vv


class DummyXC:
    no_forces = False
    xc = None

    def calculate_paw_correction(self, setup, D_sp, dH_sp):
        return 0.0


class TBSCFLoop:
    def __init__(self, hamiltonian, occ_calc, eigensolver):
        self.hamiltonian = hamiltonian
        self.occ_calc = occ_calc
        self.eigensolver = eigensolver

    def iterate(self,
                state,
                pot_calc,
                convergence=None,
                maxiter=None,
                calculate_forces=None,
                log=None):
        self.eigensolver.iterate(state, self.hamiltonian)
        state.ibzwfs.calculate_occs(self.occ_calc)
        yield
        state.potential, state.vHt_x, _ = pot_calc.calculate(
            state.density, state.vHt_x)


class DummyBasis:
    def __init__(self, natoms):
        self.my_atom_indices = np.arange(natoms)

    def add_to_density(self, nt_sR, f_asi):
        pass

    def construct_density(self, rho_MM, nt_G, q):
        pass


class TBDFTComponentsBuilder(LCAODFTComponentsBuilder):
    def check_cell(self, cell):
        pass

    def create_uniform_grids(self):
        grid = NoGrid(
            self.atoms.cell.complete() / Bohr,
            self.atoms.pbc,
            dtype=self.dtype,
            comm=self.communicators['d'])
        return grid, grid

    def get_pseudo_core_densities(self):
        return PSCoreDensities(self.grid, self.fracpos_ac)

    def create_basis_set(self):
        self.basis = DummyBasis(len(self.atoms))
        return self.basis

    def create_hamiltonian_operator(self):
        return TBHamiltonian(self.basis)

    def create_potential_calculator(self):
        xc = DummyXC()
        return TBPotentialCalculator(xc, self.setups, self.nct_R, self.atoms)

    def create_scf_loop(self):
        occ_calc = self.create_occupation_number_calculator()
        hamiltonian = self.create_hamiltonian_operator()
        eigensolver = self.create_eigensolver(hamiltonian)
        return TBSCFLoop(hamiltonian, occ_calc, eigensolver)

    def create_ibz_wave_functions(self,
                                  basis: BasisFunctions,
                                  potential,
                                  coefficients=None):
        assert self.communicators['w'].size == 1

        ibzwfs, tciexpansions = create_lcao_ibzwfs(
            basis, potential,
            self.ibz, self.communicators, self.setups,
            self.fracpos_ac, self.grid, self.dtype,
            self.nbands, self.ncomponents, self.atomdist, self.nelectrons)

        vtphit: dict[Setup, list[Spline]] = {}

        for setup in self.setups.setups.values():
            try:
                vt_r = setup.vt_g
            except AttributeError:
                vt_r = calculate_pseudo_potential(setup, self.xc.xc)[0]

            vt_r[-1] = 0.0  # ???
            vt = setup.rgd.spline(vt_r, points=300)
            vtphit_j = []
            for phit in setup.basis_functions_J:
                rc = phit.get_cutoff()
                r_g = np.linspace(0, rc, 150)
                vt_g = vt.map(r_g) / (4 * pi)**0.5
                phit_g = phit.map(r_g)
                vtphit_j.append(Spline(phit.l, rc, vt_g * phit_g))
            vtphit[setup] = vtphit_j

        vtciexpansions = TCIExpansions([s.basis_functions_J
                                        for s in self.setups],
                                       [vtphit[s] for s in self.setups],
                                       tciexpansions.I_a)

        kpt_qc = np.array([wfs.kpt_c for wfs in ibzwfs])
        manytci = vtciexpansions.get_manytci_calculator(
            self.setups, self.grid._gd, self.fracpos_ac,
            kpt_qc, self.dtype, NullTimer())

        manytci.Pindices = manytci.Mindices
        my_atom_indices = basis.my_atom_indices

        for wfs, V_MM in zip(ibzwfs, manytci.P_qIM(my_atom_indices)):
            V_MM = V_MM.toarray()
            V_MM += V_MM.T.conj().copy()
            M1 = 0
            for m in manytci.Mindices.nm_a:
                M2 = M1 + m
                V_MM[M1:M2, M1:M2] *= 0.5
                M1 = M2
            wfs.V_MM = Matrix(M2, M2, data=V_MM)

        return ibzwfs


def pairpot(atoms):
    """Simple pair-potential for testing.

    >>> from ase import Atoms
    >>> r = covalent_radii[1]
    >>> atoms = Atoms('H2', [(0, 0, 0), (0, 0, 2 * r)])
    >>> e, f, s = pairpot(atoms)
    >>> print(f'{e:.6f} eV')
    -9.677419 eV
    >>> f
    array([[0., 0., 0.],
           [0., 0., 0.]])

    """
    radii = {}
    symbol_a = atoms.symbols
    for symbol in symbol_a:
        radii[symbol] = covalent_radii[atomic_numbers[symbol]]

    r0 = {}
    for s1, r1 in radii.items():
        for s2, r2 in radii.items():
            r0[(s1, s2)] = r1 + r2
    rcutmax = 2 * max(r0.values(), default=1.0)

    energy = 0.0
    force_av = np.zeros((len(atoms), 3))
    stress_vv = np.zeros((3, 3))

    for i, j, d, D_v in zip(*neighbor_list('ijdD', atoms, rcutmax)):
        d0 = r0[(symbol_a[i], symbol_a[j])]
        e0 = 6.0 / d0
        x = d0 / d
        if x > 0.5:
            energy += 0.5 * e0 * (-5 + x * (24 + x * (-36 + 16 * x)))
            f = -0.5 * e0 * (24 + x * (-72 + 48 * x)) * d0 / d**2
            F_v = D_v * f / d
            force_av[i] += F_v
            force_av[j] -= F_v
            # print(i, j, d, D_v, F_v)
            stress_vv += np.outer(F_v, D_v)

    return energy, force_av, stress_vv


def calculate_pseudo_potential(setup: Setup, xc):
    phit_jg = np.array(setup.data.phit_jg)
    rgd = setup.rgd

    # Density:
    nt_g = np.einsum('jg, j, jg -> g',
                     phit_jg, setup.f_j, phit_jg) / (4 * pi)
    nt_g += setup.data.nct_g * (1 / (4 * pi)**0.5)

    # XC:
    vt_g = rgd.zeros()
    xc.calculate_spherical(rgd, nt_g[np.newaxis], vt_g[np.newaxis])

    # Zero-potential:
    vt_g += setup.data.vbar_g / (4 * pi)**0.5

    # Coulomb:
    g_g = setup.ghat_l[0].map(rgd.r_g)
    Q = -rgd.integrate(nt_g) / rgd.integrate(g_g)
    rhot_g = nt_g + Q * g_g
    vHtr_g = rgd.poisson(rhot_g)

    W = rgd.integrate(g_g * vHtr_g, n=-1) / (4 * pi)**0.5

    vtr_g = vt_g * rgd.r_g + vHtr_g

    vtr_g[1:] /= rgd.r_g[1:]
    vtr_g[0] = vtr_g[1]

    return vtr_g * (4 * pi)**0.5, W


def poly():
    """Polynomium used for pair potential."""
    import matplotlib.pyplot as plt
    c = np.linalg.solve([[1, 0.5, 0.25, 0.125],
                         [1, 1, 1, 1],
                         [0, 1, 1, 0.75],
                         [0, 1, 2, 3]],
                        [0, -1, 0, 0])
    print(c)
    d = np.linspace(0.5, 2, 101)
    plt.plot(d, c[0] + c[1] / d + c[2] / d**2 + c[3] / d**3)
    plt.show()


if __name__ == '__main__':
    poly()
