from __future__ import annotations
import numpy as np
from functools import partial
from gpaw.core.arrays import DistributedArrays as DA
from gpaw.setup import Setups
from gpaw.typing import Array1D, Array2D, ArrayLike1D
from gpaw.new.brillouin import IBZ
from gpaw.mpi import MPIComm
from ase.units import Ha
from gpaw.new.density import Density
import _gpaw
from gpaw.core.atom_arrays import AtomArrays


class IBZWaveFunctions:
    def __init__(self,
                 ibz: IBZ,
                 ranks: ArrayLike1D,
                 kpt_comm: MPIComm,
                 mykpts: list[WaveFunctions],
                 nelectrons: int):
        self.ibz = ibz
        self.ranks = ranks
        self.kpt_comm = kpt_comm
        self.mykpts = mykpts
        self.nelectrons = nelectrons
        self.fermi_levels = None
        self.collinear = False
        self.spin_degeneracy = 2

        self.ibz_index_to_local_index = {}
        j = 0
        for i, rank in enumerate(ranks):
            if rank == kpt_comm.rank:
                self.ibz_index_to_local_index[i] = j
                j += 1

    def __iter__(self):
        for wfs in self.mykpts:
            yield wfs

    def __getitem__(self, n: int) -> WaveFunctions:
        return self.mykpts[n]

    @classmethod
    def from_random_numbers(cls,
                            ibz,
                            band_comm,
                            kpt_comm,
                            grid,
                            setups,
                            fracpos,
                            nbands: int,
                            nelectrons: float) -> IBZWaveFunctions:
        assert len(ibz) == 1
        ranks = [0]

        mykpts = []
        for kpt, weight, rank in zip(ibz.points, ibz.weights, ranks):
            if rank != kpt_comm.rank:
                continue
            basis = grid.new(kpt=kpt)
            wfs = WaveFunctions.from_random_numbers(basis, weight,
                                                    nbands, band_comm,
                                                    setups,
                                                    fracpos)
            mykpts.append(wfs)

        return cls(ibz, ranks, kpt_comm, mykpts, nelectrons)

    def orthonormalize(self, work_array=None):
        for wfs in self:
            wfs.orthonormalize(work_array)

    def calculate_occs(self, occ_calc, fixed_fermi_level=False):
        degeneracy = self.spin_degeneracy

        occs, fermi_levels, e_entropy = occ_calc.calculate(
            nelectrons=self.nelectrons / degeneracy,
            eigenvalues=[wfs.eigs * Ha for wfs in self],
            weights=[wfs.weight for wfs in self],
            fermi_levels_guess=(None
                                if self.fermi_levels is None else
                                self.fermi_levels * Ha))

        if not fixed_fermi_level or self.fermi_levels is None:
            self.fermi_levels = np.array(fermi_levels) / Ha

        for occsk, wfs in zip(occs, self):
            wfs._occs = occsk

        self.e_entropy = e_entropy * degeneracy / Ha
        e_band = 0.0
        for wfs in self:
            e_band += wfs.occs @ wfs.eigs * wfs.weight * degeneracy
        self.e_band = self.kpt_comm.sum(e_band)

    def calculate_density(self, out: Density) -> None:
        density = out
        density.density.data[:] = 0.0
        density.density_matrices.data[:] = 0.0
        for wfs in self:
            wfs.add_to_density(density.density, density.density_matrices)

    def get_eigs_and_occs(self, i):
        assert self.ranks[i] == self.kpt_comm.rank
        wfs = self.mykpts[self.ibz_index_to_local_index[i]]
        return wfs.eigs, wfs.occs

    def forces(self, dv: AtomArrays):
        F = np.zeros((len(dv.layout.shapes), 3))
        for wfs in self:
            wfs.force_contribution(dv, F)
        return F


class WaveFunctions:
    def __init__(self,
                 wave_functions: DA,
                 spin: int | None,
                 setups: Setups,
                 positions: Array2D,
                 weight: float = 1.0,
                 spin_degeneracy: int = 2):
        self.wave_functions = wave_functions
        self.spin = spin
        self.setups = setups
        self.weight = weight
        self.spin_degeneracy = spin_degeneracy

        self._projections = None
        self.projectors = setups.create_projectors(wave_functions.layout,
                                                   positions)
        self.orthonormalized = False

        self._eigs: Array1D | None = None
        self._occs: Array1D | None = None

    @property
    def eigs(self) -> Array1D:
        if self._eigs is None:
            raise ValueError
        return self._eigs

    @property
    def occs(self) -> Array1D:
        if self._occs is None:
            raise ValueError
        return self._occs

    @property
    def projections(self):
        if self._projections is None:
            self._projections = self.projectors.integrate(self.wave_functions)
        return self._projections

    @property
    def myeigs(self):
        assert self.wave_functions.comm.size == 1
        return self.eigs

    @property
    def myoccs(self):
        assert self.wave_functions.comm.size == 1
        return self.occs

    @classmethod
    def from_random_numbers(cls, basis, weight, nbands, band_comm, setups,
                            positions):
        wfs = basis.random(nbands, band_comm)
        return cls(wfs, 0, setups, positions)

    def add_to_density(self,
                       density,
                       density_matrices: AtomArrays) -> None:
        myoccs = self.weight * self.spin_degeneracy * self.myoccs
        for f, psit in zip(myoccs, self.wave_functions.data):
            # Same as density.data += f * abs(psit)**2, but much faster:
            _gpaw.add_to_density(f, psit, density.data[self.spin])

        for D, P in zip(density_matrices.values(), self.projections.values()):
            D[:, :, self.spin] += np.einsum('in, n, jn -> ij',
                                            P.conj(), myoccs, P)

    def orthonormalize(self, work_array=None):
        if self.orthonormalized:
            return
        wfs = self.wave_functions
        domain_comm = wfs.layout.comm

        projections = self.projections

        projections2 = projections.new()
        wfs2 = wfs.new(data=work_array)

        dS = partial(self.setups.overlap_correction, out=projections2)

        S = wfs.matrix_elements(wfs, domain_sum=False)
        projections.matrix_elements(projections, function=dS,
                                    domain_sum=False, out=S, add_to_out=True)
        domain_comm.sum(S.data, 0)
        if domain_comm.rank == 0:
            S.invcholesky()
        # S now contains the inverse of the Cholesky factorization
        domain_comm.broadcast(S.data, 0)
        # cc ??????

        S.multiply(wfs, out=wfs2)
        projections.matrix.multiply(S, opb='T', out=projections2)
        wfs.data[:] = wfs2.data
        projections.data[:] = projections2.data

        self.orthonormalized = True

    def subspace_diagonalize(self,
                             Ht,
                             dH,
                             work_array=None,
                             Htpsit=None,
                             scalapack_parameters=(None, 1, 1, -1)):
        """

        Ht(in, out)::

           ~   ^   ~
           H = T + v

        dH::

            ~  ~    a    ~  ~
          <psi|p> dH    <p|psi>
              m i   ij    j   n
        """
        self.orthonormalize(work_array)
        psit = self.wave_functions
        projections = self.projections
        psit2 = psit.new(data=work_array)
        projections2 = projections.new()
        domain_comm = psit.layout.comm

        Ht = partial(Ht, out=psit2, spin=0)
        dH = partial(dH, out=projections2, spin=0)
        H = psit.matrix_elements(psit, function=Ht, domain_sum=False)
        projections.matrix_elements(projections, function=dH,
                                    domain_sum=False, out=H, add_to_out=True)
        domain_comm.sum(H.data, 0)

        if domain_comm.rank == 0:
            slcomm, r, c, b = scalapack_parameters
            if r == c == 1:
                slcomm = None
            self._eigs = H.eigh(scalapack=(slcomm, r, c, b))
            # H.data[n, :] now contains the n'th eigenvector and eps_n[n]
            # the n'th eigenvalue
        domain_comm.broadcast(H.data, 0)
        domain_comm.broadcast(self.eigs, 0)
        if Htpsit is not None:
            H.multiply(psit2, out=Htpsit)

        H.multiply(psit, out=psit2)
        psit.data[:] = psit2.data
        projections.matrix.multiply(H, opb='T', out=projections2)
        projections.data[:] = projections2.data

    def force_contribution(self, dv: AtomArrays, F_av: Array2D):
        F_ainv = self.projectors.derivative(self.wave_functions)
        myoccs = self.weight * self.spin_degeneracy * self.myoccs
        for a, F_inv in F_ainv.items():
            F_inv = F_inv.conj()
            F_inv *= myoccs[:, np.newaxis]
            dH_ii = dv[a][:, :, self.spin]
            P_in = self.projections[a]
            F_vii = np.einsum('inv, jn, jk -> vik', F_inv, P_in, dH_ii)
            F_inv *= self.myeigs[:, np.newaxis]
            dO_ii = self.setups[a].dO_ii
            F_vii -= np.einsum('inv, jn, jk -> vik', F_inv, P_in, dO_ii)
            F_av[a] += 2 * F_vii.real.trace(0, 1, 2)
