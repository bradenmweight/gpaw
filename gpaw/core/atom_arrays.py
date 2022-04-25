from __future__ import annotations

import numbers
from typing import Sequence

import numpy as np
from gpaw.mpi import MPIComm, serial_comm
from gpaw.typing import Array1D, ArrayLike1D
from gpaw.core.matrix import Matrix


class AtomArraysLayout:
    def __init__(self,
                 shapes: list[int | tuple[int, ...]],
                 atomdist: AtomDistribution | MPIComm = serial_comm,
                 dtype=float):
        """Description of layour of atom arrays.

        Parameters
        ----------
        shapes:
            Shapse of arrays - one for each atom.
        atomdist:
            Distribution of atoms.
        dtype:
            Data-type (float or complex).
        """
        self.shape_a = [shape if isinstance(shape, tuple) else (shape,)
                        for shape in shapes]
        if not isinstance(atomdist, AtomDistribution):
            atomdist = AtomDistribution(np.zeros(len(shapes), int), atomdist)
        self.atomdist = atomdist
        self.dtype = np.dtype(dtype)

        self.size = sum(np.prod(shape) for shape in self.shape_a)

        self.myindices = []
        self.mysize = 0
        I1 = 0
        for a in atomdist.indices:
            I2 = I1 + np.prod(self.shape_a[a])
            self.myindices.append((a, I1, I2))
            self.mysize += I2 - I1
            I1 = I2

    def __repr__(self):
        return (f'AtomArraysLayout({self.shape_a}, {self.atomdist}, '
                f'{self.dtype})')

    def new(self, atomdist=None):
        """Create new AtomsArrayLayout object with new atomdist."""
        return AtomArraysLayout(self.shape_a, atomdist or self.atomdist,
                                self.dtype)

    def empty(self,
              dims: int | tuple[int, ...] = (),
              comm: MPIComm = serial_comm) -> AtomArrays:
        """Create new AtomArrays object.

        parameters
        ----------
        dims:
            Extra dimensions.
        comm:
            Distribute dimensions along this communicator.
        """
        return AtomArrays(self, dims, comm)

    def sizes(self) -> tuple[list[dict[int, int]], Array1D]:
        """Compute array sizes for all ranks.

        >>> AtomArraysLayout([3, 4]).sizes()
        ([{0: 3, 1: 4}], array([7]))
        """
        comm = self.atomdist.comm
        size_ra: list[dict[int, int]] = [{} for _ in range(comm.size)]
        size_r = np.zeros(comm.size, int)
        for a, (rank, shape) in enumerate(zip(self.atomdist.rank_a,
                                              self.shape_a)):
            size = np.prod(shape)
            size_ra[rank][a] = size
            size_r[rank] += size
        return size_ra, size_r


class AtomDistribution:
    def __init__(self, ranks: ArrayLike1D, comm: MPIComm = serial_comm):
        """Atom-distribution.

        Parameters
        ----------
        ranks:
            List of ranks, one rank per atom.
        comm:
            MPI-communicator.
        """
        self.comm = comm
        self.rank_a = ranks
        self.indices = np.where(ranks == comm.rank)[0]

    @classmethod
    def from_number_of_atoms(cls,
                             natoms: int,
                             comm: MPIComm = serial_comm) -> AtomDistribution:
        """Distribute atoms evenly.

        >>> AtomDistribution.from_number_of_atoms(3).rank_a
        array([0, 0, 0])
        """
        blocksize = (natoms + comm.size - 1) // comm.size
        rank_a = np.empty(natoms, int)
        a1 = 0
        for rank in range(comm.size):
            a2 = a1 + blocksize
            rank_a[a1:a2] = rank
            if a2 >= natoms:
                break
            a1 = a2
        return cls(rank_a, comm)

    @classmethod
    def from_atom_indices(cls,
                          atom_indices: Sequence[int],
                          comm: MPIComm = serial_comm,
                          *,
                          natoms: int = None) -> AtomDistribution:
        """Create distribution from atom indices.

        >>> AtomDistribution.from_atom_indices([0, 1, 2]).rank_a
        array([0, 0, 0])
        """
        if natoms is None:
            natoms = comm.max(max(atom_indices)) + 1
        rank_a = np.zeros(natoms, int)  # type: ignore
        rank_a[atom_indices] = comm.rank
        comm.sum(rank_a)
        return cls(rank_a, comm)

    def __repr__(self):
        return (f'AtomDistribution(ranks={self.rank_a}, '
                f'comm={self.comm.rank}/{self.comm.size})')

    def gather(self):
        return AtomDistribution(np.zeros(len(self.rank_a), int))


class AtomArrays:
    def __init__(self,
                 layout: AtomArraysLayout,
                 dims: int | tuple[int, ...] = (),
                 comm: MPIComm = serial_comm,
                 data: np.ndarray = None):
        """AtomArrays object.

        parameters
        ----------
        layout:
            Layout-description.
        dims:
            Extra dimensions.
        comm:
            Distribute dimensions along this communicator.
        data:
            Data array for storage.
        """
        myshape = (layout.mysize,)
        domain_comm = layout.atomdist.comm
        dtype = layout.dtype

        self.myshape = myshape
        self.comm = comm
        self.domain_comm = domain_comm

        # convert int to tuple:
        self.dims = dims if isinstance(dims, tuple) else (dims,)

        if self.dims:
            mydims0 = (self.dims[0] + comm.size - 1) // comm.size
            d1 = min(comm.rank * mydims0, self.dims[0])
            d2 = min((comm.rank + 1) * mydims0, self.dims[0])
            mydims0 = d2 - d1
            self.mydims = (mydims0,) + self.dims[1:]
        else:
            self.mydims = ()

        fullshape = self.mydims + self.myshape

        if data is not None:
            if data.shape != fullshape:
                raise ValueError(
                    f'Bad shape for data: {data.shape} != {fullshape}')
            if data.dtype != dtype:
                raise ValueError(
                    f'Bad dtype for data: {data.dtype} != {dtype}')
        else:
            data = np.empty(fullshape, dtype)

        self.data = data
        self._matrix: Matrix | None = None

        self.layout = layout
        self._arrays = {}
        for a, I1, I2 in layout.myindices:
            self._arrays[a] = self.data[..., I1:I2].reshape(
                self.mydims + layout.shape_a[a])
        self.natoms: int = len(layout.shape_a)

    def __repr__(self):
        txt = f'AtomArrays({self.layout}, dims={self.dims}'
        if self.comm.size > 1:
            txt += f', comm={self.comm.rank}/{self.comm.size}'
        return txt + ')'

    @property
    def matrix(self) -> Matrix:
        if self._matrix is not None:
            return self._matrix

        shape = (self.dims[0], np.prod(self.dims[1:]) * np.prod(self.myshape))
        myshape = (self.mydims[0], np.prod(self.mydims[1:]) * np.prod(self.myshape))
        dist = (self.comm, -1, 1)

        data = self.data.reshape(myshape)
        self._matrix = Matrix(*shape, data=data, dist=dist)

        return self._matrix

    def new(self, layout=None, data=None):
        """Create new AtomArrays object of same kind.

        Parameters
        ----------
        layout:
            Layout-description.
        data:
            Array to use for storage.
        """
        return AtomArrays(layout or self.layout,
                          self.dims,
                          self.comm,
                          data=data)

    def __getitem__(self, a):
        if isinstance(a, numbers.Integral):
            return self._arrays[a]
        if len(self.dims) == 1:
            a_ai = AtomArrays(self.layout, data=self.data[a[1]].copy())
            return a_ai
        1 / 0

    def get(self, a):
        return self._arrays.get(a)

    def __setitem__(self, a, value):
        self._arrays[a][:] = value

    def __contains__(self, a):
        return a in self._arrays

    def items(self):
        return self._arrays.items()

    def keys(self):
        return self._arrays.keys()

    def values(self):
        return self._arrays.values()

    def gather(self, broadcast=False, copy=False) -> AtomArrays | None:
        """Gather all atoms on master."""
        comm = self.layout.atomdist.comm
        if comm.size == 1:
            if copy:
                aa = self.new()
                aa.data[:] = self.data
                return aa
            return self

        if comm.rank == 0 or broadcast:
            aa = self.new(layout=self.layout.new(atomdist=serial_comm))
        else:
            aa = None

        if comm.rank == 0:
            size_ra, size_r = self.layout.sizes()
            shape = self.mydims + (size_r.max(),)
            buffer = np.empty(shape, self.layout.dtype)
            for rank in range(1, comm.size):
                buf = buffer[..., :size_r[rank]]
                comm.receive(buf, rank)
                b1 = 0
                for a, size in size_ra[rank].items():
                    b2 = b1 + size
                    A = aa[a]
                    A[:] = buf[..., b1:b2].reshape(A.shape)
            for a, array in self._arrays.items():
                aa[a] = array
        else:
            comm.send(self.data, 0)

        if broadcast:
            comm.broadcast(aa.data, 0)

        return aa

    def scatter_from(self, data: np.ndarray = None) -> None:
        comm = self.layout.atomdist.comm
        if comm.size == 1:
            self.data[:] = data
            return

        if comm.rank != 0:
            comm.receive(self.data, 0, 42)
            return

        size_ra, size_r = self.layout.sizes()
        aa = self.new(layout=self.layout.new(atomdist=serial_comm),
                      data=data)
        requests = []
        for rank, (totsize, size_a) in enumerate(zip(size_r, size_ra)):
            if rank != 0:
                buf = np.empty(self.mydims + (totsize,), self.layout.dtype)
                b1 = 0
                for a, size in size_a.items():
                    b2 = b1 + size
                    buf[..., b1:b2] = aa[a].reshape(self.mydims + (size,))
                request = comm.send(buf, rank, 42, False)
                # Remember to store a reference to the
                # send buffer (buf) so that is isn't
                # deallocated
                requests.append((request, buf))
            else:
                for a in size_a:
                    self[a] = aa[a]

        for request, _ in requests:
            comm.wait(request)

    def to_lower_triangle(self):
        """Convert `N*N` matrices to `N*(N+1)/2` vectors.

        >>> a = AtomArraysLayout([(3, 3)]).empty()
        >>> a[0][:] = [[11, 12, 13],
        ...            [12, 22, 23],
        ...            [13, 23, 33]]
        >>> a.to_lower_triangle()[0]
        array([11., 12., 22., 13., 23., 33.])
        """
        shape_a = []
        for i1, i2 in self.layout.shape_a:
            assert i1 == i2
            shape_a.append((i1 * (i1 + 1) // 2,))
        layout = AtomArraysLayout(shape_a, self.layout.atomdist.comm)
        a_axp = layout.empty(self.dims)
        for a_xii, a_xp in zip(self.values(), a_axp.values()):
            i = a_xii.shape[-1]
            L = np.tril_indices(i)
            for a_p, a_ii in zip(a_xp.reshape((-1, i * (i + 1) // 2)),
                                 a_xii.reshape((-1, i, i))):
                a_p[:] = a_ii[L]
        return a_axp

    def to_full(self):
        r"""Convert `N(N+1)/2` vectors to `N\times N` matrices.

        >>> a = AtomArraysLayout([6]).empty()
        >>> a[0][:] = [1, 2, 3, 4, 5, 6]
        >>> a.to_full()[0]
        array([[1., 2., 4.],
               [2., 3., 5.],
               [4., 5., 6.]])
        """
        shape_a = []
        for (p,) in self.layout.shape_a:
            i = int((2 * p + 0.25)**0.5)
            shape_a.append((i, i))
        layout = AtomArraysLayout(shape_a, self.layout.atomdist.comm)
        a_axii = layout.empty(self.dims)
        for a_xp, a_xii in zip(self.values(), a_axii.values()):
            i = a_xii.shape[-1]
            a_xii[(...,) + np.tril_indices(i)] = a_xp
            u = (...,) + np.triu_indices(i, 1)
            a_xii[u] = np.swapaxes(a_xii, -1, -2)[u]
        return a_axii
