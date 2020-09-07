from typing import Tuple, Dict, Any, Sequence

import numpy as np

from gpaw import GPAW
from gpaw.projections import Projections
from gpaw.utilities.partition import AtomPartition

Array2D = Any
Array4D = Any


class WannierOverlaps:
    def __init__(self,
                 cell: Sequence[Sequence[float]],
                 monkhorst_pack_size: Sequence[int],
                 directions: Dict[Tuple[int, int, int], int],
                 overlaps: Array4D):

        self.monkhorst_pack_size = tuple(monkhorst_pack_size)
        self.cell = np.array(cell)
        self.directions = directions

        nkpts, ndirs, self.nbands, nbands = overlaps.shape
        assert nbands == self.nbands
        assert nkpts == np.prod(monkhorst_pack_size)
        assert ndirs == len(directions)

        self._overlaps = overlaps

    def overlap(self,
                bz_index: int,
                direction: Tuple[int, int, int]) -> Array2D:
        return self._overlaps[bz_index, self.directions[direction]]

    def write(self, filename):
        1 / 0


class BZRealSpaceWaveFunctions:
    def __init__(self, kd, gd, u_knR, P_k):
        self.kd = kd
        self.gd = gd
        self.u_knR = u_knR
        self.P_k = P_k

    def __getitem__(self, bz_index):
        return self.u_knR[bz_index], self.P_k[bz_index]

    @classmethod
    def from_calculation(cls,
                         calc: GPAW,
                         n1: int = 0,
                         n2: int = 0,
                         spin=0) -> 'BZRealSpaceWaveFunctions':
        wfs = calc.wfs
        kd = wfs.kd

        gd = wfs.gd.new_descriptor(comm=calc.world)

        nproj_a = wfs.kpt_qs[0][0].projections.nproj_a
        # All atoms on rank-0:
        atom_partition = AtomPartition(gd.comm, np.zeros_like(nproj_a))

        u_nR = gd.empty((n2 - n1), complex, global_array=True)

        u_knR = gd.empty((kd.nbzkpts, n2 - n1), complex)
        P_k = []
        for ibz_index in range(kd.nibzkpts):
            for n in range(n1, n2):
                u_nR[n - n1] = wfs.get_wave_function_array(n=n,
                                                           k=ibz_index,
                                                           s=spin,
                                                           periodic=True)
            gd.distribute(u_nR, u_knR[ibz_index])

            P_nI = wfs.collect_projections(ibz_index, spin)
            projections = Projections(
                nbands=n2 - n1,
                nproj_a=nproj_a,
                atom_partition=atom_partition,
                data=P_nI[n1:n2])
            P_k.append(projections)

        return BZRealSpaceWaveFunctions(kd, gd, u_knR, P_k)


def calculate_overlaps(calc: GPAW,
                       n1: int = 0,
                       n2: int = 0,
                       soc: bool = False,
                       spin: int = 0) -> WannierOverlaps:
    if n2 <= 0:
        n2 += calc.get_number_of_bands()

    bzwfs = BZRealSpaceWaveFunctions.from_calculation(calc, n1, n2, spin)

    kd = bzwfs.kd
    gd = bzwfs.gd

    directions = {(1, 0, 0): 0, (0, 1, 0): 1, (0, 0, 1): 2}
    Z_kdnn = np.empty((kd.nbzkpts, len(directions), n2 - n1, n2 - n1), complex)

    size = kd.N_c
    spos_ac = calc.spos_ac
    setups = calc.wfs.setups

    for bz_index1 in range(kd.nbzkpts):
        u1_nR, P1_ani = bzwfs[bz_index1]
        i1_c = np.unravel_index(bz_index1, size)
        for direction, d in directions.items():
            i2_c = np.array(i1_c) + direction
            bz_index2 = np.ravel_multi_index(i2_c, size, 'wrap')
            u2_nR, P2_ani = bzwfs[bz_index2]
            phase_c = (i2_c % size - i2_c) // size
            if phase_c.any():
                u2_nR = u2_nR * gd.plane_wave(phase_c)
            Z_kdnn[bz_index1, d] = gd.integrate(u1_nR, u2_nR,
                                                global_integral=False)
            for a, P1_ni in P1_ani.items():
                dO_ii = setups[a].dO_ii
                P2_ni = P2_ani[a]
                Z_nn = P1_ni.conj().dot(dO_ii).dot(P2_ni.T).astype(complex)
                if phase_c.any():
                    Z_nn *= np.exp(2j * np.pi * phase_c.dot(spos_ac[a]))
                Z_kdnn[bz_index1, d] += Z_nn

    gd.comm.sum(Z_kdnn)

    overlaps = WannierOverlaps(calc.atoms.cell,
                               kd.N_c,
                               directions,
                               Z_kdnn)
    return overlaps
