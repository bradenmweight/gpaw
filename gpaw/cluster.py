"""Extensions to the ase Atoms class

"""
import numpy as np

from ase import Atoms
from ase.io import read
from ase.build.connected import connected_indices

from gpaw.core import UGDesc


class Cluster(Atoms):
    """A class for cluster structures
    to enable simplified manipulation"""

    def __init__(self, *args, **kwargs):

        self.data = {}

        if len(args) > 0:
            filename = args[0]
            if isinstance(filename, str):
                self.read(filename, kwargs.get('filetype'))
                return
        else:
            Atoms.__init__(self, [])

        if kwargs.get('filename') is not None:
            filename = kwargs.pop('filename')
            Atoms.__init__(self, *args, **kwargs)
            self.read(filename, kwargs.get('filetype'))
        else:
            Atoms.__init__(self, *args, **kwargs)

    def extreme_positions(self):
        """get the extreme positions of the structure"""
        pos = self.get_positions()
        return np.array([np.minimum.reduce(pos), np.maximum.reduce(pos)])

    def find_connected(self, index, dmax=None, scale=1.5):
        """Find atoms connected to self[index] and return them."""
        return self[connected_indices(self, index, dmax, scale)]

    def minimal_box(self, border=0, h=None, multiple=4):
        """The box needed to fit the structure in.

        The structure is moved to fit into the box [(0,x),(0,y),(0,z)]
        with x,y,z > 0 (fitting the ASE constriction).
        The border argument can be used to add a border of empty space
        around the structure.

        If h is set, the box is extended to ensure that box/h is
        a multiple of 'multiple'.
        This ensures that GPAW uses the desired h.

        The shift applied to the structure is returned.
         """

        if len(self) == 0:
            return None

        extr = self.extreme_positions()

        # add borders
        if isinstance(border, list):
            b = border
        else:
            b = [border, border, border]
        for c in range(3):
            extr[0][c] -= b[c]
            extr[1][c] += b[c] - extr[0][c]  # shifted already

        pbc = self.pbc
        old_cell = self.cell

        if True in pbc:

            extr2 = np.zeros((3, 3))

            for ip, p in enumerate(pbc):

                if p:
                    extr[0][ip] = 0
                    extr2[ip][:] = old_cell[ip]

                else:
                    e = np.zeros(3)
                    e[ip] = extr[1][ip]
                    extr2[ip][:] = e

        # check for multiple of 4
        if h is not None:

            if not hasattr(h, '__len__'):
                h = np.array([h, h, h])

                if True in pbc:
                    grid = UGDesc.from_cell_and_grid_spacing(extr2, h, pbc)
                    H = extr2 / grid.size
                    h = np.zeros(3)
                    h1 = 0
                    i = 0
                    for ip, p in enumerate(pbc):
                        if p:
                            h1 += np.linalg.norm(H[ip])
                            i += 1
                    h = [h1 / i, h1 / i, h1 / i]

            for c in range(3):

                if True in pbc:
                    if not pbc[c]:
                        L = np.linalg.norm(extr2[c])

                        N = np.ceil(L / h[c] / multiple) * multiple

                        # correct L
                        dL = N * h[c] - L
                        extr2[c, c] += dL
                        extr[0][c] -= dL / 2

                else:
                    # apply the same as in paw.py
                    L = extr[1][c]  # shifted already
                    N = np.ceil(L / h[c] / multiple) * multiple
                    # correct L
                    dL = N * h[c] - L
                    # move accordingly
                    extr[1][c] += dL  # shifted already
                    extr[0][c] -= dL / 2.

        # move lower corner to (0, 0, 0)
        shift = tuple(-1. * np.array(extr[0]))
        self.translate(shift)

        if True in pbc:
            self.set_cell(tuple(extr2))
        else:

            self.set_cell(tuple(extr[1]))
        return shift

    def read(self, filename, format=None):
        """Read the structure from some file. The type can be given
        or it will be guessed from the filename."""

        self.__init__(read(filename, format=format))
        return len(self)
