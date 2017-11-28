import numpy as np

from gpaw.tddft.units import au_to_eV
from gpaw.tddft.units import eV_to_au


def frequencies(frequencies, foldings, widths, units='eV'):
    f_w = []
    for freq in np.array([frequencies]).ravel():
        for folding in np.array([foldings]).ravel():
            for width in np.array([widths]).ravel():
                f_w.append(Frequency(freq, folding, width, units))
    return f_w


class Frequency(object):
    def __init__(self, frequency, folding, width, units='eV'):
        self.frequency = float(frequency)

        if width is None:
            folding = None

        self.folding = folding
        if self.folding is None:
            self.width = None
        else:
            self.width = float(width)

        if units == 'eV':
            for arg in ['frequency', 'width']:
                if getattr(self, arg) is not None:
                    setattr(self, arg, getattr(self, arg) * eV_to_au)
        elif units != 'au':
            raise RuntimeError('Unknown units: %s' % units)

        if self.folding not in [None, 'Gauss', 'Lorentz']:
            raise RuntimeError('Unknown folding: %s' % self.folding)

    def envelope(self, time):
        if self.folding is None:
            return 1.0
        elif self.folding == 'Gauss':
            return np.exp(- 0.5 * self.width**2 * time**2)
        elif self.folding == 'Lorentz':
            return np.exp(- self.width * time)

    def todict(self):
        d = dict(units='au')
        for arg in ['frequency', 'folding', 'width']:
            d[arg] = getattr(self, arg)
        return d

    def __repr__(self):
        s = '%.5f eV' % (self.frequency * au_to_eV)
        if self.folding is None:
            return s
        s += ' w %s(%.5f eV)' % (self.folding, self.width * au_to_eV)
        return s
