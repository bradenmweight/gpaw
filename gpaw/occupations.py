# Copyright (C) 2003  CAMP
# Please see the accompanying LICENSE file for further information.

"""Occpation number objects."""

from ase.units import Hartree
import numpy as np

class OccupationNumbers:
    """Base class for all occupation number objects."""
    def __init__(self, fixmagmom):
        self.fixmagmom = fixmagmom        
        self.magmom = None     # magnetic moment
        self.e_entropy = None  # -ST
        self.e_band = None     # band energy (sum_n eps_n * f_n)
        self.fermilevel = None # Fermi level
        self.homo = np.nan     # HOMO eigenvalue
        self.lumo = np.nan     # LUMO eigenvalue
        self.nvalence = None   # number of electrons
        self.split = 0.0       # splitting of Fermi levels from fixmagmom=True
        self.niter = 0         # number of iterations for finding Fermi level
        
    def calculate(self, wfs):
        """Calculate everything.

        The following is calculated:

        * occupation numbers
        * magnetic moment
        * entropy
        * band energy
        * Fermi level
        * HOMO and LUMO energies
        """

        # Allocate:
        for kpt in wfs.kpt_u:
            if kpt.f_n is None:
                kpt.f_n = np.empty(wfs.mynbands)

        # Allow subclasses to adjust nvalence:
        self.set_number_of_electrons(wfs)

        # Let the master domaindo the work and broadcast results:
        data = np.empty(7)
        if wfs.gd.comm.rank == 0:
            self.calculate_occupation_numbers(wfs)
            self.calculate_band_energy(wfs)
            data[:] = [self.magmom, self.e_entropy, self.e_band,
                       self.homo, self.lumo,
                       self.fermilevel, self.split]
        wfs.world.broadcast(data, 0)
        (self.magmom, self.e_entropy,
         self.e_band, self.homo, self.lumo, self.fermilevel, self.split) = data

        for kpt in wfs.kpt_u:
            wfs.gd.comm.broadcast(kpt.f_n, 0)

    def set_number_of_electrons(self, wfs):
        self.nvalence = wfs.nvalence

    def calculate_occupation_numbers(self, wfs):
        raise NotImplementedError

    def calculate_band_energy(self, wfs):
        """Sum up all eigenvalues weighted with occupation numbers"""
        e_band = 0.0
        for kpt in wfs.kpt_u:
            e_band += np.dot(kpt.f_n, kpt.eps_n)    
        self.e_band = wfs.bd.comm.sum(wfs.kpt_comm.sum(e_band))

    def print_fermi_level(self, stream):
        pass

    def get_fermi_level(self):
        raise ValueError('Can not calculate Fermi level!')

    def set_fermi_level(self, fermilevel):
        self.fermilevel = fermilevel

def occupy(f_n, eps_n, ne, weight=1):
    """Fill in occupation numbers.

    return HOMO and LUMO energies."""

    N = len(f_n)
    if ne == N * weight:
        f_n[:] = weight
        return eps_n[-1], np.inf

    n, f = divmod(ne, weight)
    n = int(n)
    f_n[:n] = weight
    assert n < N
    f_n[n] = f
    f_n[n + 1:] = 0.0
    if f > 0.0:
        return eps_n[n], eps_n[n]
    return eps_n[n - 1], eps_n[n]

class ZeroKelvin(OccupationNumbers):
    def __init__(self, fixmagmom):
        OccupationNumbers.__init__(self, fixmagmom)
        
    def calculate_occupation_numbers(self, wfs):
        assert wfs.gamma
        if self.fixmagmom:
            self.fixed_moment(wfs)
        elif wfs.nspins == 1:
            self.spin_paired(wfs)
        else:
            self.spin_polarized(wfs)
        self.e_entropy = 0.0

    def print_fermi_level(self, stream):
        if self.fermilevel is not None and np.isfinite(self.fermilevel):
            if self.split == 0.0:
                stream.write('Fermi Level: %.5f\n' %
                             (Hartree * self.fermilevel))
            else:
                stream.write('Fermi Levels: %.5f, %.5f\n' %
                             (Hartree * (self.fermilevel + 0.5 * self.split),
                              Hartree * (self.fermilevel - 0.5 * self.split)))

    def get_fermi_level(self):
        if self.fermilevel is None or not np.isfinite(self.fermilevel):
            OccupationNumbers.get_fermi_level(self)  # fail
        else:
            return self.fermilevel

    def get_homo_lumo(self, wfs):
        if self.nvalence is None:
            self.calculate(wfs)
        if np.isfinite(self.homo) and np.isfinite(self.lumo):
            return np.array([self.homo, self.lumo])
        else:
            raise ValueError("Can't find HOMO and/or LUMO!")

    def fixed_moment(self, wfs):
        assert wfs.nspins == 2 and wfs.bd.comm.size == 1
        fermilevels = np.zeros(2)
        for kpt in wfs.kpt_u:
            eps_n = wfs.bd.collect(kpt.eps_n)
            f_n = np.empty(wfs.nbands)
            sign = 1 - kpt.s * 2
            ne = 0.5 * (self.nvalence + sign * self.magmom)
            homo, lumo = occupy(f_n, eps_n, ne)
            wfs.bd.distribute(f_n, kpt.f_n)
            fermilevels[kpt.s] = 0.5 * (homo + lumo)
        wfs.kpt_comm.sum(fermilevels)
        self.fermilevel = fermilevels.mean()
        self.split = fermilevels[0] - fermilevels[1]
        
    def spin_paired(self, wfs):
        kpt = wfs.kpt_u[0]
        eps_n = wfs.bd.collect(kpt.eps_n)
        if wfs.bd.comm.rank == 0:
            f_n = np.empty(wfs.nbands)
            self.homo, self.lumo = occupy(f_n, eps_n, self.nvalence, 2)
            self.fermilevel = 0.5 * (self.homo + self.lumo)
        else:
            f_n = None
            self.fermilevel = np.nan
        wfs.bd.distribute(f_n, kpt.f_n)
        self.magmom = 0.0
        
    def spin_polarized(self, wfs):
        eps_un = [wfs.bd.collect(kpt.eps_n) for kpt in wfs.kpt_u]
        self.fermilevel = np.nan
        if wfs.bd.comm.rank == 0:
            if wfs.kpt_comm.size == 2:
                if wfs.kpt_comm.rank == 1:
                    wfs.kpt_comm.send(eps_un[0], 0)
                else:
                    eps_sn = [eps_un[0], np.empty(wfs.nbands)]
                    wfs.kpt_comm.receive(eps_sn[1], 1)
            else:
                eps_sn = eps_un

            if wfs.kpt_comm.rank == 0:
                eps_n = np.ravel(eps_sn)
                f_n = np.empty(wfs.nbands * 2)
                nsorted = eps_n.argsort()
                self.homo, self.lumo = occupy(f_n, eps_n[nsorted],
                                              self.nvalence)
                f_sn = f_n[nsorted.argsort()].reshape((2, wfs.nbands))
                self.magmom = f_sn[0].sum() - f_sn[1].sum()
                self.fermilevel = 0.5 * (self.homo + self.lumo)

            if wfs.kpt_comm.size == 2:
                if wfs.kpt_comm.rank == 0:
                    wfs.kpt_comm.send(f_sn[1], 1)
                else:
                    f_sn = [None, np.empty(wfs.nbands)]
                    wfs.kpt_comm.receive(f_sn[1], 0)
        else:
            f_sn = [None, None]

        for kpt in wfs.kpt_u:
            wfs.bd.distribute(f_sn[kpt.s], kpt.f_n)

class SmoothDistribution(ZeroKelvin):
    """Base class for Fermi-Dirac and other smooth distributions."""
    def __init__(self, width, fixmagmom, maxiter):
        """Smooth distribution.

        Find the Fermi level by integrating in energy until
        the number of electrons is correct.

        width: float
            Width of distribution in eV.
        fixmagmom: bool
            Fix spin moment calculations.  A separate Fermi level for
            spin up and down electrons is found: self.fermilevel +
            self.split and self.fermilevel - self.split.
        maxiter: int
            Maximum number of iterations.
        """

        ZeroKelvin.__init__(self, fixmagmom)
        self.width = width / Hartree
        self.maxiter = maxiter
        
    def calculate_occupation_numbers(self, wfs):
        if self.width == 0 or self.nvalence == wfs.nbands * 2:
            ZeroKelvin.calculate_occupation_numbers(self, wfs)
            return

        if self.fermilevel is None:
            self.fermilevel = self.guess_fermi_level(wfs)

        if not self.fixmagmom:
            self.fermilevel, self.magmom, self.e_entropy = \
                             self.find_fermi_level(wfs, self.nvalence,
                                                   self.fermilevel)
            if wfs.nspins == 1:
                self.magmom = 0.0
        else:
            fermilevels = np.empty(2)
            self.e_entropy = 0.0
            for s in range(2):
                sign = 1 - s * 2
                ne = 0.5 * (self.nvalence + sign * self.magmom)
                fermilevel = self.fermilevel + 0.5 * sign * self.split
                fermilevels[s], magmom, e_entropy = \
                                self.find_fermi_level(wfs, ne, fermilevel, [s])
                self.e_entropy += e_entropy
            self.fermilevel = fermilevels.mean()
            self.split = fermilevels[0] - fermilevels[1]

    def get_homo_lumo(self, wfs):
        if self.width == 0:
            return ZeroKelvin.get_homo_lumo(self, wfs)
        
        if wfs.nspins == 2:
            raise NotImplementedError
        
        n = self.nvalence // 2
        homo = wfs.world.max(max([kpt.eps_n[n - 1] for kpt in wfs.kpt_u]))
        lumo = -wfs.world.max(-min([kpt.eps_n[n] for kpt in wfs.kpt_u]))
        return np.array([homo, lumo])
        
    def guess_fermi_level(self, wfs):
        fermilevel = 0.0
        myeps_n = np.ravel([wfs.bd.collect(kpt.eps_n) for kpt in wfs.kpt_u])
        if wfs.bd.comm.rank == 0:
            if wfs.kpt_comm.rank > 0:
                wfs.kpt_comm.gather(myeps_n, 0)
            else:
                eps_n = np.empty(wfs.nspins * wfs.nibzkpts * wfs.nbands)
                wfs.kpt_comm.gather(myeps_n, 0, eps_n)
                eps_n = eps_n.ravel()
                eps_n.sort()
                n, f = divmod(self.nvalence * wfs.nibzkpts, 3 - wfs.nspins)
                n = int(n)
                if f > 0.0:
                    fermilevel = eps_n[n]
                else:
                    fermilevel = 0.5 * (eps_n[n - 1] + eps_n[n])

        # XXX broadcast would be better!
        return wfs.bd.comm.sum(wfs.kpt_comm.sum(fermilevel))
                    
    def find_fermi_level(self, wfs, ne, fermilevel, spins=(0, 1)):
        niter = 0
        while True:
            data = np.zeros(4)
            for kpt in wfs.kpt_u:
                if kpt.s in spins:
                    data += self.distribution(kpt, fermilevel)
            wfs.kpt_comm.sum(data)
            wfs.bd.comm.sum(data)
            n, dnde, magmom, e_entropy = data
            dn = ne - n
            if abs(dn) < 1e-9:
                break
            if abs(dnde) < 1e-9:
                fermilevel = self.guess_fermi_level(wfs)
                niter += 1
                if niter > self.maxiter:
                    raise RuntimeError('Could not locate the Fermi level! ' +
                                       'See ticket #27.')
                continue
            if niter > self.maxiter:
                raise RuntimeError('Could not locate the Fermi level!')
            de = dn / dnde
            if abs(de) > self.width:
                de *= self.width / abs(de)
            fermilevel += de
            niter += 1

        self.niter = niter
        return fermilevel, magmom, e_entropy

class FermiDirac(SmoothDistribution):
    def __init__(self, width, fixmagmom=False, maxiter=1000):
        SmoothDistribution.__init__(self, width, fixmagmom, maxiter)

    def distribution(self, kpt, fermilevel):
        x = (kpt.eps_n - fermilevel) / self.width
        x = x.clip(-100, 100)
        y = np.exp(x)
        z = y + 1.0
        kpt.f_n[:] = kpt.weight / z
        n = kpt.f_n.sum()
        dnde = (n - (kpt.f_n**2).sum() / kpt.weight) / self.width        
        y *= x
        y /= z
        y -= np.log(z)
        e_entropy = -kpt.weight * y.sum() * self.width
        sign = 1 - kpt.s * 2
        return np.array([n, dnde, n * sign, e_entropy])
