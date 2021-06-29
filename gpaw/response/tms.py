import numpy as np

import gpaw.mpi as mpi
from gpaw.response.susceptibility import FourComponentSusceptibilityTensor

FCST = FourComponentSusceptibilityTensor


class TransverseMagneticSusceptibility(FCST):
    """Class calculating the transverse magnetic susceptibility
    and related physical quantities."""

    def __init__(self, *args, **kwargs):
        assert kwargs['fxc'] == 'ALDA'

        # Enable scaling to fit to Goldstone theorem
        if 'fxckwargs' in kwargs and 'fxc_scaling' in kwargs['fxckwargs']:
            self.fxc_scaling = kwargs['fxckwargs']['fxc_scaling']
        else:
            self.fxc_scaling = None

        FCST.__init__(self, *args, **kwargs)

    def get_macroscopic_component(self, spincomponent, q_c, frequencies,
                                  filename=None, txt=None):
        """Calculates the spatially averaged (macroscopic) component of the
        transverse magnetic susceptibility and writes it to a file.
        
        Parameters
        ----------
        spincomponent : str
            '+-': calculate chi+-, '-+: calculate chi-+
        q_c, frequencies, filename, txt : see gpaw.response.susceptibility

        Returns
        -------
        see gpaw.response.susceptibility
        """
        assert spincomponent in ['+-', '-+']

        return FCST.get_macroscopic_component(self, spincomponent, q_c,
                                              frequencies, filename=filename,
                                              txt=txt)

    def get_component_array(self, spincomponent, q_c, frequencies,
                            array_ecut=50, filename=None, txt=None):
        """Calculates a specific spin component of the
        transverse magnetic susceptibility and writes it to a file.
        
        Parameters
        ----------
        spincomponent : str
            '+-': calculate chi+-, '-+: calculate chi-+
        q_c, frequencies,
        array_ecut, filename, txt : see gpaw.response.susceptibility

        Returns
        -------
        see gpaw.response.susceptibility
        """
        assert spincomponent in ['+-', '-+']

        return FCST.get_component_array(self, spincomponent, q_c,
                                        frequencies, array_ecut=array_ecut,
                                        filename=filename, txt=txt)

    def _calculate_component(self, spincomponent, pd, wd):
        """Calculate a transverse magnetic susceptibility element.

        Returns
        -------
        pd, wd, chiks_wGG, chi_wGG : see gpaw.response.susceptibility
        """
        chiks_wGG = self.calculate_ks_component(spincomponent, pd,
                                                wd, txt=self.cfd)
        Kxc_GG = self.get_xc_kernel(spincomponent, pd,
                                    chiks_wGG=chiks_wGG, txt=self.cfd)

        chi_wGG = self.invert_dyson(chiks_wGG, Kxc_GG)

        return pd, wd, chiks_wGG, chi_wGG

    def get_xc_kernel(self, spincomponent, pd, chiks_wGG=None, txt=None):
        """Get the exchange correlation kernel."""
        Kxc_GG = self.fxc(spincomponent, pd, txt=self.cfd)

        fxc_scaling = self.fxc_scaling

        if fxc_scaling is not None:
            assert isinstance(fxc_scaling[0], bool)
            if fxc_scaling[0]:
                if fxc_scaling[1] is None:
                    assert pd.kd.gamma
                    print('Finding rescaling of kernel to fulfill the '
                          'Goldstone theorem', file=self.fd)
                    mode = fxc_scaling[2]
                    assert mode in ['fm', 'afm']
                    fxc_scaling[1] = get_goldstone_scaling(mode,
                                                           self.chiks.omega_w,
                                                           chiks_wGG, Kxc_GG,
                                                           world=self.world)

                assert isinstance(fxc_scaling[1], float)
                Kxc_GG *= fxc_scaling[1]

        self.fxc_scaling = fxc_scaling

        return Kxc_GG


def get_goldstone_scaling(mode, omega_w, chi0_wGG, Kxc_GG, world=mpi.world):
    """Get kernel scaling parameter fulfilling the Goldstone theorem."""
    # Find the frequency to determine the scaling from
    wgs = find_goldstone_frequency(mode, omega_w)

    # Only one rank, rgs, has the given frequency and finds the rescaling
    nw = len(omega_w)
    mynw = (nw + world.size - 1) // world.size
    rgs, mywgs = wgs // mynw, wgs % mynw
    fxcsbuf = np.empty(1, dtype=float)
    if world.rank == rgs:
        chi0_GG = chi0_wGG[mywgs]
        fxcsbuf[:] = find_goldstone_scaling(mode, chi0_GG, Kxc_GG)

    # Broadcast found rescaling
    world.broadcast(fxcsbuf, rgs)
    fxcs = fxcsbuf[0]

    return fxcs


def find_goldstone_frequency(mode, omega_w):
    """Factory function for finding the appropriate frequency to determine
    the kernel scaling from according to different Goldstone criteria."""
    assert mode in ['fm', 'afm'],\
        f"Allowed Goldstone scaling modes are 'fm', 'afm'. Got: {mode}"

    if mode == 'fm':
        return find_fm_goldstone_frequency(omega_w)
    elif mode == 'afm':
        return find_afm_goldstone_frequency(omega_w)


def find_fm_goldstone_frequency(omega_w):
    """Find omega=0. as the fm Goldstone frequency."""
    wgs = np.abs(omega_w).argmin()
    if not np.allclose(omega_w[wgs], 0., atol=1.e-8):
        raise ValueError("Frequency grid needs to include"
                         + " omega=0. to allow 'fm' Goldstone scaling")

    return wgs


def find_afm_goldstone_frequency(omega_w):
    """Find the second smallest positive frequency
    as the afm Goldstone frequency."""
    # Set omega=0. and negative frequencies to np.inf
    omega1_w = np.where(omega_w < 1.e-8, np.inf, omega_w)
    # Sort for the two smallest positive frequencies
    omega2_w = np.partition(omega1_w, 1)
    # Find original index of second smallest positive frequency
    wgs = np.abs(omega_w - omega2_w[1]).argmin()

    return wgs


def find_goldstone_scaling(mode, chi0_GG, Kxc_GG):
    """Factory function for finding the scaling of the kernel
    according to different Goldstone criteria."""
    assert mode in ['fm', 'afm'],\
        f"Allowed Goldstone scaling modes are 'fm', 'afm'. Got: {mode}"

    if mode == 'fm':
        return find_fm_goldstone_scaling(chi0_GG, Kxc_GG)
    elif mode == 'afm':
        return find_afm_goldstone_scaling(chi0_GG, Kxc_GG)


def find_fm_goldstone_scaling(chi0_GG, Kxc_GG):
    """Find goldstone scaling of the kernel by ensuring that the
    macroscopic inverse enhancement function has a root in (q=0, omega=0)."""
    fxcs = 1.
    kappaM = calculate_macroscopic_kappa(chi0_GG, Kxc_GG * fxcs)
    # If kappaM > 0, increase scaling (recall: kappaM ~ 1 - Kxc Re{chi_0})
    scaling_incr = 0.1 * np.sign(kappaM)
    while abs(kappaM) > 1.e-7 and abs(scaling_incr) > 1.e-7:
        fxcs += scaling_incr
        if fxcs <= 0.0 or fxcs >= 10.:
            raise Exception('Found an invalid fxc_scaling of %.4f' % fxcs)

        kappaM = calculate_macroscopic_kappa(chi0_GG, Kxc_GG * fxcs)

        # If kappaM changes sign, change sign and refine increment
        if np.sign(kappaM) != np.sign(scaling_incr):
            scaling_incr *= -0.2

    return fxcs


def find_afm_goldstone_scaling(chi0_GG, Kxc_GG):
    """Find goldstone scaling of the kernel by ensuring that the
    macroscopic magnon spectrum vanishes at q=0. for finite frequencies."""
    fxcs = 1.
    _, chiM = calculate_macroscopic_kappa(chi0_GG, Kxc_GG * fxcs)
    # If chi > 0., increase the scaling. If chi < 0., decrease the scaling.
    scaling_incr = 0.1 * np.sign(chiM)
    while (chiM < 0. or chiM > 1.e-7) or abs(scaling_incr) > 1.e-7:
        fxcs += scaling_incr
        if fxcs <= 0. or fxcs >= 10.:
            raise Exception('Found an invalid fxc_scaling of %.4f' % fxcs)

        _, chiM = calculate_macroscopic_kappa(chi0_GG, Kxc_GG * fxcs)

        # If chi changes sign, change sign and refine increment
        if np.sign(chiM) != np.sign(scaling_incr):
            scaling_incr *= -0.2

    return fxcs


def calculate_macroscopic_kappa(chi0_GG, Kxc_GG, return_chi=False):
    """Invert dyson equation and calculate the inverse enhancement function."""
    chi_GG = np.dot(np.linalg.inv(np.eye(len(chi0_GG)) +
                                  np.dot(chi0_GG, Kxc_GG * fxcs)),
                    chi0_GG)
    kappaM = (chi0_GG[0, 0] / chi_GG[0, 0]).real

    if not return_chi:
        return kappaM
    else:
        return kappaM, chi_GG[0, 0]
