from __future__ import annotations

from collections.abc import Sequence

import numpy as np

from gpaw.typing import Array1D
from gpaw.response import timer, ResponseContext
from gpaw.response.pair_functions import Chi
from gpaw.response.fxc_kernels import FXCKernel
from gpaw.response.goldstone import get_goldstone_scaling


class HXCScaling:
    """Helper for scaling hxc kernels."""

    def __init__(self, mode, lambd=None):
        self.mode = mode
        self._lambd = lambd

    @property
    def lambd(self):
        return self._lambd

    def calculate_scaling(self, chiks, Khxc_GG, dyson_solver):
        if chiks.spincomponent in ['+-', '-+']:
            self._lambd = get_goldstone_scaling(
                self.mode, chiks, Khxc_GG, dyson_solver)
        else:
            raise ValueError('No scaling method implemented for '
                             f'spincomponent={chiks.spincomponent}')


class HXCKernel:
    """Hartree-exchange-correlation kernel in a plane-wave basis."""

    def __init__(self,
                 Vbare_G: Array1D | None = None,
                 fxc_kernel: FXCKernel | None = None):
        """Construct the Hxc kernel."""
        self.Vbare_G = Vbare_G
        self.fxc_kernel = fxc_kernel

        if Vbare_G is None:
            assert fxc_kernel is not None
            self.nG = fxc_kernel.GG_shape[0]
        else:
            self.nG = len(Vbare_G)
            if fxc_kernel is not None:
                assert fxc_kernel.GG_shape[0] == self.nG

    def get_Khxc_GG(self):
        """Hartree-exchange-correlation kernel."""
        # Allocate array
        Khxc_GG = np.zeros((self.nG, self.nG), dtype=complex)
        if self.Vbare_G is not None:  # Add the Hartree kernel
            Khxc_GG.flat[::self.nG + 1] += self.Vbare_G
        if self.fxc_kernel is not None:  # Add the xc kernel
            # Unfold the fxc kernel into the Kxc kernel matrix
            Khxc_GG += self.fxc_kernel.get_Kxc_GG()
        return Khxc_GG


class DysonSolver:
    """Class for invertion of Dyson-like equations."""

    def __init__(self, context):
        self.context = context

    @timer('Solve Dyson-like equations')
    def __call__(self, chiks: Chi, hxc_kernel: HXCKernel,
                 hxc_scaling: HXCScaling | None = None) -> Chi:
        """Solve the dyson equation and return the many-body susceptibility."""
        dyson_equations = DysonEquations(chiks, hxc_kernel)
        return dyson_equations.invert(
            hxc_scaling=hxc_scaling, context=self.context)


class DysonEquations(Sequence):
    """Sequence of Dyson-like equations at different complex frequencies z."""

    def __init__(self, chiks: Chi, hxc_kernel: HXCKernel):
        assert chiks.distribution == 'zGG' and\
            chiks.blockdist.fully_block_distributed, \
            "DysonSolver needs chiks' frequencies to be distributed over world"
        nG = hxc_kernel.nG
        assert chiks.array.shape[1:] == (nG, nG)
        self.zblocks = chiks.blocks1d
        self.chiks = chiks
        self.Khxc_GG = hxc_kernel.get_Khxc_GG()

    def __len__(self):
        return self.zblocks.nlocal

    def __getitem__(self, z):
        return DysonEquation(self.chiks.array[z], self.Khxc_GG)

    def invert(self,
               hxc_scaling: HXCScaling | None = None,
               context: ResponseContext | None = None) -> Chi:
        """Invert Dyson equations to obtain χ(z)."""
        txtout = []
        if hxc_scaling is None:
            lambd = None  # no scaling of the self-enhancement function
        else:
            if hxc_scaling.lambd is None:  # calculate, if not already
                hxc_scaling.calculate_scaling(self)
            lambd = hxc_scaling.lambd
            txtout.append(r'Rescaling the self-enhancement function by a '
                          f'factor of  λ={lambd}')

        if context is not None:
            txtout.append('Inverting Dyson-like equation')
            context.print('\n'.join(txtout))

        chi = self.chiks.new()
        for z, dyson_equation in enumerate(self):
            chi.array[z] = dyson_equation.invert(lambd=lambd)

        return chi


class DysonEquation:
    """Dyson equation at wave vector q and frequency z.

    The Dyson equation is given in plane-wave components as

    χ(q,z) = χ_KS(q,z) + Ξ(q,z) χ(q,z),

    where the self-enhancement function encodes the electron correlations
    induced by by the effective (Hartree-exchange-correlation) interaction:

    Ξ(q,z) = χ_KS(q,z) K_hxc(q,z)

    See [to be published] for more information.
    """

    def __init__(self, chiks_GG, Khxc_GG):
        self.nG = chiks_GG.shape[0]
        self.chiks_GG = chiks_GG
        self.Khxc_GG = Khxc_GG

    @property
    def xi_GG(self):
        """Calculate the self-enhancement function."""
        return self.chiks_GG @ self.Khxc_GG

    def invert(self, lambd: float | None = None):
        """Invert the Dyson equation (with or without a rescaling of Ξ).

        χ(q,z) = [1 - λ Ξ(q,z)]^(-1) χ_KS(q,z)
        """
        if lambd is None:
            lambd = 1.  # no rescaling
        enhancement_GG = np.linalg.inv(
            np.eye(self.nG) - lambd * self.xi_GG)
        return enhancement_GG @ self.chiks_GG


class DysonEnhancer:
    """Class for applying self-enhancement functions."""
    def __init__(self, context):
        self.context = context

    def __call__(self, chiks: Chi, xi: Chi) -> Chi:
        """Solve the Dyson equation and return the many-body susceptibility."""
        assert chiks.distribution == 'zGG' and \
            chiks.blockdist.fully_block_distributed
        assert xi.distribution == 'zGG' and \
            xi.blockdist.fully_block_distributed
        assert chiks.spincomponent == xi.spincomponent
        assert np.allclose(chiks.zd.hz_z, xi.zd.hz_z)
        assert np.allclose(chiks.qpd.q_c, xi.qpd.q_c)

        chi = chiks.new()
        chi.array = self.invert_dyson(chiks.array, xi.array)

        return chi

    @timer('Invert Dyson-like equation')
    def invert_dyson(self, chiks_zGG, xi_zGG):
        r"""Invert the frequency dependent Dyson equation in plane-wave basis:
                                           __
                                           \
        χ_GG'^+-(q,z) = χ_KS,GG'^+-(q,z) + /  Ξ_GG1^++(q,z) χ_G1G'^+-(q,z)
                                           ‾‾
                                           G1
        """
        self.context.print('Inverting Dyson-like equation')
        chi_zGG = np.empty_like(chiks_zGG)
        for chi_GG, chiks_GG, xi_GG in zip(chi_zGG, chiks_zGG, xi_zGG):
            chi_GG[:] = self.invert_dyson_single_frequency(chiks_GG, xi_GG)
        return chi_zGG

    @staticmethod
    def invert_dyson_single_frequency(chiks_GG, xi_GG):
        enhancement_GG = np.linalg.inv(np.eye(len(chiks_GG)) - xi_GG)
        chi_GG = enhancement_GG @ chiks_GG
        return chi_GG
