import pickle
import warnings
from math import pi
from pathlib import Path

import numpy as np

from ase.parallel import paropen
from ase.units import Ha

from gpaw import GPAW, debug
import gpaw.mpi as mpi
from gpaw.hybrids.eigenvalues import non_self_consistent_eigenvalues
from gpaw.pw.descriptor import (count_reciprocal_vectors, PWMapping)
from gpaw.utilities.progressbar import ProgressBar

from gpaw.response import ResponseGroundStateAdapter, ResponseContext
from gpaw.response.chi0 import Chi0Calculator
from gpaw.response.pair import KPointPairFactory, phase_shifted_fft_indices
from gpaw.response.pair_functions import SingleQPWDescriptor
from gpaw.response.pw_parallelization import Blocks1D
from gpaw.response.screened_interaction import initialize_w_calculator
from gpaw.response.coulomb_kernels import CoulombKernel
from gpaw.response import timer


from ase.utils.filecache import MultiFileJSONCache as FileCache
from contextlib import ExitStack
from ase.parallel import broadcast


def compare_dicts(dict1, dict2, rel_tol=1e-14, abs_tol=1e-14):
    """
    Compare each key-value pair within dictionaries that contain nested data
    structures of arbitrary depth. If a kvp contains floats, you may specify
    the tolerance (abs or rel) to which the floats are compared. Individual
    elements within lists are not compared to floating point precision.

    :params dict1: Dictionary containing kvp to compare with other dictionary.
    :params dict2: Second dictionary.
    :params rel_tol: Maximum difference for being considered "close",
    relative to the magnitude of the input values as defined by math.isclose().
    :params abs_tol: Maximum difference for being considered "close",
    regardless of the magnitude of the input values as defined by
    math.isclose().

    :returns: bool indicating kvp's don't match (False) or do match (True)
    """
    from math import isclose
    if dict1.keys() != dict2.keys():
        return False

    for key in dict1.keys():
        val1 = dict1[key]
        val2 = dict2[key]

        if isinstance(val1, dict) and isinstance(val2, dict):
            # recursive func call to ensure nested structures are also compared
            if not compare_dicts(val1, val2, rel_tol, abs_tol):
                return False
        elif isinstance(val1, float) and isinstance(val2, float):
            if not isclose(val1, val2, rel_tol=rel_tol, abs_tol=abs_tol):
                return False
        else:
            if val1 != val2:
                return False

    return True


class Sigma:
    def __init__(self, iq, q_c, fxc, esknshape, **inputs):
        """Inputs are used for cache invalidation, and are stored for each
           file.
        """
        self.iq = iq
        self.q_c = q_c
        self.fxc = fxc
        self._buf = np.zeros((2, *esknshape))
        # self-energies and derivatives:
        self.sigma_eskn, self.dsigma_eskn = self._buf

        self.inputs = inputs

    def sum(self, comm):
        comm.sum(self._buf)

    def __iadd__(self, other):
        self.validate_inputs(other.inputs)
        self._buf += other._buf
        return self

    def validate_inputs(self, inputs):
        equals = compare_dicts(inputs, self.inputs, rel_tol=1e-14,
                               abs_tol=1e-14)
        if not equals:
            raise RuntimeError('There exists a cache with mismatching input '
                               f'parameters: {inputs} != {self.inputs}.')

    @classmethod
    def fromdict(cls, dct):
        instance = cls(dct['iq'], dct['q_c'], dct['fxc'],
                       dct['sigma_eskn'].shape, **dct['inputs'])
        instance.sigma_eskn[:] = dct['sigma_eskn']
        instance.dsigma_eskn[:] = dct['dsigma_eskn']
        return instance

    def todict(self):
        return {'iq': self.iq,
                'q_c': self.q_c,
                'fxc': self.fxc,
                'sigma_eskn': self.sigma_eskn,
                'dsigma_eskn': self.dsigma_eskn,
                'inputs': self.inputs}


class G0W0Outputs:
    def __init__(self, context, shape, ecut_e, sigma_eskn, dsigma_eskn,
                 eps_skn, vxc_skn, exx_skn, f_skn):
        self.extrapolate(context, shape, ecut_e, sigma_eskn, dsigma_eskn)
        self.Z_skn = 1 / (1 - self.dsigma_skn)

        # G0W0 single-step.
        # If we want GW0 again, we need to grab the expressions
        # from e.g. e73917fca5b9dc06c899f00b26a7c46e7d6fa749
        # or earlier and use qp correctly.
        self.qp_skn = eps_skn + self.Z_skn * (
            -vxc_skn + exx_skn + self.sigma_skn)

        self.sigma_eskn = sigma_eskn
        self.dsigma_eskn = dsigma_eskn

        self.eps_skn = eps_skn
        self.vxc_skn = vxc_skn
        self.exx_skn = exx_skn
        self.f_skn = f_skn

    def extrapolate(self, context, shape, ecut_e, sigma_eskn, dsigma_eskn):
        if len(ecut_e) == 1:
            self.sigma_skn = sigma_eskn[0]
            self.dsigma_skn = dsigma_eskn[0]
            self.sigr2_skn = None
            self.dsigr2_skn = None
            return

        from scipy.stats import linregress

        # Do linear fit of selfenergy vs. inverse of number of plane waves
        # to extrapolate to infinite number of plane waves

        context.print('', flush=False)
        context.print('Extrapolating selfenergy to infinite energy cutoff:',
                      flush=False)
        context.print('  Performing linear fit to %d points' % len(ecut_e))
        self.sigr2_skn = np.zeros(shape)
        self.dsigr2_skn = np.zeros(shape)
        self.sigma_skn = np.zeros(shape)
        self.dsigma_skn = np.zeros(shape)
        invN_i = ecut_e**(-3. / 2)
        for m in range(np.prod(shape)):
            s, k, n = np.unravel_index(m, shape)

            slope, intercept, r_value, p_value, std_err = \
                linregress(invN_i, sigma_eskn[:, s, k, n])

            self.sigr2_skn[s, k, n] = r_value**2
            self.sigma_skn[s, k, n] = intercept

            slope, intercept, r_value, p_value, std_err = \
                linregress(invN_i, dsigma_eskn[:, s, k, n])

            self.dsigr2_skn[s, k, n] = r_value**2
            self.dsigma_skn[s, k, n] = intercept

        if np.any(self.sigr2_skn < 0.9) or np.any(self.dsigr2_skn < 0.9):
            context.print('  Warning: Bad quality of linear fit for some ('
                          'n,k). ', flush=False)
            context.print('           Higher cutoff might be necessary.',
                          flush=False)

        context.print('  Minimum R^2 = %1.4f. (R^2 Should be close to 1)' %
                      min(np.min(self.sigr2_skn), np.min(self.dsigr2_skn)))

    def get_results_eV(self):
        results = {
            'f': self.f_skn,
            'eps': self.eps_skn * Ha,
            'vxc': self.vxc_skn * Ha,
            'exx': self.exx_skn * Ha,
            'sigma': self.sigma_skn * Ha,
            'dsigma': self.dsigma_skn,
            'Z': self.Z_skn,
            'qp': self.qp_skn * Ha}

        results.update(
            sigma_eskn=self.sigma_eskn * Ha,
            dsigma_eskn=self.dsigma_eskn)

        if self.sigr2_skn is not None:
            assert self.dsigr2_skn is not None
            results['sigr2_skn'] = self.sigr2_skn
            results['dsigr2_skn'] = self.dsigr2_skn

        return results


class QSymmetryOp:
    def __init__(self, symno, U_cc, sign):
        self.symno = symno
        self.U_cc = U_cc
        self.sign = sign

    def apply(self, q_c):
        return self.sign * (self.U_cc @ q_c)

    def check_q_Q_symmetry(self, Q_c, q_c):
        d_c = self.apply(q_c) - Q_c
        assert np.allclose(d_c.round(), d_c)

    def get_M_vv(self, cell_cv):
        # We'll be inverting these cells a lot.
        # Should have an object with the cell and its inverse which does this.
        return cell_cv.T @ self.U_cc.T @ np.linalg.inv(cell_cv).T

    @classmethod
    def get_symops(cls, qd, iq, q_c):
        # Loop over all k-points in the BZ and find those that are
        # related to the current IBZ k-point by symmetry
        Q1 = qd.ibz2bz_k[iq]
        done = set()
        for Q2 in qd.bz2bz_ks[Q1]:
            if Q2 >= 0 and Q2 not in done:
                time_reversal = qd.time_reversal_k[Q2]
                symno = qd.sym_k[Q2]
                Q_c = qd.bzk_kc[Q2]

                symop = cls(
                    symno=symno,
                    U_cc=qd.symmetry.op_scc[symno],
                    sign=1 - 2 * time_reversal)

                symop.check_q_Q_symmetry(Q_c, q_c)
                # Q_c, symop = QSymmetryOp.from_qd(qd, Q2, q_c)
                yield Q_c, symop
                done.add(Q2)

    @classmethod
    def get_symop_from_kpair(cls, kd, qd, kpt1, kpt2):
        # from k-point pair kpt1, kpt2 get Q_c = kpt2-kpt1, corrsponding IBZ
        # k-point q_c, indexes iQ, iq and symmetry transformation relating
        # Q_c to q_c
        Q_c = kd.bzk_kc[kpt2.K] - kd.bzk_kc[kpt1.K]
        iQ = qd.where_is_q(Q_c, qd.bzk_kc)
        iq = qd.bz2ibz_k[iQ]
        q_c = qd.ibzk_kc[iq]

        # Find symmetry that transforms Q_c into q_c
        sym = qd.sym_k[iQ]
        U_cc = qd.symmetry.op_scc[sym]
        time_reversal = qd.time_reversal_k[iQ]
        sign = 1 - 2 * time_reversal
        symop = cls(sym, U_cc, sign)
        symop.check_q_Q_symmetry(Q_c, q_c)

        return symop, iq

    def apply_symop_q(self, qpd, pawcorr, kpt1, kpt2):
        # returns necessary quantities to get symmetry transformed
        # density matrix
        Q_G = phase_shifted_fft_indices(kpt1.k_c, kpt2.k_c, qpd,
                                        coordinate_transformation=self.apply)

        qG_Gv = qpd.get_reciprocal_vectors(add_q=True)
        M_vv = self.get_M_vv(qpd.gd.cell_cv)
        mypawcorr = pawcorr.remap_by_symop(self, qG_Gv, M_vv)

        return mypawcorr, Q_G


def get_nmG(kpt1, kpt2, mypawcorr, n, qpd, I_G, pair_calc, timer=None):
    timer and timer.start('utcc and pawcorr multiply')
    ut1cc_R = kpt1.ut_nR[n].conj()
    C1_aGi = mypawcorr.multiply(kpt1.P_ani, band=n)
    timer and timer.stop('utcc and pawcorr multiply')
    n_mG = pair_calc.calculate_pair_density(
        ut1cc_R, C1_aGi, kpt2, qpd, I_G)
    return n_mG


gw_logo = """\
  ___  _ _ _
 |   || | | |
 | | || | | |
 |__ ||_____|
 |___|
"""


def get_max_nblocks(world, calc, ecut):
    nblocks = world.size
    if not isinstance(calc, (str, Path)):
        raise Exception('Using a calulator is not implemented at '
                        'the moment, load from file!')
        # nblocks_calc = calc
    else:
        nblocks_calc = GPAW(calc)
    ngmax = []
    for q_c in nblocks_calc.wfs.kd.bzk_kc:
        qpd = SingleQPWDescriptor.from_q(q_c, np.min(ecut) / Ha,
                                         nblocks_calc.wfs.gd)
        ngmax.append(qpd.ngmax)
    nG = np.min(ngmax)

    while nblocks > nG**0.5 + 1 or world.size % nblocks != 0:
        nblocks -= 1

    mynG = (nG + nblocks - 1) // nblocks
    assert mynG * (nblocks - 1) < nG
    return nblocks


def get_frequencies(frequencies, domega0, omega2):
    if domega0 is not None or omega2 is not None:
        assert frequencies is None
        frequencies = {'type': 'nonlinear',
                       'domega0': 0.025 if domega0 is None else domega0,
                       'omega2': 10.0 if omega2 is None else omega2}
        warnings.warn(f'Please use frequencies={frequencies}')
    elif frequencies is None:
        frequencies = {'type': 'nonlinear',
                       'domega0': 0.025,
                       'omega2': 10.0}
    else:
        assert frequencies['type'] == 'nonlinear'
    return frequencies


def choose_ecut_things(ecut, ecut_extrapolation):
    if ecut_extrapolation is True:
        pct = 0.8
        necuts = 3
        ecut_e = ecut * (1 + (1. / pct - 1) * np.arange(necuts)[::-1] /
                         (necuts - 1))**(-2 / 3)
    elif isinstance(ecut_extrapolation, (list, np.ndarray)):
        ecut_e = np.array(np.sort(ecut_extrapolation))
        ecut = ecut_e[-1]
    else:
        ecut_e = np.array([ecut])
    return ecut, ecut_e


def select_kpts(kpts, kd):
    """Function to process input parameters that take a list of k-points given
    in different format and returns a list of indices of the corresponding
    k-points in the IBZ."""

    if kpts is None:
        # Do all k-points in the IBZ:
        return np.arange(kd.nibzkpts)

    if np.asarray(kpts).ndim == 1:
        return kpts

    # Find k-points:
    bzk_Kc = kd.bzk_kc
    indices = []
    for k_c in kpts:
        d_Kc = bzk_Kc - k_c
        d_Kc -= d_Kc.round()
        K = abs(d_Kc).sum(1).argmin()
        if not np.allclose(d_Kc[K], 0):
            raise ValueError('Could not find k-point: {k_c}'
                             .format(k_c=k_c))
        k = kd.bz2ibz_k[K]
        indices.append(k)
    return indices


class G0W0Calculator:
    def __init__(self, filename='gw', *,
                 wd,
                 chi0calc,
                 wcalc,
                 kpts, bands, nbands=None,
                 fxc_modes,
                 eta,
                 ecut_e,
                 frequencies=None,
                 exx_vxc_calculator,
                 qcache,
                 ppa=False):
        """G0W0 calculator, initialized through G0W0 object.

        The G0W0 calculator is used to calculate the quasi
        particle energies through the G0W0 approximation for a number
        of states.

        Parameters
        ----------
        filename: str
            Base filename of output files.
        wcalc: WCalculator object
            Defines the calculator for computing the screened interaction
        kpts: list
            List of indices of the IBZ k-points to calculate the quasi particle
            energies for.
        bands:
            Range of band indices, like (n1, n2), to calculate the quasi
            particle energies for. Bands n where n1<=n<n2 will be
            calculated.  Note that the second band index is not included.
        frequencies:
            Input parameters for frequency_grid.
            Can be array of frequencies to evaluate the response function at
            or dictionary of parameters for build-in nonlinear grid
            (see :ref:`frequency grid`).
        ecut_e: array(float)
            Plane wave cut-off energies in eV. Defined with choose_ecut_things
        nbands: int
            Number of bands to use in the calculation. If None, the number will
            be determined from :ecut: to yield a number close to the number of
            plane waves used.
        do_GW_too: bool
            When carrying out a calculation including vertex corrections, it
            is possible to get the standard GW results at the same time
            (almost for free).
        ppa: bool
            Use Godby-Needs plasmon-pole approximation for screened interaction
            and self-energy
        """
        self.chi0calc = chi0calc
        self.wcalc = wcalc
        self.context = self.wcalc.context
        self.ppa = ppa
        self.qcache = qcache

        # Note: self.wd should be our only representation of the frequencies.
        # We should therefore get rid of self.frequencies.
        # It is currently only used by the restart code,
        # so should be easy to remove after some further adaptation.
        self.wd = wd
        self.frequencies = frequencies

        self.ecut_e = ecut_e / Ha

        self.context.print(gw_logo)

        self.fxc_modes = fxc_modes

        if self.fxc_modes[0] != 'GW':
            assert self.wcalc.xckernel.xc != 'RPA'

        if len(self.fxc_modes) == 2:
            # With multiple fxc_modes, we previously could do only
            # GW plus one other fxc_mode.  Now we can have any set
            # of modes, but whether things are consistent or not may
            # depend on how wcalc is configured.
            assert 'GW' in self.fxc_modes
            assert self.wcalc.xckernel.xc != 'RPA'

        self.filename = filename
        self.eta = eta / Ha

        self.kpts = kpts
        self.bands = bands

        b1, b2 = self.bands
        self.shape = (self.wcalc.gs.nspins, len(self.kpts), b2 - b1)

        self.nbands = nbands

        if self.wcalc.gs.nspins != 1:
            for fxc_mode in self.fxc_modes:
                if fxc_mode != 'GW':
                    raise RuntimeError('Including a xc kernel does not '
                                       'currently work for spin-polarized '
                                       f'systems. Invalid fxc_mode {fxc_mode}.'
                                       )

        self.pair_distribution = \
            self.chi0calc.kptpair_factory.distribute_k_points_and_bands(
                b1, b2, self.chi0calc.gs.kd.ibz2bz_k[self.kpts])

        self.print_parameters(kpts, b1, b2)

        self.exx_vxc_calculator = exx_vxc_calculator

        if self.ppa:
            self.context.print('Using Godby-Needs plasmon-pole approximation:')
            self.context.print('  Fitting energy: i*E0, E0 = %.3f Hartee'
                               % self.wd.omega_w[1].imag)
        else:
            self.context.print('Using full frequency integration')

    def print_parameters(self, kpts, b1, b2):
        isl = ['',
               'Quasi particle states:']
        if kpts is None:
            isl.append('All k-points in IBZ')
        else:
            kptstxt = ', '.join([f'{k:d}' for k in self.kpts])
            isl.append(f'k-points (IBZ indices): [{kptstxt}]')
        isl.extend([f'Band range: ({b1:d}, {b2:d})',
                    '',
                    'Computational parameters:'])
        if len(self.ecut_e) == 1:
            isl.append(
                'Plane wave cut-off: '
                f'{self.chi0calc.chi0_body_calc.ecut * Ha:g} eV')
        else:
            assert len(self.ecut_e) > 1
            isl.append('Extrapolating to infinite plane wave cut-off using '
                       'points at:')
            for ec in self.ecut_e:
                isl.append(f'  {ec * Ha:.3f} eV')
        isl.extend([f'Number of bands: {self.nbands:d}',
                    f'Coulomb cutoff: {self.wcalc.coulomb.truncation}',
                    f'Broadening: {self.eta * Ha:g} eV',
                    '',
                    f'fxc modes: {", ".join(sorted(self.fxc_modes))}',
                    f'Kernel: {self.wcalc.xckernel.xc}'])
        self.context.print('\n'.join(isl))

    def get_eps_and_occs(self):
        eps_skn = np.empty(self.shape)  # KS-eigenvalues
        f_skn = np.empty(self.shape)  # occupation numbers

        nspins = self.wcalc.gs.nspins
        b1, b2 = self.bands
        for i, k in enumerate(self.kpts):
            for s in range(nspins):
                u = s + k * nspins
                kpt = self.wcalc.gs.kpt_u[u]
                eps_skn[s, i] = kpt.eps_n[b1:b2]
                f_skn[s, i] = kpt.f_n[b1:b2] / kpt.weight

        return eps_skn, f_skn

    @timer('G0W0')
    def calculate(self):
        """Starts the G0W0 calculation.

        Returns a dict with the results with the following key/value pairs:

        ===========  =============================================
        key          value
        ===========  =============================================
        ``f``        Occupation numbers
        ``eps``      Kohn-Sham eigenvalues in eV
        ``vxc``      Exchange-correlation
                     contributions in eV
        ``exx``      Exact exchange contributions in eV
        ``sigma``    Self-energy contributions in eV
        ``dsigma``   Self-energy derivatives
        ``sigma_e``  Self-energy contributions in eV
                     used for ecut extrapolation
        ``Z``        Renormalization factors
        ``qp``       Quasi particle (QP) energies in eV
        ``iqp``      GW0/GW: QP energies for each iteration in eV
        ===========  =============================================

        All the values are ``ndarray``'s of shape
        (spins, IBZ k-points, bands)."""

        # Loop over q in the IBZ:
        self.context.print('Summing all q:')
        self.calculate_all_q_points()
        sigmas = self.read_sigmas()
        self.all_results = self.postprocess(sigmas)
        # Note: self.results is a pointer pointing to one of the results,
        # for historical reasons.

        self.savepckl()
        return self.results

    def postprocess(self, sigmas):
        all_results = {}
        for fxc_mode, sigma in sigmas.items():
            all_results[fxc_mode] = self.postprocess_single(fxc_mode, sigma)

        self.print_results(all_results)
        return all_results

    def read_sigmas(self):
        if self.context.comm.rank == 0:
            sigmas = self._read_sigmas()
        else:
            sigmas = None

        return broadcast(sigmas, comm=self.context.comm)

    def _read_sigmas(self):
        assert self.context.comm.rank == 0

        # Integrate over all q-points, and accumulate the quasiparticle shifts
        for iq, q_c in enumerate(self.wcalc.qd.ibzk_kc):
            key = str(iq)

            sigmas_contrib = self.get_sigmas_dict(key)

            if iq == 0:
                sigmas = sigmas_contrib
            else:
                for fxc_mode in self.fxc_modes:
                    sigmas[fxc_mode] += sigmas_contrib[fxc_mode]

        return sigmas

    def get_sigmas_dict(self, key):
        assert self.context.comm.rank == 0
        return {fxc_mode: Sigma.fromdict(sigma)
                for fxc_mode, sigma in self.qcache[key].items()}

    def postprocess_single(self, fxc_name, sigma):
        output = self.calculate_g0w0_outputs(sigma)
        return output.get_results_eV()

    def savepckl(self):
        """Save outputs to pckl files and return paths to those files."""
        # Note: this is always called, but the paths aren't returned
        # to the caller.  Calling it again then overwrites the files.
        #
        # TODO:
        #  * Replace with JSON
        #  * Save to different files or same file?
        #  * Move this functionality to g0w0 result object
        paths = {}
        for fxc_mode in self.fxc_modes:
            path = Path(f'{self.filename}_results_{fxc_mode}.pckl')
            with paropen(path, 'wb', comm=self.context.comm) as fd:
                pickle.dump(self.all_results[fxc_mode], fd, 2)
            paths[fxc_mode] = path

        # Do not return paths to caller before we know they all exist:
        self.context.comm.barrier()
        return paths

    def calculate_q(self, ie, k, kpt1, kpt2, qpd, Wdict,
                    *, symop, sigmas, blocks1d, pawcorr):
        """Calculates the contribution to the self-energy and its derivative
        for a given set of k-points, kpt1 and kpt2."""
        mypawcorr, I_G = symop.apply_symop_q(qpd, pawcorr, kpt1, kpt2)
        if debug:
            N_c = qpd.gd.N_c
            i_cG = symop.apply(np.unravel_index(qpd.Q_qG[0], N_c))
            bzk_kc = self.wcalc.gs.kd.bzk_kc
            Q_c = bzk_kc[kpt2.K] - bzk_kc[kpt1.K]
            shift0_c = Q_c - symop.apply(qpd.q_c)
            self.check(ie, i_cG, shift0_c, N_c, Q_c, mypawcorr)

        for n in range(kpt1.n2 - kpt1.n1):
            eps1 = kpt1.eps_n[n]
            n_mG = get_nmG(kpt1, kpt2, mypawcorr,
                           n, qpd, I_G, self.chi0calc.pair_calc)

            if symop.sign == 1:
                n_mG = n_mG.conj()

            f_m = kpt2.f_n
            deps_m = eps1 - kpt2.eps_n

            nn = kpt1.n1 + n - self.bands[0]

            assert set(Wdict) == set(sigmas)
            for fxc_mode in self.fxc_modes:
                sigma = sigmas[fxc_mode]
                Wmodel = Wdict[fxc_mode]
                # m is band index of all (both unoccupied and occupied) wave
                # functions in G
                for m, (deps, f, n_G) in enumerate(zip(deps_m, f_m, n_mG)):
                    # 2 * f - 1 will be used to select the branch of Hilbert
                    # transform, see get_HW of screened_interaction.py
                    # at FullFrequencyHWModel class.
                    S_GG, dSdw_GG = Wmodel.get_HW(deps, 2 * f - 1)
                    if S_GG is None:
                        continue

                    nc_G = n_G.conj()

                    # ie: ecut index for extrapolation
                    # kpt1.s: spin index of *
                    # k: k-point index of *
                    # nn: band index of *
                    # * wave function, where the sigma expectation value is
                    # evaluated
                    slot = ie, kpt1.s, k, nn
                    myn_G = n_G[blocks1d.myslice]
                    sigma.sigma_eskn[slot] += (myn_G @ S_GG @ nc_G).real
                    sigma.dsigma_eskn[slot] += (myn_G @ dSdw_GG @ nc_G).real

    def check(self, ie, i_cG, shift0_c, N_c, Q_c, pawcorr):
        # Can we delete this check? XXX
        assert np.allclose(shift0_c.round(), shift0_c)
        shift0_c = shift0_c.round().astype(int)
        I0_G = np.ravel_multi_index(i_cG - shift0_c[:, None], N_c, 'wrap')
        qpd = SingleQPWDescriptor.from_q(Q_c, self.ecut_e[ie],
                                         self.wcalc.gs.gd)
        G_I = np.empty(N_c.prod(), int)
        G_I[:] = -1
        I1_G = qpd.Q_qG[0]
        G_I[I1_G] = np.arange(len(I0_G))
        G_G = G_I[I0_G]
        # This indexing magic should definitely be moved to a method.
        # What on earth is it really?

        assert len(I0_G) == len(I1_G)
        assert (G_G >= 0).all()
        pairden_paw_corr = self.wcalc.gs.pair_density_paw_corrections
        pawcorr_wcalc1 = pairden_paw_corr(qpd)
        assert pawcorr.almost_equal(pawcorr_wcalc1, G_G)

    def calculate_all_q_points(self):
        """Main loop over irreducible Brillouin zone points.
        Handles restarts of individual qpoints using FileCache from ASE,
        and subsequently calls calculate_q."""

        pb = ProgressBar(self.context.fd)

        self.context.timer.start('W')
        self.context.print('\nCalculating screened Coulomb potential')
        self.context.print(self.wcalc.coulomb.description())

        chi0calc = self.chi0calc
        self.context.print(self.wd)

        # Find maximum size of chi-0 matrices:
        nGmax = max(count_reciprocal_vectors(chi0calc.chi0_body_calc.ecut,
                                             self.wcalc.gs.gd, q_c)
                    for q_c in self.wcalc.qd.ibzk_kc)
        nw = len(self.wd)

        size = self.chi0calc.chi0_body_calc.integrator.blockcomm.size

        mynGmax = (nGmax + size - 1) // size
        mynw = (nw + size - 1) // size

        # some memory sizes...
        if self.context.comm.rank == 0:
            siz = (nw * mynGmax * nGmax +
                   max(mynw * nGmax, nw * mynGmax) * nGmax) * 16
            sizA = (nw * nGmax * nGmax + nw * nGmax * nGmax) * 16
            self.context.print(
                '  memory estimate for chi0: local=%.2f MB, global=%.2f MB'
                % (siz / 1024**2, sizA / 1024**2))

        # Need to pause the timer in between iterations
        self.context.timer.stop('W')
        if self.context.comm.rank == 0:
            for key, sigmas in self.qcache.items():
                sigmas = {fxc_mode: Sigma.fromdict(sigma)
                          for fxc_mode, sigma in sigmas.items()}
                for fxc_mode, sigma in sigmas.items():
                    sigma.validate_inputs(self.get_validation_inputs())

        self.context.comm.barrier()
        for iq, q_c in enumerate(self.wcalc.qd.ibzk_kc):
            with ExitStack() as stack:
                if self.context.comm.rank == 0:
                    qhandle = stack.enter_context(self.qcache.lock(str(iq)))
                    skip = qhandle is None
                else:
                    skip = False

                skip = broadcast(skip, comm=self.context.comm)

                if skip:
                    continue

                result = self.calculate_q_point(iq, q_c, pb, chi0calc)

                if self.context.comm.rank == 0:
                    qhandle.save(result)
        pb.finish()

    def calculate_q_point(self, iq, q_c, pb, chi0calc):
        # Reset calculation
        sigmashape = (len(self.ecut_e), *self.shape)
        sigmas = {fxc_mode: Sigma(iq, q_c, fxc_mode, sigmashape,
                                  **self.get_validation_inputs())
                  for fxc_mode in self.fxc_modes}

        chi0 = chi0calc.create_chi0(q_c)

        m1 = chi0calc.gs.nocc1
        for ie, ecut in enumerate(self.ecut_e):
            self.context.timer.start('W')

            # First time calculation
            if ecut == chi0.qpd.ecut:
                # Nothing to cut away:
                m2 = self.nbands
            else:
                m2 = int(self.wcalc.gs.volume * ecut**1.5
                         * 2**0.5 / 3 / pi**2)
                if m2 > self.nbands:
                    raise ValueError(f'Trying to extrapolate ecut to'
                                     f'larger number of bands ({m2})'
                                     f' than there are bands '
                                     f'({self.nbands}).')
            qpdi, Wdict, blocks1d, pawcorr = self.calculate_w(
                chi0calc, q_c, chi0,
                m1, m2, ecut, iq)
            m1 = m2

            self.context.timer.stop('W')

            for nQ, (bzq_c, symop) in enumerate(QSymmetryOp.get_symops(
                    self.wcalc.qd, iq, q_c)):

                for (progress, kpt1, kpt2)\
                        in self.pair_distribution.kpt_pairs_by_q(bzq_c, 0, m2):
                    pb.update((nQ + progress) / self.wcalc.qd.mynk)

                    k1 = self.wcalc.gs.kd.bz2ibz_k[kpt1.K]
                    i = self.kpts.index(k1)

                    self.calculate_q(ie, i, kpt1, kpt2, qpdi, Wdict,
                                     symop=symop,
                                     sigmas=sigmas,
                                     blocks1d=blocks1d,
                                     pawcorr=pawcorr)

        for sigma in sigmas.values():
            sigma.sum(self.context.comm)

        return sigmas

    def get_validation_inputs(self):
        return {'kpts': self.kpts,
                'bands': list(self.bands),
                'nbands': self.nbands,
                'ecut_e': list(self.ecut_e),
                'frequencies': self.frequencies,
                'fxc_modes': self.fxc_modes,
                'integrate_gamma': self.wcalc.integrate_gamma}

    @timer('calculate_w')
    def calculate_w(self, chi0calc, q_c, chi0,
                    m1, m2, ecut,
                    iq):
        """Calculates the screened potential for a specified q-point."""

        chi0calc.chi0_body_calc.print_info(chi0.qpd)
        chi0calc.update_chi0(chi0, m1, m2, range(self.wcalc.gs.nspins))

        Wdict = {}

        for fxc_mode in self.fxc_modes:
            rqpd = chi0.qpd.copy_with(ecut=ecut)  # reduced qpd
            rchi0 = chi0.copy_with_reduced_pd(rqpd)
            Wdict[fxc_mode] = self.wcalc.get_HW_model(rchi0,
                                                      fxc_mode=fxc_mode)
            if (chi0calc.chi0_body_calc.pawcorr is not None and
                    rqpd.ecut < chi0.qpd.ecut):
                assert not self.ppa, """In previous master, PPA with ecut
                extrapolation was not working. Now it would work, but
                disabling it here still for sake of it is not tested."""

                pw_map = PWMapping(rqpd, chi0.qpd)
                # This is extremely bad behaviour! G0W0Calculator
                # should not change properties on the
                # Chi0BodyCalculator! Change in the future! XXX
                chi0calc.chi0_body_calc.pawcorr = \
                    chi0calc.chi0_body_calc.pawcorr.reduce_ecut(pw_map.G2_G1)

        # Create a blocks1d for the reduced plane-wave description
        blocks1d = Blocks1D(chi0.body.blockdist.blockcomm, rqpd.ngmax)

        return rqpd, Wdict, blocks1d, chi0calc.chi0_body_calc.pawcorr

    @timer('calcualte_vxc_and_exx')
    def calculate_vxc_and_exx(self):
        return self.exx_vxc_calculator.calculate(
            n1=self.bands[0], n2=self.bands[1],
            kpt_indices=self.kpts)

    def print_results(self, results):
        description = ['f:      Occupation numbers',
                       'eps:     KS-eigenvalues [eV]',
                       'vxc:     KS vxc [eV]',
                       'exx:     Exact exchange [eV]',
                       'sigma:   Self-energies [eV]',
                       'dsigma:  Self-energy derivatives',
                       'Z:       Renormalization factors',
                       'qp:      QP-energies [eV]']

        self.context.print('\nResults:')
        for line in description:
            self.context.print(line)

        b1, b2 = self.bands
        names = [line.split(':', 1)[0] for line in description]
        ibzk_kc = self.wcalc.gs.kd.ibzk_kc
        for s in range(self.wcalc.gs.nspins):
            for i, ik in enumerate(self.kpts):
                self.context.print(
                    '\nk-point ' + '{} ({}): ({:.3f}, {:.3f}, '
                    '{:.3f})'.format(i, ik, *ibzk_kc[ik]) +
                    '                ' + self.fxc_modes[0])
                self.context.print('band' + ''.join(f'{name:>8}'
                                                    for name in names))

                def actually_print_results(resultset):
                    for n in range(b2 - b1):
                        self.context.print(
                            f'{n + b1:4}' +
                            ''.join('{:8.3f}'.format(
                                resultset[name][s, i, n]) for name in names))

                for fxc_mode in results:
                    self.context.print(fxc_mode.rjust(69))
                    actually_print_results(results[fxc_mode])

        self.context.write_timer()

    def calculate_g0w0_outputs(self, sigma):
        eps_skn, f_skn = self.get_eps_and_occs()
        vxc_skn, exx_skn = self.calculate_vxc_and_exx()
        kwargs = dict(
            context=self.context,
            shape=self.shape,
            ecut_e=self.ecut_e,
            eps_skn=eps_skn,
            vxc_skn=vxc_skn,
            exx_skn=exx_skn,
            f_skn=f_skn)

        return G0W0Outputs(sigma_eskn=sigma.sigma_eskn,
                           dsigma_eskn=sigma.dsigma_eskn,
                           **kwargs)


def choose_bands(bands, relbands, nvalence, nocc):
    if bands is not None and relbands is not None:
        raise ValueError('Use bands or relbands!')

    if relbands is not None:
        bands = [nvalence // 2 + b for b in relbands]

    if bands is None:
        bands = [0, nocc]

    return bands


class G0W0(G0W0Calculator):
    def __init__(self, calc, filename='gw',
                 ecut=150.0,
                 ecut_extrapolation=False,
                 xc='RPA',
                 ppa=False,
                 E0=Ha,
                 eta=0.1,
                 nbands=None,
                 bands=None,
                 relbands=None,
                 frequencies=None,
                 domega0=None,  # deprecated
                 omega2=None,  # deprecated
                 nblocks=1,
                 nblocksmax=False,
                 kpts=None,
                 world=mpi.world,
                 timer=None,
                 fxc_mode='GW',
                 truncation=None,
                 integrate_gamma=0,
                 q0_correction=False,
                 do_GW_too=False,
                 **kwargs):
        """G0W0 calculator wrapper.

        The G0W0 calculator is used to calculate the quasi
        particle energies through the G0W0 approximation for a number
        of states.

        Parameters
        ----------
        calc:
            Filename of saved calculator object.
        filename: str
            Base filename of output files.
        kpts: list
            List of indices of the IBZ k-points to calculate the quasi particle
            energies for.
        bands:
            Range of band indices, like (n1, n2), to calculate the quasi
            particle energies for. Bands n where n1<=n<n2 will be
            calculated.  Note that the second band index is not included.
        relbands:
            Same as *bands* except that the numbers are relative to the
            number of occupied bands.
            E.g. (-1, 1) will use HOMO+LUMO.
        frequencies:
            Input parameters for frequency_grid.
            Can be an array of frequencies to evaluate the response function at
            or dictionary of parameters for build-in nonlinear grid
            (see :ref:`frequency grid`).
        ecut: float
            Plane wave cut-off energy in eV.
        ecut_extrapolation: bool or list
            If set to True an automatic extrapolation of the selfenergy to
            infinite cutoff will be performed based on three points
            for the cutoff energy.
            If an array is given, the extrapolation will be performed based on
            the cutoff energies given in the array.
        nbands: int
            Number of bands to use in the calculation. If None, the number will
            be determined from :ecut: to yield a number close to the number of
            plane waves used.
        ppa: bool
            Sets whether the Godby-Needs plasmon-pole approximation for the
            dielectric function should be used.
        xc: str
            Kernel to use when including vertex corrections.
        fxc_mode: str
            Where to include the vertex corrections; polarizability and/or
            self-energy. 'GWP': Polarizability only, 'GWS': Self-energy only,
            'GWG': Both.
        do_GW_too: bool
            When carrying out a calculation including vertex corrections, it
            is possible to get the standard GW results at the same time
            (almost for free).
        truncation: str
            Coulomb truncation scheme. Can be either 2D, 1D, or 0D.
        integrate_gamma: int
            Method to integrate the Coulomb interaction. 1 is a numerical
            integration at all q-points with G=[0,0,0] - this breaks the
            symmetry slightly. 0 is analytical integration at q=[0,0,0] only -
            this conserves the symmetry. integrate_gamma=2 is the same as 1,
            but the average is only carried out in the non-periodic directions.
        E0: float
            Energy (in eV) used for fitting in the plasmon-pole approximation.
        q0_correction: bool
            Analytic correction to the q=0 contribution applicable to 2D
            systems.
        nblocks: int
            Number of blocks chi0 should be distributed in so each core
            does not have to store the entire matrix. This is to reduce
            memory requirement. nblocks must be less than or equal to the
            number of processors.
        nblocksmax: bool
            Cuts chi0 into as many blocks as possible to reduce memory
            requirements as much as possible.
        """
        # We pass a serial communicator because the parallel handling
        # is somewhat wonky, we'd rather do that ourselves:
        try:
            qcache = FileCache(f'qcache_{filename}',
                               comm=mpi.SerialCommunicator())
        except TypeError as err:
            raise RuntimeError(
                'File cache requires ASE master '
                'from September 20 2022 or newer.  '
                'You may need to pull newest ASE.') from err
        if world.rank == 0:
            qcache.strip_empties()
        mode = 'a' if qcache.filecount() > 1 else 'w'

        frequencies = get_frequencies(frequencies, domega0, omega2)

        # (calc can not actually be a calculator at all.)
        gpwfile = Path(calc)

        context = ResponseContext(txt=filename + '.txt',
                                  comm=world, timer=timer)
        gs = ResponseGroundStateAdapter.from_gpw_file(gpwfile,
                                                      context=context)

        # Check if nblocks is compatible, adjust if not
        if nblocksmax:
            nblocks = get_max_nblocks(context.comm, gpwfile, ecut)

        kptpair_factory = KPointPairFactory(gs, context, nblocks=nblocks)

        kpts = list(select_kpts(kpts, gs.kd))

        if nbands is None:
            nbands = int(gs.volume * (ecut / Ha)**1.5 * 2**0.5 / 3 / pi**2)
        else:
            if ecut_extrapolation:
                raise RuntimeError(
                    'nbands cannot be supplied with ecut-extrapolation.')

        ecut, ecut_e = choose_ecut_things(ecut, ecut_extrapolation)

        if ppa:
            # use small imaginary frequency to avoid dividing by zero:
            frequencies = [1e-10j, 1j * E0]

            parameters = {'eta': 0,
                          'hilbert': False,
                          'timeordered': False}
        else:
            # frequencies = self.frequencies
            parameters = {'eta': eta,
                          'hilbert': True,
                          'timeordered': True}

        from gpaw.response.chi0 import new_frequency_descriptor
        wcontext = context.with_txt(filename + '.w.txt', mode=mode)
        wd = new_frequency_descriptor(gs, wcontext, nbands, frequencies)

        chi0calc = Chi0Calculator(
            wd=wd, kptpair_factory=kptpair_factory,
            nbands=nbands,
            ecut=ecut,
            intraband=False,
            context=wcontext,
            **parameters)

        bands = choose_bands(bands, relbands, gs.nvalence, chi0calc.gs.nocc2)

        coulomb = CoulombKernel.from_gs(gs, truncation=truncation)
        # XXX eta needs to be converted to Hartree here,
        # XXX and it is also converted to Hartree at superclass constructor
        # XXX called below. This needs to be cleaned up.
        wcalc = initialize_w_calculator(chi0calc, wcontext,
                                        ppa=ppa,
                                        xc=xc,
                                        E0=E0, eta=eta / Ha, coulomb=coulomb,
                                        integrate_gamma=integrate_gamma,
                                        q0_correction=q0_correction)

        fxc_modes = [fxc_mode]
        if do_GW_too:
            fxc_modes.append('GW')

        exx_vxc_calculator = EXXVXCCalculator(
            gpwfile,
            snapshotfile_prefix=filename)

        super().__init__(filename=filename,
                         wd=wd,
                         chi0calc=chi0calc,
                         wcalc=wcalc,
                         ecut_e=ecut_e,
                         eta=eta,
                         fxc_modes=fxc_modes,
                         nbands=nbands,
                         bands=bands,
                         frequencies=frequencies,
                         kpts=kpts,
                         exx_vxc_calculator=exx_vxc_calculator,
                         qcache=qcache,
                         ppa=ppa,
                         **kwargs)

    @property
    def results_GW(self):
        # Compatibility with old "do_GW_too" behaviour
        if 'GW' in self.fxc_modes and self.fxc_modes[0] != 'GW':
            return self.all_results['GW']

    @property
    def results(self):
        return self.all_results[self.fxc_modes[0]]


class EXXVXCCalculator:
    """EXX and Kohn-Sham XC contribution."""

    def __init__(self, gpwfile, snapshotfile_prefix):
        self._gpwfile = gpwfile
        self._snapshotfile_prefix = snapshotfile_prefix

    def calculate(self, n1, n2, kpt_indices):
        _, vxc_skn, exx_skn = non_self_consistent_eigenvalues(
            self._gpwfile,
            'EXX',
            n1, n2,
            kpt_indices=kpt_indices,
            snapshot=f'{self._snapshotfile_prefix}-vxc-exx.json',
        )
        return vxc_skn / Ha, exx_skn / Ha
