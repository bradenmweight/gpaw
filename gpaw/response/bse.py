from dataclasses import dataclass
from datetime import timedelta
from functools import cached_property
from time import time, ctime

from ase.units import Hartree
from ase.dft import monkhorst_pack
import numpy as np
from scipy.linalg import eigh

from gpaw.blacs import BlacsGrid, Redistributor
from gpaw.kpt_descriptor import KPointDescriptor
from gpaw.mpi import world, serial_comm
from gpaw.response import ResponseContext
from gpaw.response.groundstate import CellDescriptor
from gpaw.response.chi0 import Chi0Calculator
from gpaw.response.context import timer
from gpaw.response.coulomb_kernels import CoulombKernel
from gpaw.response.df import write_response_function
from gpaw.response.frequencies import FrequencyDescriptor
from gpaw.response.pair import KPointPairFactory, get_gs_and_context
from gpaw.response.pair_functions import SingleQPWDescriptor
from gpaw.response.screened_interaction import initialize_w_calculator
from gpaw.response.g0w0 import validate_integrate_gamma


def decide_whether_tammdancoff(val_sn, con_sn):
    for n in val_sn[0]:
        if n in con_sn[0]:
            return False
    if len(val_sn) == 2:
        for n in val_sn[1]:
            if n in con_sn[1]:
                return False
    return True


@dataclass
class BSEMatrix:
    rhoG0_S: np.ndarray
    df_S: np.ndarray
    H_sS: np.ndarray

    def diagonalize_nontamdancoff(self, bse):
        df_S = self.df_S
        H_sS = self.H_sS
        rhoG0_S = self.rhoG0_S

        excludef_S = np.where(np.abs(df_S) < 0.001)[0]
        bse.context.print('  Using numpy.linalg.eig...')
        bse.context.print('  Eliminated %s pair orbitals' % len(
            excludef_S))

        H_SS = bse.collect_A_SS(H_sS)
        w_T = np.zeros(bse.nS - len(excludef_S), complex)
        if world.rank == 0:
            H_SS = np.delete(H_SS, excludef_S, axis=0)
            H_SS = np.delete(H_SS, excludef_S, axis=1)
            w_T, v_ST = np.linalg.eig(H_SS)
        world.broadcast(w_T, 0)

        # Here the eigenvectors are represented as complex conjugated rows

        df_S = np.delete(df_S, excludef_S)
        rhoG0_S = np.delete(rhoG0_S, excludef_S)

        C_T = np.zeros(bse.nS - len(excludef_S), complex)
        if world.rank == 0:
            A_T = np.dot(rhoG0_S, v_ST)
            B_T = np.dot(rhoG0_S * df_S, v_ST)
            tmp = np.dot(v_ST.conj().T, v_ST)
            overlap_tt = np.linalg.inv(tmp)
            C_T = np.dot(B_T.conj(), overlap_tt.T) * A_T
        world.broadcast(C_T, 0)
        return w_T, C_T

    def diagonalize_tamdancoff(self, bse):
        df_S = self.df_S
        H_sS = self.H_sS
        rhoG0_S = self.rhoG0_S

        nS = bse.nS
        ns = bse.ns

        if world.size == 1:
            bse.context.print('  Using lapack...')
            w_T, v_St = eigh(H_sS)
        else:
            bse.context.print('  Using scalapack...')
            assert ns == (
                -(-bse.kd.nbzkpts // world.size) * (
                    bse.nv * bse.nc *
                    bse.nspins *
                    (bse.spinors + 1)**2))

            # XXX We don't need to create new BLACS grids all the time
            # (also: remove the one further down)
            grid = BlacsGrid(world, world.size, 1)
            desc = grid.new_descriptor(nS, nS, ns, nS)

            desc2 = grid.new_descriptor(nS, nS, 2, 2)
            H_tmp = desc2.zeros(dtype=complex)
            r = Redistributor(world, desc, desc2)
            r.redistribute(H_sS, H_tmp)

            w_T = np.empty(nS)
            v_tmp = desc2.empty(dtype=complex)
            desc2.diagonalize_dc(H_tmp, v_tmp, w_T)

            r = Redistributor(grid.comm, desc2, desc)
            v_St = desc.zeros(dtype=complex)
            r.redistribute(v_tmp, v_St)
            v_St = v_St.conj().T

        A_t = np.dot(rhoG0_S, v_St)
        B_t = np.dot(rhoG0_S * df_S, v_St)

        if world.size == 1:
            C_T = B_t.conj() * A_t
        else:
            grid = BlacsGrid(world, world.size, 1)
            desc = grid.new_descriptor(nS, 1, ns, 1)
            C_t = desc.empty(dtype=complex)
            C_t[:, 0] = B_t.conj() * A_t
            C_T = desc.collect_on_master(C_t)[:, 0]
            if world.rank != 0:
                C_T = np.empty(nS, dtype=complex)
            world.broadcast(C_T, 0)
        return w_T, C_T


@dataclass
class ScreenedPotential:
    pawcorr_q: list
    W_qGG: list
    qpd_q: list


class SpinorData:
    def __init__(self, con_sn, val_sn, e_mk, v0_kmn, v1_kmn):
        self.e_mk = e_mk
        self.v0_kmn = v0_kmn
        self.v1_kmn = v1_kmn

        self.vi_s = [2 * val_sn[0, 0] - val_sn[0, -1] - 1]
        self.vf_s = [2 * con_sn[0, -1] - con_sn[0, 0] + 2]
        if self.vi_s[0] < 0:
            self.vi_s[0] = 0
        self.ci_s, self.cf_s = self.vi_s, self.vf_s
        self.ni, self.nf = self.vi_s[0], self.vf_s[0]

        self.valence_slice = slice(2 * val_sn[0, 0], 2 * (val_sn[0, -1] + 1))
        self.conduction_slice = slice(
            2 * con_sn[0, 0], 2 * (con_sn[0, -1] + 1))

    def _transform_rho(self, rho_mnG, K1, K2, slice1, slice2):
        nslice = slice(self.ni, self.nf)
        vec0_mn = self.v0_kmn[K1, slice1, nslice]
        vec1_mn = self.v1_kmn[K1, slice1, nslice]
        vec2_mn = self.v0_kmn[K2, slice2, nslice]
        vec3_mn = self.v1_kmn[K2, slice2, nslice]
        rho_0mnG = np.dot(vec0_mn.conj(), np.dot(vec2_mn, rho_mnG))
        rho_1mnG = np.dot(vec1_mn.conj(), np.dot(vec3_mn, rho_mnG))
        return rho_0mnG + rho_1mnG

    def rho_valence_valence(self, rho_mnG, K1, K2):
        return self._transform_rho(rho_mnG, K1, K2, self.valence_slice,
                                   self.valence_slice)

    def rho_conduction_conduction(self, rho_mnG, K1, K2):
        return self._transform_rho(rho_mnG, K1, K2, self.conduction_slice,
                                   self.conduction_slice)

    def rho_valence_conduction(self, rho_mnG, K1, K2):
        return self._transform_rho(rho_mnG, K1, K2, self.valence_slice,
                                   self.conduction_slice)

    def get_deps(self, K1, K2):
        epsv_m = self.e_mk[self.valence_slice, K1]
        epsc_n = self.e_mk[self.conduction_slice, K2]
        return -(epsv_m[:, np.newaxis] - epsc_n)


class BSEBackend:
    def __init__(self, *, gs, context,
                 valence_bands, conduction_bands,
                 spinors=False,
                 ecut=10.,
                 scale=1.0,
                 nbands=None,
                 eshift=None,
                 gw_skn=None,
                 truncation=None,
                 integrate_gamma='reciprocal',
                 mode='BSE',
                 q_c=[0.0, 0.0, 0.0],
                 direction=0):

        integrate_gamma = validate_integrate_gamma(integrate_gamma)

        self.gs = gs
        self.q_c = q_c
        self.direction = direction
        self.context = context

        self.spinors = spinors
        self.scale = scale

        assert mode in ['RPA', 'BSE']

        self.ecut = ecut / Hartree
        self.nbands = nbands
        self.mode = mode

        if integrate_gamma['type'] == 'sphere' and truncation is not None:
            self.context.print('***WARNING*** Analytical Coulomb integration' +
                               ' is not expected to work with Coulomb ' +
                               'truncation. ' +
                               'Use integrate_gamma=\'reciprocal\'')
        self.integrate_gamma = integrate_gamma

        # Find q-vectors and weights in the IBZ:
        self.kd = self.gs.kd
        if -1 in self.kd.bz2bz_ks:
            self.context.print('***WARNING*** Symmetries may not be right. '
                               'Use gamma-centered grid to be sure')
        offset_c = 0.5 * ((self.kd.N_c + 1) % 2) / self.kd.N_c
        bzq_qc = monkhorst_pack(self.kd.N_c) + offset_c
        self.qd = KPointDescriptor(bzq_qc)
        self.qd.set_symmetry(self.gs.atoms, self.kd.symmetry)

        # bands
        self.nspins = self.gs.nspins
        if self.nspins == 2:
            if self.spinors:
                self.spinors = False
                self.context.print('***WARNING*** Presently the spinor ' +
                                   'version does not work for spin-polarized' +
                                   ' calculations. Performing scalar ' +
                                   'calculation')

        self.val_sn = self.parse_bands(valence_bands, band_type='valence')
        self.con_sn = self.parse_bands(conduction_bands,
                                       band_type='conduction')

        self.use_tammdancoff = decide_whether_tammdancoff(self.val_sn,
                                                          self.con_sn)

        self.nv = len(self.val_sn[0])
        self.nc = len(self.con_sn[0])
        if eshift is not None:
            eshift /= Hartree
        if gw_skn is not None:
            assert self.nv + self.nc == len(gw_skn[0, 0])
            assert self.kd.nibzkpts == len(gw_skn[0])
            gw_skn = gw_skn[:, self.kd.bz2ibz_k]
            # assert self.kd.nbzkpts == len(gw_skn[0])
            gw_skn /= Hartree
        self.gw_skn = gw_skn
        self.eshift = eshift

        # Number of pair orbitals
        self.nS = self.kd.nbzkpts * self.nv * self.nc * self.nspins
        self.nS *= (self.spinors + 1)**2

        self.coulomb = CoulombKernel.from_gs(self.gs, truncation=truncation)
        self.context.print(self.coulomb.description())

        self.print_initialization(self.use_tammdancoff, self.eshift,
                                  self.gw_skn)

        self.Nv = self.nv * (self.spinors + 1)
        self.Nc = self.nc * (self.spinors + 1)
        self.ns = (-(-self.kd.nbzkpts // world.size)
                   * self.Nv * self.Nc * self.nspins)

        # Parallelization stuff
        self.nK = self.kd.nbzkpts
        self.myKrange, self.myKsize, self.mySsize = \
            self.parallelisation_sizes()

        # Setup bands
        if self.spinors:
            # Calculate spinors. Here m is index of eigenvalues with SOC
            # and n is the basis of eigenstates without SOC. Below m is used
            # for unoccupied states and n is used for occupied states so be
            # careful!
            self.spinors_data = self._spinordata()

            # Get all pair densities to allow for SOC mixing
            # Use twice as many no-SOC states as BSE bands to allow mixing
            # For example: 2 valence, 3 conduction, then
            # actually use 2 * 2 + 2 * 3 = 10 total bands
            # and then one calculates all matrix elements
            # (vv, vc, cv, and cc) in this 10x10 basis
            # from which they are then transformed into
            # to SOC basis.
            self.vi_s = self.spinors_data.vi_s
            self.vf_s = self.spinors_data.vf_s
            self.ci_s = self.spinors_data.ci_s
            self.cf_s = self.spinors_data.cf_s
        else:
            self.vi_s, self.vf_s = self.val_sn[:, 0], self.val_sn[:, -1] + 1
            self.ci_s, self.cf_s = self.con_sn[:, 0], self.con_sn[:, -1] + 1

    def parse_bands(self, bands, band_type='valence'):
        """Helper function that checks whether bands are correctly specified,
         and brings them to the format used later in the code.

        If the calculation is spin-polarized, band indices must
        be given explicitly as lists/arrays of shape (2,nbands) where the first
        index is for spin.

        If the calculation is not spin-polarized, either an integer (number of
        desired bands) or lists of band indices must be provided.

        band_type is an optional parameter that is only when a desired number
        of bands is given (rather than a list) to help figure out the correct
        band indices.
        """
        if hasattr(bands, '__iter__'):
            if self.nspins == 2:
                if len(bands) != 2 or (len(bands[0]) != len(bands[1])):
                    raise ValueError('For a spin-polarized calculation, '
                                     'the same number of bands must be '
                                     'specified for each spin! valence and '
                                     'conduction bands must be lists of shape '
                                     '(2,n)')

            bands_sn = np.atleast_2d(bands)
            return bands_sn

        # if we get here, bands is not iterable
        # check that the specified input is valid

        if self.nspins == 2:
            raise NotImplementedError('For a spin-polarized calculation, '
                                      'bands must be specified as lists '
                                      'of shape (2,n)')

        n_fully_occupied_bands, n_partially_occupied_bands = \
            self.gs.count_occupied_bands()

        if n_fully_occupied_bands != n_partially_occupied_bands:
            raise NotImplementedError('Automatic band generation is currently '
                                      'not implemented for metallic systems. '
                                      'Please specify band indices manually.')

        if band_type == 'valence':
            bands_sn = range(n_fully_occupied_bands - bands,
                             n_fully_occupied_bands)
        elif band_type == 'conduction':
            bands_sn = range(n_fully_occupied_bands,
                             n_fully_occupied_bands + bands)
        else:
            raise ValueError(f'Invalid band type: {band_type}')

        bands_sn = np.atleast_2d(bands_sn)
        return bands_sn

    def _spinordata(self):
        self.context.print('Diagonalizing spin-orbit Hamiltonian')
        soc = self.gs.soc_eigenstates(scale=self.scale)
        e_mk = soc.eigenvalues().T
        v_kmn = soc.eigenvectors()
        e_mk /= Hartree
        return SpinorData(self.con_sn, self.val_sn, e_mk,
                          v_kmn[:, :, ::2], v_kmn[:, :, 1::2])

    @timer('BSE calculate')
    def calculate(self, optical):
        # Calculate exchange interaction
        qpd0 = SingleQPWDescriptor.from_q(self.q_c, self.ecut, self.gs.gd)
        self.ikq_k = self.kd.find_k_plus_q(self.q_c)
        self.v_G = self.coulomb.V(qpd=qpd0, q_v=None)

        if optical:
            self.v_G[0] = 0.0

        kptpair_factory = KPointPairFactory(
            gs=self.gs,
            context=ResponseContext(txt='pair.txt', timer=self.context.timer,
                                    comm=serial_comm))
        pair_calc = kptpair_factory.pair_calculator()
        pawcorr = self.gs.pair_density_paw_corrections(qpd0)

        if self.mode != 'RPA':
            screened_potential = self.calculate_screened_potential()
        else:
            screened_potential = None

        # Calculate pair densities, eigenvalues and occupations
        self.context.timer.start('Pair densities')
        so = self.spinors + 1
        rhoex_KsmnG = np.zeros((self.nK, self.nspins, self.Nv,
                                self.Nc, len(self.v_G)), complex)
        df_Ksmn = np.zeros((self.nK, self.nspins, self.Nv,
                            self.Nc), float)  # -(ev - ec)
        deps_ksmn = np.zeros((self.myKsize, self.nspins, self.Nv,
                              self.Nc), float)  # -(fv - fc)

        optical_limit = np.allclose(self.q_c, 0.0)

        get_pair = kptpair_factory.get_kpoint_pair
        get_pair_density = pair_calc.get_pair_density

        # Calculate all properties diagonal in k-point
        # These include the indirect (exchange) kernel,
        # pseudo-energies, and occupation numbers
        for ik, iK in enumerate(self.myKrange):
            for s in range(self.nspins):
                pair = get_pair(qpd0, s, iK,
                                self.vi_s[s], self.vf_s[s],
                                self.ci_s[s], self.cf_s[s])
                m_m = np.arange(self.vi_s[s], self.vf_s[s])
                n_n = np.arange(self.ci_s[s], self.cf_s[s])
                if self.gw_skn is not None:
                    iKq = self.gs.kd.find_k_plus_q(self.q_c, [iK])[0]
                    epsv_m = self.gw_skn[s, iK, :self.nv]
                    epsc_n = self.gw_skn[s, iKq, self.nv:]
                    deps_ksmn[ik, s] = -(epsv_m[:, np.newaxis] - epsc_n)
                elif self.spinors:
                    iKq = self.gs.kd.find_k_plus_q(self.q_c, [iK])[0]
                    deps_ksmn[ik, s] = self.spinors_data.get_deps(iK, iKq)
                else:
                    deps_ksmn[ik, s] = -pair.get_transition_energies()

                rho_mnG = get_pair_density(qpd0, pair, m_m, n_n,
                                           pawcorr=pawcorr)
                if optical_limit:
                    n_mnv = pair_calc.get_optical_pair_density_head(
                        qpd0, pair, m_m, n_n)
                    rho_mnG[:, :, 0] = n_mnv[:, :, self.direction]
                if self.spinors:
                    if optical_limit:
                        deps0_mn = -pair.get_transition_energies()
                        rho_mnG[:, :, 0] *= deps0_mn

                    # This recreates the old behaviour of
                    # get_occupation_differences(self.val_sn[s],self.con_sn[s])
                    df_mn = (pair.kpt1.f_n[self.val_sn[s] -
                                           pair.kpt1.n1][:, np.newaxis] -
                             pair.kpt2.f_n[self.con_sn[s] - pair.kpt2.n1])

                    df_Ksmn[iK, s, ::2, ::2] = df_mn
                    df_Ksmn[iK, s, ::2, 1::2] = df_mn
                    df_Ksmn[iK, s, 1::2, ::2] = df_mn
                    df_Ksmn[iK, s, 1::2, 1::2] = df_mn

                    rhoex_KsmnG[iK, s] = \
                        self.spinors_data.rho_valence_conduction(
                        rho_mnG, iK, iKq)
                    if optical_limit:
                        rhoex_KsmnG[iK, s, :, :, 0] /= deps_ksmn[ik, s]
                else:
                    df_Ksmn[iK, s] = pair.get_occupation_differences()
                    rhoex_KsmnG[iK, s] = rho_mnG

        if self.eshift is not None:
            deps_ksmn[np.where(df_Ksmn[self.myKrange] > 1e-3)] += self.eshift
            deps_ksmn[np.where(df_Ksmn[self.myKrange] < -1e-3)] -= self.eshift

        world.sum(df_Ksmn)
        world.sum(rhoex_KsmnG)

        rhoG0_S = np.reshape(rhoex_KsmnG[:, :, :, :, 0], -1)
        self.context.timer.stop('Pair densities')

        # Calculate Hamiltonian
        self.context.timer.start('Calculate Hamiltonian')
        t0 = time()

        def update_progress(iK1):
            dt = time() - t0
            tleft = dt * self.myKsize / (iK1 + 1) - dt

            self.context.print(
                '  Finished %s pair orbitals in %s - Estimated %s left'
                % ((iK1 + 1) * self.Nv * self.Nc * self.nspins * world.size,
                    timedelta(seconds=round(dt)),
                    timedelta(seconds=round(tleft))))

        self.context.print('Calculating {} matrix elements at q_c = {}'.format(
            self.mode, self.q_c))

        # Hamiltonian buffer array
        H_ksmnKsmn = np.zeros((self.myKsize, self.nspins, self.Nv, self.Nc,
                               self.nK, self.nspins, self.Nv, self.Nc),
                              complex)

        # Add kernels to buffer array
        self.add_indirect_kernel(kptpair_factory, rhoex_KsmnG, H_ksmnKsmn)
        if self.mode != 'RPA':
            self.add_direct_kernel(kptpair_factory, pair_calc,
                                   screened_potential, update_progress,
                                   H_ksmnKsmn)

        H_ksmnKsmn /= self.gs.volume
        self.context.timer.stop('Calculate Hamiltonian')

        # XXX Why do we define a new mySsize?
        # is it different from self.mySsize,
        # from the tests it doesnt seem so.
        mySsize = self.myKsize * self.Nv * self.Nc * self.nspins
        if self.myKsize > 0:
            iS0 = self.myKrange[0] * self.Nv * self.Nc * self.nspins

        df_S = np.reshape(df_Ksmn, -1)
        # multiply by 2 when spin-paired and no SOC
        df_S *= 2.0 / self.nK / self.nspins / so
        deps_s = np.reshape(deps_ksmn, -1)
        H_sS = np.reshape(H_ksmnKsmn, (mySsize, self.nS))
        for iS in range(mySsize):
            # Multiply by occupations and adiabatic coupling
            H_sS[iS] *= df_S[iS0 + iS]
            # add bare transition energies
            H_sS[iS, iS0 + iS] += deps_s[iS]

        return BSEMatrix(rhoG0_S, df_S, H_sS)

    @timer('add_direct_kernel')
    def add_direct_kernel(self, kptpair_factory, pair_calc, screened_potential,
                          update_progress, H_ksmnKsmn):
        for ik1, iK1 in enumerate(self.myKrange):
            for s1 in range(self.nspins):
                kptv1 = kptpair_factory.get_k_point(
                    s1, iK1, self.vi_s[s1], self.vf_s[s1])
                kptc1 = kptpair_factory.get_k_point(
                    s1, self.ikq_k[iK1], self.ci_s[s1], self.cf_s[s1])

                for Q_c in self.qd.bzk_kc:
                    iK2 = self.kd.find_k_plus_q(Q_c, [kptv1.K])[0]

                    kptv2 = kptpair_factory.get_k_point(
                        s1, iK2, self.vi_s[s1], self.vf_s[s1])
                    kptc2 = kptpair_factory.get_k_point(
                        s1, self.ikq_k[iK2], self.ci_s[s1], self.cf_s[s1])

                    rho3_mmG, iq = self.get_density_matrix(
                        pair_calc, screened_potential, kptv1, kptv2)

                    rho4_nnG, iq = self.get_density_matrix(
                        pair_calc, screened_potential, kptc1, kptc2)

                    if self.spinors:
                        rho3_mmG = self.spinors_data.rho_valence_valence(
                            rho3_mmG, kptv1.K, kptv2.K)

                        rho4_nnG = self.spinors_data.rho_conduction_conduction(
                            rho4_nnG, kptc1.K, kptc2.K)

                    self.context.timer.start('Screened exchange')
                    W_mnmn = np.einsum(
                        'ijk,km,pqm->ipjq',
                        rho3_mmG.conj(),
                        screened_potential.W_qGG[iq],
                        rho4_nnG,
                        optimize='optimal')
                    W_mnmn *= self.nspins * (self.spinors + 1)
                    H_ksmnKsmn[ik1, s1, :, :, iK2, s1] -= 0.5 * W_mnmn
                    self.context.timer.stop('Screened exchange')

            if iK1 % (self.myKsize // 5 + 1) == 0:
                update_progress(iK1=iK1)

    @timer('add_indirect_kernel')
    def add_indirect_kernel(self, kptpair_factory, rhoex_KsmnG, H_ksmnKsmn):
        for ik1, iK1 in enumerate(self.myKrange):
            for s1 in range(self.nspins):
                kptv1 = kptpair_factory.get_k_point(
                    s1, iK1, self.vi_s[s1], self.vf_s[s1])
                rho1_mnG = rhoex_KsmnG[iK1, s1]
                # rhoex_KsnmG

                rho1ccV_mnG = rho1_mnG.conj()[:, :] * self.v_G
                for s2 in range(self.nspins):
                    for Q_c in self.qd.bzk_kc:
                        iK2 = self.kd.find_k_plus_q(Q_c, [kptv1.K])[0]
                        rho2_mnG = rhoex_KsmnG[iK2, s2]
                        self.context.timer.start('Coulomb')
                        H_ksmnKsmn[
                            ik1, s1, :, :, iK2, s2, :, :] += np.einsum(
                                'ijG,mnG->ijmn', rho1ccV_mnG, rho2_mnG,
                                optimize='optimal')
                        self.context.timer.stop('Coulomb')

    @timer('get_density_matrix')
    def get_density_matrix(self, pair_calc, screened_potential, kpt1, kpt2):
        self.context.timer.start('Symop')
        from gpaw.response.g0w0 import QSymmetryOp, get_nmG
        symop, iq = QSymmetryOp.get_symop_from_kpair(self.kd, self.qd,
                                                     kpt1, kpt2)
        qpd = screened_potential.qpd_q[iq]
        nG = qpd.ngmax
        pawcorr0 = screened_potential.pawcorr_q[iq]
        pawcorr, I_G = symop.apply_symop_q(qpd, pawcorr0, kpt1, kpt2)
        self.context.timer.stop('Symop')

        rho_mnG = np.zeros((len(kpt1.eps_n), len(kpt2.eps_n), nG), complex)
        for m in range(len(rho_mnG)):
            rho_mnG[m] = get_nmG(kpt1, kpt2, pawcorr, m, qpd, I_G,
                                 pair_calc, timer=self.context.timer)

        return rho_mnG, iq

    @cached_property
    def _chi0calc(self):
        return Chi0Calculator(
            self.gs, self.context.with_txt('chi0.txt'),
            wd=FrequencyDescriptor([0.0]),
            eta=0.001,
            ecut=self.ecut * Hartree,
            intraband=False,
            hilbert=False,
            nbands=self.nbands)

    @cached_property
    def blockcomm(self):
        return self._chi0calc.chi0_body_calc.blockcomm

    @cached_property
    def wcontext(self):
        return ResponseContext(txt='w.txt', comm=world)

    @cached_property
    def _wcalc(self):
        return initialize_w_calculator(
            self._chi0calc, self.wcontext,
            coulomb=self.coulomb,
            integrate_gamma=self.integrate_gamma)

    @timer('calculate_screened_potential')
    def calculate_screened_potential(self):
        """Calculate W_GG(q)."""

        pawcorr_q = []
        W_qGG = []
        qpd_q = []

        t0 = time()
        self.context.print('Calculating screened potential')
        for iq, q_c in enumerate(self.qd.ibzk_kc):
            chi0 = self._chi0calc.calculate(q_c)
            W_wGG = self._wcalc.calculate_W_wGG(chi0)
            W_GG = W_wGG[0]
            # This is such a terrible way to access the paw
            # corrections. Attributes should not be groped like
            # this... Change in the future! XXX
            pawcorr_q.append(self._chi0calc.chi0_body_calc.pawcorr)
            qpd_q.append(chi0.qpd)
            W_qGG.append(W_GG)

            if iq % (self.qd.nibzkpts // 5 + 1) == 2:
                dt = time() - t0
                tleft = dt * self.qd.nibzkpts / (iq + 1) - dt
                self.context.print(
                    '  Finished {} q-points in {} - Estimated {} left'.format(
                        iq + 1, timedelta(seconds=round(dt)), timedelta(
                            seconds=round(tleft))))

        return ScreenedPotential(pawcorr_q, W_qGG, qpd_q)

    @timer('diagonalize')
    def diagonalize_bse_matrix(self, bsematrix):
        self.context.print('Diagonalizing Hamiltonian')
        if self.use_tammdancoff:
            return bsematrix.diagonalize_tamdancoff(self)
        else:
            return bsematrix.diagonalize_nontamdancoff(self)

    @timer('get_bse_matrix')
    def get_bse_matrix(self, optical=True):
        """Calculate BSE matrix."""
        return self.calculate(optical=optical)

    @timer('get_vchi')
    def get_vchi(self, w_w=None, eta=0.1, optical=True, write_eig=None):
        """Returns v * chi where v is the bare Coulomb interaction"""

        bsematrix = self.get_bse_matrix(optical=optical)

        self.context.print('Calculating response function at %s frequency '
                           'points' % len(w_w))
        vchi_w = np.zeros(len(w_w), dtype=complex)

        w_T, C_T = self.diagonalize_bse_matrix(bsematrix)

        eta /= Hartree
        for iw, w in enumerate(w_w / Hartree):
            tmp_T = 1. / (w - w_T + 1j * eta)
            vchi_w[iw] += np.dot(tmp_T, C_T)
        vchi_w *= 4 * np.pi / self.gs.volume

        if not np.allclose(self.q_c, 0.0):
            cell_cv = self.gs.gd.cell_cv
            B_cv = 2 * np.pi * np.linalg.inv(cell_cv).T
            q_v = np.dot(self.q_c, B_cv)
            vchi_w /= np.dot(q_v, q_v)

        """Check f-sum rule."""
        nv = self.gs.nvalence
        dw_w = (w_w[1:] - w_w[:-1]) / Hartree
        wchi_w = (w_w[1:] * vchi_w[1:] + w_w[:-1] * vchi_w[:-1]) / Hartree / 2
        N = -np.dot(dw_w, wchi_w.imag) * self.gs.volume / (2 * np.pi**2)
        self.context.print('', flush=False)
        self.context.print('Checking f-sum rule:', flush=False)
        self.context.print(f'  Valence = {nv}, N = {N:f}', flush=False)
        self.context.print('')

        if write_eig is not None:
            assert isinstance(write_eig, str)
            filename = write_eig
            if world.rank == 0:
                write_bse_eigenvalues(filename, self.mode,
                                      w_T * Hartree, C_T)

        return vchi_w

    def get_dielectric_function(self, *args, filename='df_bse.csv', **kwargs):
        vchi = self.vchi(*args, optical=True, **kwargs)
        return vchi.dielectric_function(filename=filename)

    def get_eels_spectrum(self, *args, filename='df_bse.csv', **kwargs):
        vchi = self.vchi(*args, optical=False, **kwargs)
        return vchi.eels_spectrum(filename=filename)

    def get_polarizability(self, *args, filename='pol_bse.csv', **kwargs):
        # Previously it was
        # optical = (self.coulomb.truncation is None)
        # I.e. if a truncated kernel is used optical = False.
        # The reason it was set to False with Coulomb
        # truncation is that for q=0 V(G=0) is already
        # set to zero with the truncated coulomb Kernel.
        # However for finite q V(G=0) is different from zero.
        # Therefore the absorption spectra for 2D materials
        # calculated with the previous code was only correct for q=0.
        # See Issue #1055, the related MR and comments therein
        # For simplicity we set it to true for all cases here.
        vchi = self.vchi(*args, optical=True, **kwargs)
        return vchi.polarizability(filename=filename)

    def vchi(self, w_w=None, eta=0.1, write_eig='eig.dat',
             optical=True):
        vchi_w = self.get_vchi(w_w=w_w, eta=eta, optical=optical,
                               write_eig=write_eig)
        return VChi(self.gs.cd, self.context, w_w, vchi_w, optical=optical)

    def collect_A_SS(self, A_sS):
        if world.rank == 0:
            A_SS = np.zeros((self.nS, self.nS), dtype=complex)
            A_SS[:len(A_sS)] = A_sS
            Ntot = len(A_sS)
            for rank in range(1, world.size):
                nkr, nk, ns = self.parallelisation_sizes(rank)
                buf = np.empty((ns, self.nS), dtype=complex)
                world.receive(buf, rank, tag=123)
                A_SS[Ntot:Ntot + ns] = buf
                Ntot += ns
        else:
            world.send(A_sS, 0, tag=123)
        world.barrier()
        if world.rank == 0:
            return A_SS

    def parallelisation_sizes(self, rank=None):
        if rank is None:
            rank = world.rank
        nK = self.kd.nbzkpts
        myKsize = -(-nK // world.size)
        myKrange = range(rank * myKsize,
                         min((rank + 1) * myKsize, nK))
        myKsize = len(myKrange)
        mySsize = myKsize * self.nv * self.nc * self.nspins
        mySsize *= (1 + self.spinors)**2
        return myKrange, myKsize, mySsize

    def print_initialization(self, td, eshift, gw_skn):
        isl = ['----------------------------------------------------------',
               f'{self.mode} Hamiltonian',
               '----------------------------------------------------------',
               f'Started at:  {ctime()}', '',
               'Atoms                          : '
               f'{self.gs.atoms.get_chemical_formula(mode="hill")}',
               f'Ground state XC functional     : {self.gs.xcname}',
               f'Valence electrons              : {self.gs.nvalence}',
               f'Spinor calculations            : {self.spinors}',
               f'Number of bands                : {self.gs.bd.nbands}',
               f'Number of spins                : {self.gs.nspins}',
               f'Number of k-points             : {self.kd.nbzkpts}',
               f'Number of irreducible k-points : {self.kd.nibzkpts}',
               f'Number of q-points             : {self.qd.nbzkpts}',
               f'Number of irreducible q-points : {self.qd.nibzkpts}', '']

        for q in self.qd.ibzk_kc:
            isl.append(f'    q: [{q[0]:1.4f} {q[1]:1.4f} {q[2]:1.4f}]')
        isl.append('')
        if gw_skn is not None:
            isl.append('User specified BSE bands')
        isl.extend([f'Response PW cutoff             : {self.ecut * Hartree} '
                    f'eV',
                    f'Screening bands included       : {self.nbands}'])
        if len(self.val_sn) == 1:
            isl.extend([f'Valence bands                  : {self.val_sn[0]}',
                        f'Conduction bands               : {self.con_sn[0]}'])
        else:
            isl.extend([f'Valence bands                  : {self.val_sn[0]}'
                        f' {self.val_sn[1]}',
                        f'Conduction bands               : {self.con_sn[0]}'
                        f' {self.con_sn[1]}'])
        if eshift is not None:
            isl.append(f'Scissors operator              : {eshift * Hartree}'
                       f'eV')
        isl.extend([
            f'Tamm-Dancoff approximation     : {td}',
            f'Number of pair orbitals        : {self.nS}',
            '',
            f'Truncation of Coulomb kernel   : {self.coulomb.truncation}'])
        if self.integrate_gamma['type'] == 'sphere':
            isl.append(
                'Coulomb integration scheme     : Analytical - gamma only')
        elif self.integrate_gamma['type'] == 'reciprocal':
            isl.append(
                'Coulomb integration scheme     : Numerical - all q-points')
        else:
            isl.append(
                f'Coulomb integration scheme     : {self.integrate_gamma}')
        isl.extend([
            '',
            '----------------------------------------------------------',
            '----------------------------------------------------------',
            '',
            f'Parallelization - Total number of CPUs   : {world.size}',
            '  Screened potential',
            f'    K-point/band decomposition           : {world.size}',
            '  Hamiltonian',
            f'    Pair orbital decomposition           : {world.size}'])
        self.context.print('\n'.join(isl))


class BSE(BSEBackend):
    def __init__(self, calc=None, timer=None, txt='-', **kwargs):
        """Creates the BSE object

        calc: str or calculator object
            The string should refer to the .gpw file contaning KS orbitals
        ecut: float
            Plane wave cutoff energy (eV)
        nbands: int
            Number of bands used for the screened interaction
        valence_bands: list
            Valence bands used in the BSE Hamiltonian
        conduction_bands: list
            Conduction bands used in the BSE Hamiltonian
        eshift: float
            Scissors operator opening the gap (eV)
         q_c: list of three floats
-            Wavevector in reduced units on which the response is calculated
        direction: int
            if q_c = [0, 0, 0] this gives the direction in cartesian
            coordinates - 0=x, 1=y, 2=z
        gw_skn: list / array
            List or array defining the gw quasiparticle energies used in
            the BSE Hamiltonian. Should match spin, k-points and
            valence/conduction bands
        truncation: str or None
            Coulomb truncation scheme. Can be None or 2D.
        integrate_gamma: dict
        txt: str
            txt output
        mode: str
            Theory level used. can be RPA TDHF or BSE. Only BSE is screened.
        """
        gs, context = get_gs_and_context(
            calc, txt, world=world, timer=timer)

        super().__init__(gs=gs, context=context, **kwargs)


def write_bse_eigenvalues(filename, mode, w_w, C_w):
    with open(filename, 'w') as fd:
        print('# %s eigenvalues (in eV) and weights' % mode, file=fd)
        print('# Number   eig   weight', file=fd)
        for iw, (w, C) in enumerate(zip(w_w, C_w)):
            print('%8d %12.6f %12.16f' % (iw, w.real, C.real),
                  file=fd)


def read_bse_eigenvalues(filename):
    _, w_w, C_w = np.loadtxt(filename, unpack=True)
    return w_w, C_w


def write_spectrum(filename, w_w, A_w):
    with open(filename, 'w') as fd:
        for w, A in zip(w_w, A_w):
            print(f'{w:.9f}, {A:.9f}', file=fd)


def read_spectrum(filename):
    w_w, A_w = np.loadtxt(filename, delimiter=',',
                          unpack=True)
    return w_w, A_w


@dataclass
class VChi:
    cd: CellDescriptor
    context: ResponseContext
    w_w: np.ndarray
    vchi_w: np.ndarray
    optical: bool

    def epsilon(self):
        assert self.optical
        return -self.vchi_w + 1.0

    def eels(self):
        assert not self.optical
        return -self.vchi_w.imag

    def alpha(self):
        assert self.optical
        L = self.cd.nonperiodic_hypervolume
        return -L * self.vchi_w / (4 * np.pi)

    def dielectric_function(self, filename='df_bse.csv'):
        """Returns and writes real and imaginary part of the dielectric
        function.

        w_w: list of frequencies (eV)
            Dielectric function is calculated at these frequencies
        eta: float
            Lorentzian broadening of the spectrum (eV)
        filename: str
            data file on which frequencies, real and imaginary part of
            dielectric function is written
        write_eig: str
            File on which the BSE eigenvalues are written
        """

        return self._hackywrite(self.epsilon(), filename)

    # XXX The default filename clashes with that of dielectric function!
    def eels_spectrum(self, filename='df_bse.csv'):
        """Returns and writes real and imaginary part of the dielectric
        function.

        w_w: list of frequencies (eV)
            Dielectric function is calculated at these frequencies
        eta: float
            Lorentzian broadening of the spectrum (eV)
        filename: str
            data file on which frequencies, real and imaginary part of
            dielectric function is written
        write_eig: str
            File on which the BSE eigenvalues are written
        """
        return self._hackywrite(self.eels(), filename)

    def polarizability(self, filename='pol_bse.csv'):
        r"""Calculate the polarizability alpha.
        In 3D the imaginary part of the polarizability is related to the
        dielectric function by Im(eps_M) = 4 pi * Im(alpha). In systems
        with reduced dimensionality the converged value of alpha is
        independent of the cell volume. This is not the case for eps_M,
        which is ill defined. A truncated Coulomb kernel will always give
        eps_M = 1.0, whereas the polarizability maintains its structure.
        pbs should be a list of booleans giving the periodic directions.

        By default, generate a file 'pol_bse.csv'. The three colomns are:
        frequency (eV), Real(alpha), Imag(alpha). The dimension of alpha
        is \AA to the power of non-periodic directions.
        """
        return self._hackywrite(self.alpha(), filename)

    def _hackywrite(self, array, filename):
        if world.rank == 0 and filename is not None:
            if array.dtype == complex:
                write_response_function(filename, self.w_w, array.real,
                                        array.imag)
            else:
                assert array.dtype == float
                write_spectrum(filename, self.w_w, array)

        world.barrier()

        self.context.print('Calculation completed at:', ctime(), flush=False)
        self.context.print('')

        return self.w_w, array
