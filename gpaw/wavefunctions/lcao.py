import numpy as np
import scipy.sparse as sparse
from ase.units import Bohr

from gpaw.lfc import BasisFunctions
from gpaw.utilities import unpack
from gpaw.utilities.tools import tri2full
from gpaw import debug
from gpaw.lcao.overlap import NewTwoCenterIntegrals as NewTCI
from gpaw.utilities.blas import gemm, gemmdot
from gpaw.wavefunctions.base import WaveFunctions
from gpaw.lcao.atomic_correction import get_atomic_correction
from gpaw.wavefunctions.mode import Mode


class LCAO(Mode):
    name = 'lcao'

    def __init__(self, atomic_correction=None, interpolation=3,
                 force_complex_dtype=False):
        self.atomic_correction = atomic_correction
        self.interpolation = interpolation
        Mode.__init__(self, force_complex_dtype)

    def __call__(self, *args, **kwargs):
        return LCAOWaveFunctions(*args,
                                 atomic_correction=self.atomic_correction,
                                 **kwargs)

    def todict(self):
        dct = Mode.todict(self)
        dct['interpolation'] = self.interpolation
        return dct


# replace by class to make data structure perhaps a bit less confusing
def get_r_and_offsets(nl, spos_ac, cell_cv):
    r_and_offset_aao = {}

    def add(a1, a2, R_c, offset):
        if not (a1, a2) in r_and_offset_aao:
            r_and_offset_aao[(a1, a2)] = []
        r_and_offset_aao[(a1, a2)].append((R_c, offset))

    for a1, spos1_c in enumerate(spos_ac):
        a2_a, offsets = nl.get_neighbors(a1)
        for a2, offset in zip(a2_a, offsets):
            spos2_c = spos_ac[a2] + offset

            R_c = np.dot(spos2_c - spos1_c, cell_cv)
            add(a1, a2, R_c, offset)
            if a1 != a2 or offset.any():
                add(a2, a1, -R_c, -offset)

    return r_and_offset_aao


class LCAOWaveFunctions(WaveFunctions):
    mode = 'lcao'

    def __init__(self, ksl, gd, nvalence, setups, bd,
                 dtype, world, kd, kptband_comm, timer,
                 atomic_correction=None, collinear=True):
        WaveFunctions.__init__(self, gd, nvalence, setups, bd,
                               dtype, collinear, world, kd,
                               kptband_comm, timer)
        self.ksl = ksl
        self.S_qMM = None
        self.T_qMM = None
        self.P_aqMi = None
        self.debug_tci = False

        if atomic_correction is None:
            if ksl.using_blacs:
                atomic_correction = 'distributed'
            else:
                atomic_correction = 'dense'
        if isinstance(atomic_correction, str):
            atomic_correction = get_atomic_correction(atomic_correction)
        self.atomic_correction = atomic_correction

        self.timer.start('TCI: Evaluate splines')
        self.tci = NewTCI(gd.cell_cv, gd.pbc_c, setups, kd.ibzk_qc, kd.gamma)
        self.timer.stop('TCI: Evaluate splines')

        self.basis_functions = BasisFunctions(gd,
                                              [setup.phit_j
                                               for setup in setups],
                                              kd,
                                              dtype=dtype,
                                              cut=True)

    def set_orthonormalized(self, o):
        pass

    def empty(self, n=(), global_array=False, realspace=False):
        if realspace:
            return self.gd.empty(n, self.dtype, global_array)
        else:
            if isinstance(n, int):
                n = (n,)
            nao = self.setups.nao
            return np.empty(n + (nao,), self.dtype)

    def __str__(self):
        s = 'Wave functions: LCAO\n'
        s += '  Diagonalizer: %s\n' % self.ksl.get_description()
        s += '  Atomic Correction: %s\n' % self.atomic_correction.description
        s += '  Datatype: %s\n' % self.dtype.__name__
        return s

    def set_eigensolver(self, eigensolver):
        WaveFunctions.set_eigensolver(self, eigensolver)
        if eigensolver:
            eigensolver.initialize(self.gd, self.dtype, self.setups.nao,
                                   self.ksl)

    def set_positions(self, spos_ac, atom_partition=None, move_wfs=False):
        with self.timer('Basic WFS set positions'):
            WaveFunctions.set_positions(self, spos_ac, atom_partition)

        with self.timer('Basis functions set positions'):
            self.basis_functions.set_positions(spos_ac)

        if self.ksl is not None:
            self.basis_functions.set_matrix_distribution(self.ksl.Mstart,
                                                         self.ksl.Mstop)

        nq = len(self.kd.ibzk_qc)
        nao = self.setups.nao
        mynbands = self.bd.mynbands
        Mstop = self.ksl.Mstop
        Mstart = self.ksl.Mstart
        mynao = Mstop - Mstart

        #if self.ksl.using_blacs:  # XXX
            # S and T have been distributed to a layout with blacs, so
            # discard them to force reallocation from scratch.
            #
            # TODO: evaluate S and T when they *are* distributed, thus saving
            # memory and avoiding this problem
        for kpt in self.kpt_u:
            kpt.S_MM = None
            kpt.T_MM = None

        self.S_qMM = self.T_qMM = None
        self.P_aqMi = newP_aqMi = {}


        for kpt in self.kpt_u:
            if kpt.C_nM is None:
                kpt.C_nM = np.empty((mynbands, nao), self.dtype)


        if self.debug_tci:
            if self.ksl.using_blacs:
                self.tci.set_matrix_distribution(Mstart, mynao)
            oldS_qMM = np.empty((nq, mynao, nao), self.dtype)
            oldT_qMM = np.empty((nq, mynao, nao), self.dtype)


            oldP_aqMi = {}
            for a in self.basis_functions.my_atom_indices:
                ni = self.setups[a].ni
                oldP_aqMi[a] = np.empty((nq, nao, ni), self.dtype)

                # Calculate lower triangle of S and T matrices:
                self.tci.calculate(spos_ac, oldS_qMM, oldT_qMM,
                                   oldP_aqMi)

        from gpaw.lcao.tci import TCI
        # Establish arbitrary enumeration of setups:
        I_setup = {}
        setups_I = list(self.setups.setups.values())
        for I, setup in enumerate(setups_I):
            I_setup[setup] = I
        I_a = [I_setup[setup] for setup in self.setups]

        newtci = TCI([s.phit_j for s in setups_I],
                     [s.pt_j for s in setups_I],
                     I_a, self.gd.cell_cv, spos_ac=spos_ac,
                     pbc_c=self.gd.pbc_c, ibzk_qc=self.kd.ibzk_qc,
                     dtype=self.dtype)

        newS_qMM = np.zeros((nq, mynao, nao), self.dtype)
        newT_qMM = np.zeros((nq, mynao, nao), self.dtype)
        natoms = len(spos_ac)

        Mindices = self.setups.basis_indices()
        Pindices = self.setups.projector_indices()

        for a1 in self.basis_functions.my_atom_indices:
            I1, I2 = Pindices[a1]
            assert I2 - I1 == self.setups[a1].ni
            P_qMi = np.empty((nq, nao, I2 - I1), dtype=self.dtype)

            for a2 in range(natoms):
                N1, N2 = Mindices[a2]
                P_qmi = P_qMi[:, N1:N2, :]
                o = newtci.calculate(a1, a2, P=True)
                if o.P_qim is None:
                    P_qmi[:] = 0.0
                else:
                    P_qmi[:] = o.P_qim.transpose(0, 2, 1).conj()
            newP_aqMi[a1] = P_qMi

        for a1 in range(natoms):
            M1, M2 = Mindices[a1]
            if M2 < Mstart or M1 >= Mstop:  # > or >=?  Test this
                continue

            M1local = max(M1 - Mstart, 0)
            M2local = min(M2 - Mstart, mynao)
            nM = M2local - M1local

            if not nM:
                continue

            a2max = a1 + 1 if self.ksl.using_blacs else natoms

            for a2 in range(self.gd.comm.rank, a2max, self.gd.comm.size):
                o = newtci.calculate(a1, a2, OT=True)
                N1, N2 = Mindices[a2]
                if o.O_qmm is not None:
                    m1 = max(Mstart - M1, 0)
                    m2 = m1 + nM
                    newS_qMM[:, M1local:M2local, N1:N2] = o.O_qmm[:, m1:m2]
                    newT_qMM[:, M1local:M2local, N1:N2] = o.T_qmm[:, m1:m2]


        if self.debug_tci:
            assert len(oldP_aqMi) == len(newP_aqMi)
            Perr = 0.0
            for a in oldP_aqMi:
                P_err = max(Perr, np.abs(oldP_aqMi[a] - newP_aqMi[a]).max())
            assert Perr < 1e-15, Perr

        # TODO
        #   OK complex/conj, periodic images
        #   OK scalapack
        #   derivatives/forces
        #   sparse
        #   use symmetry/conj tricks to reduce calculations
        #   enable caching of spherical harmonics

        if self.atomic_correction.name != 'dense':
            from gpaw.lcao.newoverlap import newoverlap
            self.P_neighbors_a, self.P_aaqim = newoverlap(self, spos_ac)
        self.atomic_correction.gobble_data(self)

        self.atomic_correction.add_overlap_correction(self, newS_qMM)
        if self.debug_tci:
            self.atomic_correction.add_overlap_correction(self, oldS_qMM)

        if self.atomic_correction.implements_distributed_projections():
            my_atom_indices = self.atomic_correction.get_a_values()
        else:
            my_atom_indices = self.basis_functions.my_atom_indices
        self.allocate_arrays_for_projections(my_atom_indices)

        #S_MM = None  # allow garbage collection of old S_qMM after redist
        if self.debug_tci:
            oldS_qMM = self.ksl.distribute_overlap_matrix(oldS_qMM, root=-1)
            oldT_qMM = self.ksl.distribute_overlap_matrix(oldT_qMM, root=-1)

        newS_qMM = self.ksl.distribute_overlap_matrix(newS_qMM, root=-1)
        newT_qMM = self.ksl.distribute_overlap_matrix(newT_qMM, root=-1)

        #if (debug and self.bd.comm.size == 1 and self.gd.comm.rank == 0 and
        #    nao > 0 and not self.ksl.using_blacs):
            # S and T are summed only on comm master, so check only there
        #    from numpy.linalg import eigvalsh
        #    self.timer.start('Check positive definiteness')
        #    for S_MM in S_qMM:
        #        tri2full(S_MM, UL='L')
        #        smin = eigvalsh(S_MM).real.min()
        #        if smin < 0:
        #            raise RuntimeError('Overlap matrix has negative '
        #                               'eigenvalue: %e' % smin)
        #    self.timer.stop('Check positive definiteness')
        self.positions_set = True

        if self.debug_tci:
            Serr = np.abs(newS_qMM - oldS_qMM).max()
            Terr = np.abs(newT_qMM - oldT_qMM).max()
            print('S maxerr', Serr)
            print('T maxerr', Terr)
            assert Serr < 1e-15, Serr
            assert Terr < 1e-15, Terr

            assert len(oldP_aqMi) == len(newP_aqMi)
            for a in oldP_aqMi:
                Perr = np.abs(oldP_aqMi[a] - newP_aqMi[a]).max()
                assert Perr < 1e-15, (a, Perr)

        for kpt in self.kpt_u:
            q = kpt.q
            kpt.S_MM = newS_qMM[q]
            kpt.T_MM = newT_qMM[q]
        self.S_qMM = newS_qMM
        self.T_qMM = newT_qMM

    def initialize(self, density, hamiltonian, spos_ac):
        # Note: The above line exists also in set_positions.
        # This is guaranteed to be correct, but we can probably remove one.
        # Of course no human can understand the initialization process,
        # so this will be some other day.
        self.timer.start('LCAO WFS Initialize')
        if density.nt_sG is None:
            if self.kpt_u[0].f_n is None or self.kpt_u[0].C_nM is None:
                density.initialize_from_atomic_densities(self.basis_functions)
            else:
                # We have the info we need for a density matrix, so initialize
                # from that instead of from scratch.  This will be the case
                # after set_positions() during a relaxation
                density.initialize_from_wavefunctions(self)
            # Initialize GLLB-potential from basis function orbitals
            if hamiltonian.xc.type == 'GLLB':
                hamiltonian.xc.initialize_from_atomic_orbitals(
                    self.basis_functions)

        else:
            # After a restart, nt_sg doesn't exist yet, so we'll have to
            # make sure it does.  Of course, this should have been taken care
            # of already by this time, so we should improve the code elsewhere
            density.calculate_normalized_charges_and_mix()

        hamiltonian.update(density)
        self.timer.stop('LCAO WFS Initialize')

        return 0, 0

    def initialize_wave_functions_from_lcao(self):
        """Fill the calc.wfs.kpt_[u].psit_nG arrays with useful data.

        Normally psit_nG is NOT used in lcao mode, but some extensions
        (like ase.dft.wannier) want to have it.
        This code is adapted from fd.py / initialize_from_lcao_coefficients()
        and fills psit_nG with data constructed from the current lcao
        coefficients (kpt.C_nM).

        (This may or may not work in band-parallel case!)
        """
        from gpaw.wavefunctions.arrays import UniformGridWaveFunctions
        bfs = self.basis_functions
        for kpt in self.mykpts:
            kpt.psit = UniformGridWaveFunctions(
                self.bd.nbands, self.gd, self.dtype, kpt=kpt.q, dist=None,
                spin=kpt.s, collinear=True)
            kpt.psit_nG[:] = 0.0
            bfs.lcao_to_grid(kpt.C_nM, kpt.psit_nG[:self.bd.mynbands], kpt.q)

    def initialize_wave_functions_from_restart_file(self):
        """Dummy function to ensure compatibility to fd mode"""
        self.initialize_wave_functions_from_lcao()

    def add_orbital_density(self, nt_G, kpt, n):
        rank, u = self.kd.get_rank_and_index(kpt.s, kpt.k)
        assert rank == self.kd.comm.rank
        assert self.kpt_u[u] is kpt
        psit_G = self._get_wave_function_array(u, n, realspace=True)
        self.add_realspace_orbital_to_density(nt_G, psit_G)

    def calculate_density_matrix(self, f_n, C_nM, rho_MM=None):
        self.timer.start('Calculate density matrix')
        rho_MM = self.ksl.calculate_density_matrix(f_n, C_nM, rho_MM)
        self.timer.stop('Calculate density matrix')
        return rho_MM

        if 1:
            # XXX Should not conjugate, but call gemm(..., 'c')
            # Although that requires knowing C_Mn and not C_nM.
            # that also conforms better to the usual conventions in literature
            Cf_Mn = C_nM.T.conj() * f_n
            self.timer.start('gemm')
            gemm(1.0, C_nM, Cf_Mn, 0.0, rho_MM, 'n')
            self.timer.stop('gemm')
            self.timer.start('band comm sum')
            self.bd.comm.sum(rho_MM)
            self.timer.stop('band comm sum')
        else:
            # Alternative suggestion. Might be faster. Someone should test this
            from gpaw.utilities.blas import r2k
            C_Mn = C_nM.T.copy()
            r2k(0.5, C_Mn, f_n * C_Mn, 0.0, rho_MM)
            tri2full(rho_MM)

    def calculate_atomic_density_matrices_with_occupation(self, D_asp, f_un):
        ac = self.atomic_correction
        if ac.implements_distributed_projections():
            D2_asp = ac.redistribute(self, D_asp, type='asp', op='forth')
            WaveFunctions.calculate_atomic_density_matrices_with_occupation(
                self, D2_asp, f_un)
            D3_asp = ac.redistribute(self, D2_asp, type='asp', op='back')
            for a in D_asp:
                D_asp[a][:] = D3_asp[a]
        else:
            WaveFunctions.calculate_atomic_density_matrices_with_occupation(
                self, D_asp, f_un)

    def calculate_density_matrix_delta(self, d_nn, C_nM, rho_MM=None):
        self.timer.start('Calculate density matrix')
        rho_MM = self.ksl.calculate_density_matrix_delta(d_nn, C_nM, rho_MM)
        self.timer.stop('Calculate density matrix')
        return rho_MM

    def add_to_density_from_k_point_with_occupation(self, nt_sG, kpt, f_n):
        """Add contribution to pseudo electron-density. Do not use the standard
        occupation numbers, but ones given with argument f_n."""
        # Custom occupations are used in calculation of response potential
        # with GLLB-potential
        if kpt.rho_MM is None:
            rho_MM = self.calculate_density_matrix(f_n, kpt.C_nM)
            if hasattr(kpt, 'c_on'):
                assert self.bd.comm.size == 1
                d_nn = np.zeros((self.bd.mynbands, self.bd.mynbands),
                                dtype=kpt.C_nM.dtype)
                for ne, c_n in zip(kpt.ne_o, kpt.c_on):
                    assert abs(c_n.imag).max() < 1e-14
                    d_nn += ne * np.outer(c_n.conj(), c_n).real
                rho_MM += self.calculate_density_matrix_delta(d_nn, kpt.C_nM)
        else:
            rho_MM = kpt.rho_MM
        self.timer.start('Construct density')
        self.basis_functions.construct_density(rho_MM, nt_sG[kpt.s], kpt.q)
        self.timer.stop('Construct density')

    def add_to_kinetic_density_from_k_point(self, taut_G, kpt):
        raise NotImplementedError('Kinetic density calculation for LCAO '
                                  'wavefunctions is not implemented.')

    def calculate_forces(self, hamiltonian, F_av):
        self.timer.start('LCAO forces')

        ksl = self.ksl
        nao = ksl.nao
        mynao = ksl.mynao
        nq = len(self.kd.ibzk_qc)
        dtype = self.dtype
        tci = self.tci
        gd = self.gd
        bfs = self.basis_functions

        Mstart = ksl.Mstart
        Mstop = ksl.Mstop

        from gpaw.kohnsham_layouts import BlacsOrbitalLayouts
        isblacs = isinstance(ksl, BlacsOrbitalLayouts)  # XXX

        if not isblacs:
            self.timer.start('TCI derivative')
            dThetadR_qvMM = np.empty((nq, 3, mynao, nao), dtype)
            dTdR_qvMM = np.empty((nq, 3, mynao, nao), dtype)
            dPdR_aqvMi = {}
            for a in self.basis_functions.my_atom_indices:
                ni = self.setups[a].ni
                dPdR_aqvMi[a] = np.empty((nq, 3, nao, ni), dtype)
            tci.calculate_derivative(self.spos_ac, dThetadR_qvMM, dTdR_qvMM,
                                     dPdR_aqvMi)
            gd.comm.sum(dThetadR_qvMM)
            gd.comm.sum(dTdR_qvMM)
            self.timer.stop('TCI derivative')

            my_atom_indices = bfs.my_atom_indices
            atom_indices = bfs.atom_indices

            def _slices(indices):
                for a in indices:
                    M1 = bfs.M_a[a] - Mstart
                    M2 = M1 + self.setups[a].nao
                    if M2 > 0:
                        yield a, max(0, M1), M2

            def slices():
                return _slices(atom_indices)

            def my_slices():
                return _slices(my_atom_indices)

        dH_asp = hamiltonian.dH_asp
        vt_sG = hamiltonian.vt_sG

        #
        #         -----                    -----
        #          \    -1                  \    *
        # E      =  )  S     H    rho     =  )  c     eps  f  c
        #  mu nu   /    mu x  x z    z nu   /    n mu    n  n  n nu
        #         -----                    -----
        #          x z                       n
        #
        # We use the transpose of that matrix.  The first form is used
        # if rho is given, otherwise the coefficients are used.
        self.timer.start('Initial')

        rhoT_uMM = []
        ET_uMM = []

        if not isblacs:
            if self.kpt_u[0].rho_MM is None:
                self.timer.start('Get density matrix')
                for kpt in self.kpt_u:
                    rhoT_MM = ksl.get_transposed_density_matrix(kpt.f_n,
                                                                kpt.C_nM)
                    rhoT_uMM.append(rhoT_MM)
                    ET_MM = ksl.get_transposed_density_matrix(kpt.f_n *
                                                              kpt.eps_n,
                                                              kpt.C_nM)
                    ET_uMM.append(ET_MM)

                    if hasattr(kpt, 'c_on'):
                        # XXX does this work with BLACS/non-BLACS/etc.?
                        assert self.bd.comm.size == 1
                        d_nn = np.zeros((self.bd.mynbands, self.bd.mynbands),
                                        dtype=kpt.C_nM.dtype)
                        for ne, c_n in zip(kpt.ne_o, kpt.c_on):
                                d_nn += ne * np.outer(c_n.conj(), c_n)
                        rhoT_MM += ksl.get_transposed_density_matrix_delta(
                            d_nn, kpt.C_nM)
                        ET_MM += ksl.get_transposed_density_matrix_delta(
                            d_nn * kpt.eps_n, kpt.C_nM)
                self.timer.stop('Get density matrix')
            else:
                rhoT_uMM = []
                ET_uMM = []
                for kpt in self.kpt_u:
                    H_MM = self.eigensolver.calculate_hamiltonian_matrix(
                        hamiltonian, self, kpt)
                    tri2full(H_MM)
                    S_MM = kpt.S_MM.copy()
                    tri2full(S_MM)
                    ET_MM = np.linalg.solve(S_MM, gemmdot(H_MM,
                                                          kpt.rho_MM)).T.copy()
                    del S_MM, H_MM
                    rhoT_MM = kpt.rho_MM.T.copy()
                    rhoT_uMM.append(rhoT_MM)
                    ET_uMM.append(ET_MM)
        self.timer.stop('Initial')

        if isblacs:  # XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
            from gpaw.blacs import BlacsGrid, Redistributor

            def get_density_matrix(f_n, C_nM, redistributor):
                rho1_mm = ksl.calculate_blocked_density_matrix(f_n,
                                                               C_nM).conj()
                rho_mm = redistributor.redistribute(rho1_mm)
                return rho_mm

            pcutoff_a = [max([pt.get_cutoff() for pt in setup.pt_j])
                         for setup in self.setups]
            phicutoff_a = [max([phit.get_cutoff() for phit in setup.phit_j])
                           for setup in self.setups]

            # XXX should probably use bdsize x gdsize instead
            # That would be consistent with some existing grids
            grid = BlacsGrid(ksl.block_comm, self.gd.comm.size,
                             self.bd.comm.size)

            blocksize1 = -(-nao // grid.nprow)
            blocksize2 = -(-nao // grid.npcol)
            # XXX what are rows and columns actually?
            desc = grid.new_descriptor(nao, nao, blocksize1, blocksize2)

            rhoT_umm = []
            ET_umm = []
            redistributor = Redistributor(grid.comm, ksl.mmdescriptor, desc)
            Fpot_av = np.zeros_like(F_av)
            for u, kpt in enumerate(self.kpt_u):
                self.timer.start('Get density matrix')
                rhoT_mm = get_density_matrix(kpt.f_n, kpt.C_nM, redistributor)
                rhoT_umm.append(rhoT_mm)
                self.timer.stop('Get density matrix')

                self.timer.start('Potential')
                rhoT_mM = ksl.distribute_to_columns(rhoT_mm, desc)

                vt_G = vt_sG[kpt.s]
                Fpot_av += bfs.calculate_force_contribution(vt_G, rhoT_mM,
                                                            kpt.q)
                del rhoT_mM
                self.timer.stop('Potential')

            self.timer.start('Get density matrix')
            for kpt in self.kpt_u:
                ET_mm = get_density_matrix(kpt.f_n * kpt.eps_n, kpt.C_nM,
                                           redistributor)
                ET_umm.append(ET_mm)
            self.timer.stop('Get density matrix')

            M1start = blocksize1 * grid.myrow
            M2start = blocksize2 * grid.mycol

            M1stop = min(M1start + blocksize1, nao)
            M2stop = min(M2start + blocksize2, nao)

            m1max = M1stop - M1start
            m2max = M2stop - M2start

        if not isblacs:
            # Kinetic energy contribution
            #
            #           ----- d T
            #  a         \       mu nu
            # F += 2 Re   )   -------- rho
            #            /    d R         nu mu
            #           -----    mu nu
            #        mu in a; nu
            #
            Fkin_av = np.zeros_like(F_av)
            for u, kpt in enumerate(self.kpt_u):
                dEdTrhoT_vMM = (dTdR_qvMM[kpt.q] *
                                rhoT_uMM[u][np.newaxis]).real
                for a, M1, M2 in my_slices():
                    Fkin_av[a, :] += \
                        2.0 * dEdTrhoT_vMM[:, M1:M2].sum(-1).sum(-1)
            del dEdTrhoT_vMM

            # Density matrix contribution due to basis overlap
            #
            #            ----- d Theta
            #  a          \           mu nu
            # F  += -2 Re  )   ------------  E
            #             /        d R        nu mu
            #            -----        mu nu
            #         mu in a; nu
            #
            Ftheta_av = np.zeros_like(F_av)
            for u, kpt in enumerate(self.kpt_u):
                dThetadRE_vMM = (dThetadR_qvMM[kpt.q] *
                                 ET_uMM[u][np.newaxis]).real
                for a, M1, M2 in my_slices():
                    Ftheta_av[a, :] += \
                        -2.0 * dThetadRE_vMM[:, M1:M2].sum(-1).sum(-1)
            del dThetadRE_vMM

        if isblacs:
            from gpaw.lcao.overlap import TwoCenterIntegralCalculator
            self.timer.start('Prepare TCI loop')
            M_a = bfs.M_a

            Fkin2_av = np.zeros_like(F_av)
            Ftheta2_av = np.zeros_like(F_av)

            overlapcalc = TwoCenterIntegralCalculator(self.kd.ibzk_qc,
                                                      derivative=False)

            # XXX this is not parallel *AT ALL*.
            self.timer.start('Get neighbors')
            nl = tci.atompairs.pairs.neighbors
            r_and_offset_aao = get_r_and_offsets(nl, self.spos_ac, gd.cell_cv)
            atompairs = sorted(r_and_offset_aao.keys())
            self.timer.stop('Get neighbors')

            T_expansions = tci.T_expansions
            Theta_expansions = tci.Theta_expansions
            P_expansions = tci.P_expansions
            nq = len(self.kd.ibzk_qc)

            self.timer.start('broadcast dH')
            alldH_asp = {}
            for a in range(len(self.setups)):
                gdrank = bfs.sphere_a[a].rank
                if gdrank == gd.rank:
                    dH_sp = dH_asp[a]
                else:
                    ni = self.setups[a].ni
                    dH_sp = np.empty((self.nspins, ni * (ni + 1) // 2))
                gd.comm.broadcast(dH_sp, gdrank)
                # okay, now everyone gets copies of dH_sp
                alldH_asp[a] = dH_sp
            self.timer.stop('broadcast dH')

            # This will get sort of hairy.  We need to account for some
            # three-center overlaps, such as:
            #
            #         a1
            #      Phi   ~a3    a3  ~a3     a2     a2,a1
            #   < ----  |p  > dH   <p   |Phi  > rho
            #      dR
            #
            # To this end we will loop over all pairs of atoms (a1, a3),
            # and then a sub-loop over (a3, a2).
            from gpaw.lcao.overlap import DerivativeAtomicDisplacement

            class Displacement(DerivativeAtomicDisplacement):
                def __init__(self, a1, a2, R_c, offset):
                    phases = overlapcalc.phaseclass(overlapcalc.ibzk_qc,
                                                    offset)
                    DerivativeAtomicDisplacement.__init__(self, None, a1, a2,
                                                          R_c, offset, phases)

            # Cache of Displacement objects with spherical harmonics with
            # evaluated spherical harmonics.
            disp_aao = {}

            def get_displacements(a1, a2, maxdistance):
                # XXX the way maxdistance is handled it can lead to
                # bad caching when different maxdistances are passed
                # to subsequent calls with same pair of atoms
                disp_o = disp_aao.get((a1, a2))
                if disp_o is None:
                    disp_o = []
                    for R_c, offset in r_and_offset_aao[(a1, a2)]:
                        if np.linalg.norm(R_c) > maxdistance:
                            continue
                        disp = Displacement(a1, a2, R_c, offset)
                        disp_o.append(disp)
                    disp_aao[(a1, a2)] = disp_o
                return [disp for disp in disp_o if disp.r < maxdistance]

            self.timer.stop('Prepare TCI loop')
            self.timer.start('Not so complicated loop')

            for (a1, a2) in atompairs:
                if a1 >= a2:
                    # Actually this leads to bad load balance.
                    # We should take a1 > a2 or a1 < a2 equally many times.
                    # Maybe decide which of these choices
                    # depending on whether a2 % 1 == 0
                    continue

                m1start = M_a[a1] - M1start
                m2start = M_a[a2] - M2start
                if m1start >= blocksize1 or m2start >= blocksize2:
                    continue  # (we have only one block per CPU)

                T_expansion = T_expansions.get(a1, a2)
                Theta_expansion = Theta_expansions.get(a1, a2)
                nm1, nm2 = T_expansion.shape

                m1stop = min(m1start + nm1, m1max)
                m2stop = min(m2start + nm2, m2max)

                if m1stop <= 0 or m2stop <= 0:
                    continue

                m1start = max(m1start, 0)
                m2start = max(m2start, 0)
                J1start = max(0, M1start - M_a[a1])
                J2start = max(0, M2start - M_a[a2])
                M1stop = J1start + m1stop - m1start
                J2stop = J2start + m2stop - m2start

                dTdR_qvmm = T_expansion.zeros((nq, 3), dtype=dtype)
                dThetadR_qvmm = Theta_expansion.zeros((nq, 3), dtype=dtype)

                disp_o = get_displacements(a1, a2,
                                           phicutoff_a[a1] + phicutoff_a[a2])
                for disp in disp_o:
                    disp.evaluate_overlap(T_expansion, dTdR_qvmm)
                    disp.evaluate_overlap(Theta_expansion, dThetadR_qvmm)

                for u, kpt in enumerate(self.kpt_u):
                    rhoT_mm = rhoT_umm[u][m1start:m1stop, m2start:m2stop]
                    ET_mm = ET_umm[u][m1start:m1stop, m2start:m2stop]
                    Fkin_v = 2.0 * (dTdR_qvmm[kpt.q][:, J1start:M1stop,
                                                     J2start:J2stop] *
                                    rhoT_mm[np.newaxis]).real.sum(-1).sum(-1)
                    Ftheta_v = 2.0 * (dThetadR_qvmm[kpt.q][:, J1start:M1stop,
                                                           J2start:J2stop] *
                                      ET_mm[np.newaxis]).real.sum(-1).sum(-1)
                    Fkin2_av[a1] += Fkin_v
                    Fkin2_av[a2] -= Fkin_v
                    Ftheta2_av[a1] -= Ftheta_v
                    Ftheta2_av[a2] += Ftheta_v

            Fkin_av = Fkin2_av
            Ftheta_av = Ftheta2_av
            self.timer.stop('Not so complicated loop')

            dHP_and_dSP_aauim = {}

            a2values = {}
            for (a2, a3) in atompairs:
                if a3 not in a2values:
                    a2values[a3] = []
                a2values[a3].append(a2)

            Fatom_av = np.zeros_like(F_av)
            Frho_av = np.zeros_like(F_av)
            self.timer.start('Complicated loop')
            for a1, a3 in atompairs:
                if a1 == a3:
                    # Functions reside on same atom, so their overlap
                    # does not change when atom is displaced
                    continue
                m1start = M_a[a1] - M1start
                if m1start >= blocksize1:
                    continue

                P_expansion = P_expansions.get(a1, a3)
                nm1 = P_expansion.shape[0]
                m1stop = min(m1start + nm1, m1max)
                if m1stop <= 0:
                    continue

                m1start = max(m1start, 0)
                J1start = max(0, M1start - M_a[a1])
                J1stop = J1start + m1stop - m1start

                disp_o = get_displacements(a1, a3,
                                           phicutoff_a[a1] + pcutoff_a[a3])
                if len(disp_o) == 0:
                    continue

                dPdR_qvmi = P_expansion.zeros((nq, 3), dtype=dtype)
                for disp in disp_o:
                    disp.evaluate_overlap(P_expansion, dPdR_qvmi)

                dPdR_qvmi = dPdR_qvmi[:, :, J1start:J1stop, :].copy()
                for a2 in a2values[a3]:
                    m2start = M_a[a2] - M2start
                    if m2start >= blocksize2:
                        continue

                    P_expansion2 = P_expansions.get(a2, a3)
                    nm2 = P_expansion2.shape[0]
                    m2stop = min(m2start + nm2, m2max)
                    if m2stop <= 0:
                        continue

                    disp_o = get_displacements(a2, a3,
                                               phicutoff_a[a2] + pcutoff_a[a3])
                    if len(disp_o) == 0:
                        continue

                    m2start = max(m2start, 0)
                    J2start = max(0, M2start - M_a[a2])
                    J2stop = J2start + m2stop - m2start

                    if (a2, a3) in dHP_and_dSP_aauim:
                        dHP_uim, dSP_uim = dHP_and_dSP_aauim[(a2, a3)]
                    else:
                        P_qmi = P_expansion2.zeros((nq,), dtype=dtype)
                        for disp in disp_o:
                            disp.evaluate_direct(P_expansion2, P_qmi)
                        P_qmi = P_qmi[:, J2start:J2stop].copy()
                        dH_sp = alldH_asp[a3]
                        dS_ii = self.setups[a3].dO_ii

                        dHP_uim = []
                        dSP_uim = []
                        for u, kpt in enumerate(self.kpt_u):
                            dH_ii = unpack(dH_sp[kpt.s])
                            dHP_im = np.dot(P_qmi[kpt.q], dH_ii).T.conj()
                            # XXX only need nq of these
                            dSP_im = np.dot(P_qmi[kpt.q], dS_ii).T.conj()
                            dHP_uim.append(dHP_im)
                            dSP_uim.append(dSP_im)
                            dHP_and_dSP_aauim[(a2, a3)] = dHP_uim, dSP_uim

                    for u, kpt in enumerate(self.kpt_u):
                        rhoT_mm = rhoT_umm[u][m1start:m1stop, m2start:m2stop]
                        ET_mm = ET_umm[u][m1start:m1stop, m2start:m2stop]
                        dPdRdHP_vmm = np.dot(dPdR_qvmi[kpt.q], dHP_uim[u])
                        dPdRdSP_vmm = np.dot(dPdR_qvmi[kpt.q], dSP_uim[u])

                        Fatom_c = 2.0 * (dPdRdHP_vmm *
                                         rhoT_mm).real.sum(-1).sum(-1)
                        Frho_c = 2.0 * (dPdRdSP_vmm *
                                        ET_mm).real.sum(-1).sum(-1)
                        Fatom_av[a1] += Fatom_c
                        Fatom_av[a3] -= Fatom_c

                        Frho_av[a1] -= Frho_c
                        Frho_av[a3] += Frho_c

            self.timer.stop('Complicated loop')

        if not isblacs:
            # Potential contribution
            #
            #           -----      /  d Phi  (r)
            #  a         \        |        mu    ~
            # F += -2 Re  )       |   ---------- v (r)  Phi  (r) dr rho
            #            /        |     d R                nu          nu mu
            #           -----    /         a
            #        mu in a; nu
            #
            self.timer.start('Potential')
            Fpot_av = np.zeros_like(F_av)

            for u, kpt in enumerate(self.kpt_u):
                vt_G = vt_sG[kpt.s]
                Fpot_av += bfs.calculate_force_contribution(vt_G, rhoT_uMM[u],
                                                            kpt.q)
            self.timer.stop('Potential')

            # Density matrix contribution from PAW correction
            #
            #           -----                        -----
            #  a         \      a                     \     b
            # F +=  2 Re  )    Z      E        - 2 Re  )   Z      E
            #            /      mu nu  nu mu          /     mu nu  nu mu
            #           -----                        -----
            #           mu nu                    b; mu in a; nu
            #
            # with
            #                  b*
            #         -----  dP
            #   b      \       i mu    b   b
            #  Z     =  )   -------- dS   P
            #   mu nu  /     dR        ij  j nu
            #         -----    b mu
            #           ij
            #
            self.timer.start('Paw correction')
            Frho_av = np.zeros_like(F_av)
            for u, kpt in enumerate(self.kpt_u):
                work_MM = np.zeros((mynao, nao), dtype)
                ZE_MM = None
                for b in my_atom_indices:
                    setup = self.setups[b]
                    dO_ii = np.asarray(setup.dO_ii, dtype)
                    dOP_iM = np.zeros((setup.ni, nao), dtype)
                    gemm(1.0, self.P_aqMi[b][kpt.q], dO_ii, 0.0, dOP_iM, 'c')
                    for v in range(3):
                        gemm(1.0, dOP_iM,
                             dPdR_aqvMi[b][kpt.q][v][Mstart:Mstop],
                             0.0, work_MM, 'n')
                        ZE_MM = (work_MM * ET_uMM[u]).real
                        for a, M1, M2 in slices():
                            dE = 2 * ZE_MM[M1:M2].sum()
                            Frho_av[a, v] -= dE  # the "b; mu in a; nu" term
                            Frho_av[b, v] += dE  # the "mu nu" term
            del work_MM, ZE_MM
            self.timer.stop('Paw correction')

            # Atomic density contribution
            #            -----                         -----
            #  a          \     a                       \     b
            # F  += -2 Re  )   A      rho       + 2 Re   )   A      rho
            #             /     mu nu    nu mu          /     mu nu    nu mu
            #            -----                         -----
            #            mu nu                     b; mu in a; nu
            #
            #                  b*
            #         ----- d P
            #  b       \       i mu   b   b
            # A     =   )   ------- dH   P
            #  mu nu   /    d R       ij  j nu
            #         -----    b mu
            #           ij
            #
            self.timer.start('Atomic Hamiltonian force')
            Fatom_av = np.zeros_like(F_av)
            for u, kpt in enumerate(self.kpt_u):
                for b in my_atom_indices:
                    H_ii = np.asarray(unpack(dH_asp[b][kpt.s]), dtype)
                    HP_iM = gemmdot(H_ii,
                                    np.ascontiguousarray(
                                        self.P_aqMi[b][kpt.q].T.conj()))
                    for v in range(3):
                        dPdR_Mi = dPdR_aqvMi[b][kpt.q][v][Mstart:Mstop]
                        ArhoT_MM = (gemmdot(dPdR_Mi, HP_iM) * rhoT_uMM[u]).real
                        for a, M1, M2 in slices():
                            dE = 2 * ArhoT_MM[M1:M2].sum()
                            Fatom_av[a, v] += dE  # the "b; mu in a; nu" term
                            Fatom_av[b, v] -= dE  # the "mu nu" term
            self.timer.stop('Atomic Hamiltonian force')

        F_av += Fkin_av + Fpot_av + Ftheta_av + Frho_av + Fatom_av
        self.timer.start('Wait for sum')
        ksl.orbital_comm.sum(F_av)
        if self.bd.comm.rank == 0:
            self.kd.comm.sum(F_av, 0)
        self.timer.stop('Wait for sum')
        self.timer.stop('LCAO forces')

    def _get_wave_function_array(self, u, n, realspace=True, periodic=False):
        # XXX Taking kpt is better than taking u
        kpt = self.kpt_u[u]
        C_M = kpt.C_nM[n]

        if realspace:
            psit_G = self.gd.zeros(dtype=self.dtype)
            self.basis_functions.lcao_to_grid(C_M, psit_G, kpt.q)
            if periodic and self.dtype == complex:
                k_c = self.kd.ibzk_kc[kpt.k]
                return self.gd.plane_wave(-k_c) * psit_G
            return psit_G
        else:
            return C_M

    def write(self, writer, write_wave_functions=False):
        WaveFunctions.write(self, writer)
        if write_wave_functions:
            self.write_wave_functions(writer)

    def write_wave_functions(self, writer):
        writer.add_array(
            'coefficients',
            (self.nspins, self.kd.nibzkpts, self.bd.nbands, self.setups.nao),
            dtype=self.dtype)
        for s in range(self.nspins):
            for k in range(self.kd.nibzkpts):
                C_nM = self.collect_array('C_nM', k, s)
                writer.fill(C_nM * Bohr**-1.5)

    def read(self, reader):
        WaveFunctions.read(self, reader)
        r = reader.wave_functions
        if 'coefficients' in r:
            self.read_wave_functions(r)

    def read_wave_functions(self, reader):
        for kpt in self.kpt_u:
            C_nM = reader.proxy('coefficients', kpt.s, kpt.k)
            kpt.C_nM = self.bd.empty(self.setups.nao, dtype=self.dtype)
            for myn, C_M in enumerate(kpt.C_nM):
                n = self.bd.global_index(myn)
                # XXX number of bands could have been rounded up!
                if n >= len(C_nM):
                    break
                C_M[:] = C_nM[n] * Bohr**1.5

    def estimate_memory(self, mem):
        nq = len(self.kd.ibzk_qc)
        nao = self.setups.nao
        ni_total = sum([setup.ni for setup in self.setups])
        itemsize = mem.itemsize[self.dtype]
        mem.subnode('C [qnM]', nq * self.bd.mynbands * nao * itemsize)
        nM1, nM2 = self.ksl.get_overlap_matrix_shape()
        mem.subnode('S, T [2 x qmm]', 2 * nq * nM1 * nM2 * itemsize)
        mem.subnode('P [aqMi]', nq * nao * ni_total // self.gd.comm.size)
        self.tci.estimate_memory(mem.subnode('TCI'))
        self.basis_functions.estimate_memory(mem.subnode('BasisFunctions'))
        self.eigensolver.estimate_memory(mem.subnode('Eigensolver'),
                                         self.dtype)
