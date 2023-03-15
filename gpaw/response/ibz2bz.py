import numpy as np


class IBZ2BZMapper:
    """Functionality to map data from k-points in the IBZ to the full BZ."""

    def __init__(self, kd, spos_ac, R_asii):
        """Construct the IBZ2BZMapper.

        Parameters
        ----------
        kd : KPointDescriptor
        spos_ac : np.array
            Scaled atomic positions
        R_asii : list
            Atomic symmetry rotations
        """
        self.kd = kd
        self.spos_ac = spos_ac
        self.R_asii = R_asii

    @classmethod
    def from_calculator(cls, calc):
        R_asii = [setup.R_sii for setup in calc.setups]
        return cls(calc.wfs.kd, calc.spos_ac, R_asii)

    def get_ik_c(self, K):
        ik = self.kd.bz2ibz_k[K]
        ik_c = self.kd.ibzk_kc[ik]
        return ik_c

    def get_rotation_matrix(self, K):
        """Coordinate rotation matrix, mapping IBZ -> K."""
        s = self.kd.sym_k[K]
        U_cc = self.kd.symmetry.op_scc[s]
        return U_cc

    def get_time_reversal(self, K):
        """Does the mapping IBZ -> K involve time reversal?"""
        time_reversal = self.kd.time_reversal_k[K]
        return time_reversal

    def get_atomic_rotation_matrices(self, K):
        """Atomic permutation and rotation involved in the IBZ -> K mapping.

        Returns
        -------
        b_a : list
            Atomic permutations (atom b is mapped onto atom a)
        U_aii : list
            Atomic rotation matrices for the PAW projections
        """
        s = self.kd.sym_k[K]
        U_cc = self.get_rotation_matrix(K)
        ik_c = self.get_ik_c(K)

        b_a = []
        U_aii = []
        for a, R_sii in enumerate(self.R_asii):
            b = self.kd.symmetry.a_sa[s, a]
            S_c = np.dot(self.spos_ac[a], U_cc) - self.spos_ac[b]
            x = np.exp(2j * np.pi * np.dot(ik_c, S_c))
            U_ii = R_sii[s].T * x
            b_a.append(b)
            U_aii.append(U_ii)

        return b_a, U_aii

    def map_kpoint(self, K):
        """Get the relative k-point coordinates after the IBZ -> K mapping.

        NB: The mapped k-point can lie outside the BZ, but will always be
        related to self.kd.bzk_kc[K] by a reciprocal lattice vector.
        """
        U_cc = self.get_rotation_matrix(K)
        time_reversal = self.get_time_reversal(K)

        # Apply symmetry operations to the irreducible k-point
        ik_c = self.get_ik_c(K)
        sign = 1 - 2 * time_reversal
        k_c = sign * U_cc @ ik_c

        return k_c

    def map_pseudo_wave(self, K, ut_R):
        """Map the periodic part of the pseudo wave from the IBZ -> K.

        The mapping takes place on the coarse real-space grid.

        NB: The k-point corresponding to the output ut_R does not necessarily
        lie within the BZ, see map_kpoint().
        """
        U_cc = self.get_rotation_matrix(K)
        time_reversal = self.get_time_reversal(K)

        # Apply symmetry operations to the periodic part of the pseudo wave
        if not (U_cc == np.eye(3)).all():
            N_c = ut_R.shape
            i_cr = np.dot(U_cc.T, np.indices(N_c).reshape((3, -1)))
            i = np.ravel_multi_index(i_cr, N_c, 'wrap')
            utout_R = ut_R.ravel()[i].reshape(N_c)
        else:
            utout_R = ut_R.copy()
        if time_reversal:
            utout_R = utout_R.conj()

        assert utout_R is not ut_R,\
            "We don't want the output array to point back at the input array"

        return utout_R

    def map_projections(self, K, projections):
        """Perform IBZ -> K mapping of the PAW projections.

        NB: The projections of atom b may be mapped onto *another* atom a.
        """
        time_reversal = self.get_time_reversal(K)
        b_a, U_aii = self.get_atomic_rotation_matrices(K)

        mapped_projections = projections.new()
        for a, (b, U_ii) in enumerate(zip(b_a, U_aii)):
            # Map projections
            Pin_ni = projections[b]
            Pout_ni = Pin_ni @ U_ii
            if time_reversal:
                Pout_ni = np.conj(Pout_ni)

            # Store output projections
            I1, I2 = mapped_projections.map[a]
            mapped_projections.array[..., I1:I2] = Pout_ni

        return mapped_projections
