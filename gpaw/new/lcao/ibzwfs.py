from gpaw.new.ibzwfs import IBZWaveFunctions


class LCAOIBZWaveFunction(IBZWaveFunctions):
    def move(self, fracpos_ac, atomdist):
        from gpaw.new.lcao.builder import tci_helper

        super().move(fracpos_ac, atomdist)

        for wfs in self:
            basis = wfs.basis
            setups = wfs.setups
            break
        basis.set_positions(fracpos_ac)
        myM = (basis.Mmax + self.band_comm.size - 1) // self.band_comm.size
        basis.set_matrix_distribution(
            min(self.band_comm.rank * myM, basis.Mmax),
            min((self.band_comm.rank + 1) * myM, basis.Mmax))

        S_qMM, T_qMM, P_qaMi, tciexpansions, tci_derivatives = tci_helper(
            basis, self.ibz, self.domain_comm, self.band_comm, self.kpt_comm,
            fracpos_ac, atomdist,
            self.grid, self.dtype, setups)

        for wfs in self:
            wfs.tci_derivatives = tci_derivatives
            wfs.S_MM = S_qMM[wfs.q]
            wfs.T_MM = T_qMM[wfs.q]
            wfs.P_aMi = P_qaMi[wfs.q]
