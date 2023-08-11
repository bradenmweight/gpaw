from math import pi

import numpy as np

from gpaw.core import UniformGrid
from gpaw.new import zips
from gpaw.new.pot_calc import PotentialCalculator


class UniformGridPotentialCalculator(PotentialCalculator):
    def __init__(self,
                 wf_grid: UniformGrid,
                 fine_grid: UniformGrid,
                 setups,
                 xc,
                 poisson_solver,
                 *,
                 nct_aR, nct_R,
                 tauct_aR=None, tauct_R=None,
                 interpolation_stencil_range=3,
                 xp=np):
        self.fine_grid = fine_grid
        self.grid = wf_grid
        self.nct_aR = nct_aR
        self.tauct_aR = tauct_aR

        fracpos_ac = nct_aR.fracpos_ac
        atomdist = nct_aR.atomdist

        self.vbar_ar = setups.create_local_potentials(fine_grid, fracpos_ac,
                                                      atomdist, xp=xp)
        self.ghat_aLr = setups.create_compensation_charges(fine_grid,
                                                           fracpos_ac,
                                                           atomdist,
                                                           xp=xp)

        self.vbar_r = fine_grid.empty(xp=xp)
        self.vbar_ar.to_uniform_grid(out=self.vbar_r)

        n = interpolation_stencil_range
        self.interpolation_stencil_range = n
        self.interpolate = wf_grid.transformer(fine_grid, n, xp=xp)
        self.restrict = fine_grid.transformer(wf_grid, n, xp=xp)

        super().__init__(xc, poisson_solver, setups,
                         nct_R=nct_R, tauct_R=tauct_R,
                         fracpos_ac=fracpos_ac)
        self.interpolation_domain = nct_aR.grid

    def __str__(self):
        txt = super().__str__()
        degree = self.interpolation_stencil_range * 2 - 1
        name = ['linear', 'cubic', 'quintic', 'heptic'][degree // 2]
        txt += (f'interpolation: tri-{name}' +
                f' # {degree}. degree polynomial\n')
        return txt

    def calculate_charges(self, vHt_r):
        return self.ghat_aLr.integrate(vHt_r)

    def calculate_non_selfconsistent_exc(self, xc, nt_sR, ibzwfs):
        nt_sr, _, _ = self._interpolate_density(nt_sR)
        vxct_sr = nt_sr.desc.zeros(nt_sr.dims)
        e_xc, _ = xc.calculate(nt_sr, vxct_sr, ibzwfs,
                               interpolate=self.interpolate,
                               restrict=self.restrict)
        return e_xc

    def _interpolate_density(self, nt_sR):
        nt_sr = self.interpolate(nt_sR)
        if not nt_sR.desc.pbc_c.all():
            Nt1_s = nt_sR.integrate()
            Nt2_s = nt_sr.integrate()
            for Nt1, Nt2, nt_r in zips(Nt1_s, Nt2_s, nt_sr):
                if Nt2 > 1e-14:
                    nt_r.data *= Nt1 / Nt2
        return nt_sr, None, None

    def calculate_pseudo_potential(self, density, ibzwfs, vHt_r):
        nt_sr, _, _ = self._interpolate_density(density.nt_sR)
        grid2 = nt_sr.desc

        if density.taut_sR is not None:
            taut_sr = self.interpolate(density.taut_sR)
        else:
            taut_sr = None

        e_xc, vxct_sr, dedtaut_sr = self.xc.calculate(nt_sr, taut_sr)

        charge_r = grid2.empty()
        charge_r.data[:] = nt_sr.data[:density.ndensities].sum(axis=0)
        e_zero = self.vbar_r.integrate(charge_r)

        ccc_aL = density.calculate_compensation_charge_coefficients()

        # Normalize: (LCAO basis functions may extend outside box)
        comp_charge = (4 * pi)**0.5 * sum(ccc_L[0]
                                          for ccc_L in ccc_aL.values())
        comp_charge = ccc_aL.layout.atomdist.comm.sum(comp_charge)
        pseudo_charge = charge_r.integrate()
        charge_r.data *= -(comp_charge + density.charge) / pseudo_charge

        self.ghat_aLr.add_to(charge_r, ccc_aL)

        if vHt_r is None:
            vHt_r = grid2.zeros()
        self.poisson_solver.solve(vHt_r, charge_r)
        e_coulomb = 0.5 * vHt_r.integrate(charge_r)

        e_kinetic = 0.0
        vt_sr = vxct_sr
        vt_sr.data += vHt_r.data + self.vbar_r.data
        vt_sR = self.restrict(vt_sr)
        for spin, (vt_R, nt_R) in enumerate(zips(vt_sR, density.nt_sR)):
            e_kinetic -= vt_R.integrate(nt_R)
            if spin < density.ndensities:
                e_kinetic += vt_R.integrate(self.nct_R)

        if dedtaut_sr is not None:
            dedtaut_sR = self.restrict(dedtaut_sr)
            for dedtaut_R, taut_R in zips(dedtaut_sR,
                                          density.taut_sR):
                e_kinetic -= dedtaut_R.integrate(taut_R)
                e_kinetic += dedtaut_R.integrate(self.tauct_R)

        e_external = 0.0

        return {'kinetic': e_kinetic,
                'coulomb': e_coulomb,
                'zero': e_zero,
                'xc': e_xc,
                'external': e_external}, vt_sR, dedtaut_sR, vHt_r

    def _move(self, fracpos_ac, atomdist, ndensities):
        self.ghat_aLr.move(fracpos_ac, atomdist)
        self.vbar_ar.move(fracpos_ac, atomdist)
        self.vbar_ar.to_uniform_grid(out=self.vbar_r)
        self.nct_aR.move(fracpos_ac, atomdist)
        self.nct_aR.to_uniform_grid(out=self.nct_R, scale=1.0 / ndensities)

    def force_contributions(self, state):
        density = state.density
        potential = state.potential
        nt_R = density.nt_sR[0]
        vt_R = potential.vt_sR[0]
        if density.ndensities > 1:
            nt_R = nt_R.desc.empty()
            nt_R.data[:] = density.nt_sR.data[:density.ndensities].sum(axis=0)
            vt_R = vt_R.desc.empty()
            vt_R.data[:] = (
                potential.vt_sR.data[:density.ndensities].sum(axis=0) /
                density.ndensities)

        nt_r = self.interpolate(nt_R)
        if not nt_r.desc.pbc_c.all():
            scale = nt_R.integrate() / nt_r.integrate()
            nt_r.data *= scale

        return (self.ghat_aLr.derivative(state.vHt_x),
                self.nct_aR.derivative(vt_R),
                self.vbar_ar.derivative(nt_r))
