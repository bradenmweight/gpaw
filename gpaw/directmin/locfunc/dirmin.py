from gpaw.directmin.fdpw.inner_loop import InnerLoop
from gpaw.directmin.fdpw.inner_loop_exst import InnerLoop as ILEXST
import numpy as np


class DirectMinLocalize:

    def __init__(self, obj_f, wfs, maxiter=300, g_tol=1.0e-3, randval=0.01):

        self.obj_f = obj_f
        self.randval = randval

        if obj_f.name == 'Zero':
            self.iloop = ILEXST(
                obj_f, wfs, 'all', 5.0e-4, maxiter,
                g_tol=g_tol, useprec=True)
        else:
            self.iloop = InnerLoop(
                obj_f, wfs, maxiter=maxiter,
                g_tol=g_tol)

    def run(self, wfs, dens, log=None, max_iter=None,
            g_tol=None, rewritepsi=True, ham=None, seed=None):

        if g_tol is not None:
            self.iloop.tol = g_tol
        if max_iter is not None:
            self.iloop.max_iter = max_iter
        # if log is None:
        #     log = parprint

        wfs.timer.start('Inner loop')

        if self.obj_f.name == 'Zero':
            etotal, counter = self.iloop.run(
                0.0, wfs, dens, log, 0,
                small_random=False,
                ham=ham)
        else:
            counter = self.iloop.run(0.0, wfs, dens, log, 0,
                                     randvalue=self.randval, seed=seed)

        if rewritepsi:
            for kpt in wfs.kpt_u:
                k = self.iloop.n_kps * kpt.s + kpt.q
                # n_occ = self.n_occ[k]
                dim1 = self.iloop.U_k[k].shape[0]
                if wfs.mode == 'lcao':
                    kpt.C_nM[:dim1] = \
                        self.iloop.U_k[k].T @ kpt.C_nM[:dim1]
                    wfs.atomic_correction.calculate_projections(
                        wfs, kpt)
                else:
                    kpt.psit_nG[:dim1] = np.tensordot(
                        self.iloop.U_k[k].T, kpt.psit_nG[:dim1], axes=1)
                    wfs.pt.integrate(kpt.psit_nG, kpt.P_ani, kpt.q)

        wfs.timer.stop('Inner loop')

        return counter
