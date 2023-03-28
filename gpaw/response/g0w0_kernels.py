from ase.units import Ha
import numpy as np
from gpaw.xc.fxc import KernelWave, XCFlags, FXCCache
from gpaw.xc.rpa import GCut

from gpaw.response.pair_functions import SingleQPWDescriptor
from gpaw.pw.descriptor import PWMapping


class G0W0Kernel:
    def __init__(self, xc, context, **kwargs):
        self.xc = xc
        self.context = context
        self.xcflags = XCFlags(xc)
        self._kwargs = kwargs

    def calculate(self, qpd):
        if self.xc == 'RPA':
            return np.eye(qpd.ngmax)

        return calculate_spinkernel(
            qpd=qpd,
            xcflags=self.xcflags,
            context=self.context,
            **self._kwargs)


def calculate_spinkernel(*, ecut, xcflags, gs, qd, ns, qpd, context):
    assert xcflags.spin_kernel
    xc = xcflags.xc

    ibzq_qc = qd.ibzk_kc
    iq = np.argmin(np.linalg.norm(ibzq_qc - qpd.q_c[np.newaxis], axis=1))
    assert np.allclose(ibzq_qc[iq], qpd.q_c)

    ecut_max = ecut * Ha  # XXX very ugly this

    cache = FXCCache(tag=gs.atoms.get_chemical_formula(mode='hill'),
                     xc=xc, ecut=ecut_max)
    handle = cache.handle(iq)

    if not handle.exists():
        # Somehow we calculated many q even though this function
        # only works on one q?  Very confusing.
        kernel = KernelWave(
            q_empty=iq, ibzq_qc=qd.ibzk_kc,
            xc=xcflags.xc,
            ecut=ecut_max, gs=gs,
            context=context)

        for iq_calculated, array in kernel.calculate_fhxc():
            if context.comm.rank == 0:
                # assert array is not None  # XXXXXXX here
                cache.handle(iq_calculated).write(array)
            context.comm.barrier()

    context.comm.barrier()
    fv = handle.read()
    assert fv is not None

    # If we want a reduced plane-wave description, create qpd mapping
    if qpd.ecut < ecut:
        # Recreate nonreduced plane-wave description corresponding to ecut_max
        qpdnr = SingleQPWDescriptor.from_q(qpd.q_c, ecut, qpd.gd,
                                           gammacentered=qpd.gammacentered)
        pw_map = PWMapping(qpd, qpdnr)
        gcut = GCut(pw_map.G2_G1)
        fv = gcut.spin_cut(fv, ns=ns)

    return fv
