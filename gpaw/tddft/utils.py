# Written by Lauri Lehtovaara 2008

from gpaw.utilities.blas import axpy, scal
from gpaw.utilities.mblas import multi_axpy, multi_dotc, multi_dotu, multi_scal

import _gpaw
import gpaw.cuda


class MultiBlas:
    def __init__(self, gd, timer=None):
        self.gd = gd
        self.timer = None

    # Multivector ZAXPY: a x + y => y
    def multi_zaxpy(self, a, x, y):
        assert type(x) == type(y)
        if self.timer is not None:
            self.timer.start('Multi zaxpy')
        if isinstance(a, (float, complex)):
            axpy(complex(a), x, y)
        else:
            multi_axpy(a.astype(complex), x, y)
        if self.timer is not None:
            self.timer.stop('Multi zaxpy')

    # Multivector dot product, a^H b, where ^H is transpose
    def multi_zdotc(self, x, y, s=None):
        if self.timer is not None:
            self.timer.start('Multi zdotc')
        ss = multi_dotc(x, y, s)
        if self.gd.comm.size > 1:
            if isinstance(ss, gpaw.cuda.gpuarray.GPUArray):
                s_cpu = ss.get(pagelocked=True)
                self.gd.comm.sum(s_cpu)
                ss.set(s_cpu)
            else:
                self.gd.comm.sum(ss)
        if self.timer is not None:
            self.timer.stop('Multi zdotc')
        return ss

    def multi_zdotu(self, x, y, s = None):
        if self.timer is not None:
            self.timer.start('Multi zdotu')
        ss = multi_dotu(x, y, s)
        if self.gd.comm.size > 1:
            if isinstance(s, gpaw.cuda.gpuarray.GPUArray):
                s_cpu = ss.get(pagelocked=True)
                self.gd.comm.sum(s_cpu)
                ss.set(s_cpu)
            else:
                self.gd.comm.sum(ss)
        if self.timer is not None:
            self.timer.stop('Multi zdotu')
        return ss

    # Multiscale: a x => x
    def multi_scale(self, a, x):
        if self.timer is not None:
            self.timer.start('Multi scale')
        if isinstance(a, (float, complex)):
            scal(a,x)
        else:
            multi_scal(a,x)
        if self.timer is not None:
            self.timer.stop('Multi scale')

# -------------------------------------------------------------------

import numpy as np

class BandPropertyMonitor:
    def __init__(self, wfs, name, interval=1):
        self.niter = 0
        self.interval = interval

        self.wfs = wfs

        self.name = name

    def __call__(self):
        self.update(self.wfs)
        self.niter += self.interval

    def update(self, wfs):
        #strictly serial XXX!
        data_un = []

        for u, kpt in enumerate(wfs.kpt_u):
            data_n = getattr(kpt, self.name)

            data_un.append(data_n)

        self.write(np.array(data_un))

    def write(self, data):
        pass

class BandPropertyWriter(BandPropertyMonitor):
    def __init__(self, filename, wfs, name, interval=1):
        BandPropertyMonitor.__init__(self, wfs, name, interval)
        self.fileobj = open(filename,'w')

    def write(self, data):
        self.fileobj.write(data.tostring())
        self.fileobj.flush()

    def __del__(self):
        self.fileobj.close()

# -------------------------------------------------------------------


class StaticOverlapMonitor:
    def __init__(self, wfs, wf_u, P_aui, interval=1):
        self.niter = 0
        self.interval = interval

        self.wfs = wfs

        self.wf_u = wf_u
        self.P_aui = P_aui

    def __call__(self):
        self.update(self.wfs)
        self.niter += self.interval

    def update(self, wfs, calculate_P_ani=False):

        #strictly serial XXX!
        Porb_un = []

        for u, kpt in enumerate(wfs.kpt_u):
            swf = self.wf_u[u].ravel()

            psit_n = kpt.psit_nG.reshape((len(kpt.f_n),-1))
            Porb_n = np.dot(psit_n.conj(), swf) * wfs.gd.dv

            P_ani = kpt.P_ani

            if calculate_P_ani:
                #wfs.pt.integrate(psit_nG, P_ani, kpt.q)
                raise NotImplementedError(
                    'In case you were wondering, TODO XXX')

            for a, P_ni in P_ani.items():
                sP_i = self.P_aui[a][u]
                for n in range(wfs.bd.nbands):
                    for i in range(len(P_ni[0])):
                        for j in range(len(P_ni[0])):
                            Porb_n[n] += (P_ni[n][i].conj() *
                                       wfs.setups[a].dO_ii[i][j] *
                                       sP_i[j])

            Porb_un.append(Porb_n)

        self.write(np.array(Porb_un))

    def write(self, data):
        pass

class StaticOverlapWriter(StaticOverlapMonitor):
    def __init__(self, filename, wfs, overlap, interval=1):
        StaticOverlapMonitor.__init__(self, wfs, overlap, interval)
        self.fileobj = open(filename,'w')

    def write(self, data):
        self.fileobj.write(data.tostring())
        self.fileobj.flush()

    def __del__(self):
        self.fileobj.close()

# -------------------------------------------------------------------


class DynamicOverlapMonitor:
    def __init__(self, wfs, overlap, interval=1):
        self.niter = 0
        self.interval = interval

        self.setups = overlap.setups
        self.operator = overlap.operator
        self.wfs = wfs

    def __call__(self):
        self.update(self.wfs)
        self.niter += self.interval

    def update(self, wfs, calculate_P_ani=False):

        #strictly serial XXX!
        S_unn = []

        for kpt in wfs.kpt_u:
            psit_nG = kpt.psit_nG
            P_ani = kpt.P_ani

            if calculate_P_ani:
                #wfs.pt.integrate(psit_nG, P_ani, kpt.q)
                raise NotImplementedError(
                    'In case you were wondering, TODO XXX')

            # Construct the overlap matrix:
            S = lambda x: x
            dS_aii = dict([(a, self.setups[a].dO_ii) for a in P_ani])
            S_nn = self.operator.calculate_matrix_elements(psit_nG, P_ani,
                                                           S, dS_aii)
            S_unn.append(S_nn)

        self.write(np.array(S_unn))

    def write(self, data):
        pass

class DynamicOverlapWriter(DynamicOverlapMonitor):
    def __init__(self, filename, wfs, overlap, interval=1):
        DynamicOverlapMonitor.__init__(self, wfs, overlap, interval)
        self.fileobj = open(filename,'w')

    def write(self, data):
        self.fileobj.write(data.tostring())
        self.fileobj.flush()

    def __del__(self):
        self.fileobj.close()
