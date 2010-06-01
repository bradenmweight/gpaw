import numpy as np
from ase.dft.kpoints import get_monkhorst_shape

def find_kq(bzk_kv, q):
    """Find the index of k+q for all kpoints in BZ."""

    nkpt = bzk_kv.shape[0]
    kq = np.zeros(nkpt, dtype=int)
    nkptxyz = get_monkhorst_shape(bzk_kv)
        
    dk = 1. / nkptxyz 
    kmax = (nkptxyz - 1) * dk / 2.
    N = np.zeros(3, dtype=int)

    for k in range(nkpt):
        kplusq = bzk_kv[k] + q
        for dim in range(3):
            if kplusq[dim] > 0.5: # 
                kplusq[dim] -= 1.
            elif kplusq[dim] < -0.5:
                kplusq[dim] += 1.

            N[dim] = int(np.round((kplusq[dim] + kmax[dim])/dk[dim]))

        kq[k] = N[2] + N[1] * nkptxyz[2] + N[0] * nkptxyz[2]* nkptxyz[1]

        tmp = bzk_kv[kq[k]]
        if (abs(kplusq - tmp)).sum() > 1e-8:
            print k, kplusq, tmp
            raise ValueError('k+q index not correct!')

    return kq


def find_ibzkpt(symrel, ibzk_kv, bzk_v):
    """Given a certain kpoint, find its index in IBZ and related symmetry operations."""
    find = False
    ibzkpt = 0
    iop = 0
    timerev = False

    for ioptmp, op in enumerate(symrel):
        for i, ibzk in enumerate(ibzk_kv):
            diff_c = np.dot(ibzk, op.T) - bzk_v
            if (np.abs(diff_c - diff_c.round()) < 1e-8).all():
                ibzkpt = i
                iop = ioptmp
                find = True
                break

            diff_c = np.dot(ibzk, op.T) + bzk_v
            if (np.abs(diff_c - diff_c.round()) < 1e-8).all():
                ibzkpt = i
                iop = ioptmp
                find = True
                timerev = True
                break
            
        if find == True:
            break
        
    if find == False:        
        print bzk_v
        print ibzk_kv
        raise ValueError('Cant find corresponding IBZ kpoint!')

    return ibzkpt, iop, timerev


def symmetrize_wavefunction(a_g, op_cc, kpt0, kpt1, timerev):

    if (np.abs(op_cc - np.eye(3,dtype=int)) < 1e-10).all():
        if timerev:
            return a_g.conj()
        else:
            return a_g
    elif (np.abs(op_cc + np.eye(3,dtype=int)) < 1e-10).all():
        return a_g.conj()
    else:
        import _gpaw
        b_g = np.zeros_like(a_g)
        if timerev:
            _gpaw.symmetrize_wavefunction(a_g, b_g, op_cc.T.copy(), kpt0, -kpt1)
            return b_g.conj()
        else:
            _gpaw.symmetrize_wavefunction(a_g, b_g, op_cc.T.copy(), kpt0, kpt1)
            return b_g
