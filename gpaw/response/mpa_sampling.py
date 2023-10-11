import numpy as np

def sampling_branches(w_dist, ps='2l', varpi=1, shift0=1.e-5, shift=0.1):

    if len(w_dist) == 1:  # the value of ps is irrelevant in the case of a single pole
        assert (shift0 >= 0)
        assert (varpi > shift0)
        w_grid = np.array([w_dist + 1j*shift0, w_dist + 1j*varpi], dtype=complex)
        return w_grid

    if ps == '1l':  # only one branch
        assert (varpi >= 0)
        w_grid = np.array(w_dist + 1j*varpi, dtype=complex)
        return w_grid

    if ps == '2l':  # two branches
        assert (shift0 >= 0)
        assert (shift >= 0)
        assert (varpi > shift0 and varpi > shift)
        w_grid = np.concatenate((w_dist[0] + 1j*shift0, w_dist[1:] + 1j*shift, w_dist + 1j*varpi))
        return w_grid


def frequency_distribution(npol, wrange, alpha=1):

    if npol == 1:
        w_grid = 0
        return w_grid

    assert(wrange[0] >= 0)
    assert(wrange[1] > wrange[0])

    if alpha == 0:  # homogeneous distribution
        w_grid = np.linspace(wrange[0], wrange[1], 2*npol)
        return w_grid

    ws = wrange[1] - wrange[0]
    partition = np.ones(npol)
    for i in range(1, npol+1):
        partition[i] = pivot_slice(i, 1/3)  # semi-homogeneous partition
    wgrid = wrange[0] + ws*partition**alpha
    

def pivot_slice(npol, pivot):

    assert(0 < pivot < 1)

    pivot_slice(1, pivot) = 0
    pivot_slice(2, pivot) = 1
    pivot_slice(3, pivot) = pivot
    if(npol % 2 == 0):
        pivot_slice(npol, pivot) = (pivot_slice(npol-1, pivot) + pivot_slice(npol-3, pivot))/2
    else:
        pivot_slice(npol, pivot) = (pivot_slice(npol-2, pivot) + pivot_slice(npol-3, pivot))/2
    pivot_slice(2, pivot) = 1


    assert(npol > 1)
    lp = int(np.log(npol - 1) / np.log(2))
    r = int((npol - 1) % (2**lp))
       


def mpa_frequency_sampling(npol, w0, d, ps='2l', alpha=1):
    """
    This function creates a frequency grid in the complex plane.
    The grid can have 1 or 2 branches with points non homogeneously
    distributed along the real frequency axis.
    See Fig. 1 and Eq. (18) of Ref. [1], and Eq. (11) of Ref. [2]
    [1] DA. Leon et al, PRB 104, 115157 (2021)
    [2] DA. Leon et al, PRB 107, 155130 (2023)

    Parameters
    ----------
    npol : numper of poles (half the number of frequency points)
    w0 : array of two complex numbers defining the sampling range
    d : array of two real numbers defining the damping range
    ps : string of length 2 defining a sampling with 1 or 2 lines
    alpha : exponent of the distribution of points along the real axis
    ______________________________________________________________________
                  Example: double parallel sampling with 9 poles
    ----------------------------------------------------------------------
                            complex frequency plane
    ----------------------------------------------------------------------
    |(w0[0].real, w0[0].imag) .. . . . . . .   . (w0[1].real, w0[1].imag)|
    |                                                                    |
    |     (w0[0].real, d[0])  .. . . . . . .   . (w0[1].real, d[1])      |
    ______________________________________________________________________
    """
 
    w0 = np.array(w0)
    assert (w0.real >= 0).all()
    assert (w0.imag >= 0).all()
    assert w0[1].real >= w0[0].real # positive interval

    if npol == 1:
        w_grid = np.array(w0, dtype=complex)
        return w_grid

    if ps == '1l':  # DALV: We could use a match-case function
        if alpha == 0:
            return np.linspace(w0[0], w0[1], 2 * npol)
        raise ValueError("If ps = '1l', only alpha = 0 is implemented")

    if ps == '2l':
        assert (w0.imag > d).all() # two branches
        if alpha == 0: # homogeneous distribution
            w_grid = np.concatenate((np.linspace(complex(np.real(w0[0]),
                                     d[1]), complex(np.real(w0[1]), d[1]),
                                    npol), np.linspace(w0[0], w0[1],
                                                       npol)))
            w_grid[0] = complex(np.real(w0[0]), d[0])
            return w_grid
        # inhomogeneous distribution in powers of 1/2^alpha
        ws = w0[1] - w0[0] # complex length
        w_grid = np.ones(2 * npol, dtype=complex) # number of freqs = 2 npol
        w_grid[0] = complex(np.real(w0[0]), d[0])
        w_grid[npol - 1] = complex(np.real(w0[1]), d[1])
        w_grid[npol] = w0[0]
        w_grid[2 * npol - 1] = w0[1]
        lp = int(np.log(npol - 1) / np.log(2))
        r = int((npol - 1) % (2**lp))

        if r > 0:
            for i in range(1, 2 * r):
                w_grid[npol + i] = w0[0] + ws * (i / 2.**(lp + 1)
                                                 )**alpha
                w_grid[i] = complex(np.real(w_grid[npol + i]), d[1])
            for i in range(2 * r, npol):
                w_grid[npol + i] = w0[0] + ws * ((i - r) / 2.**(lp)
                                                 )**alpha
                w_grid[i] = complex(np.real(w_grid[npol + i]), d[1])
        else:
            w_grid[npol + 1] = w0[0] + ws / (2.**(lp + 1))**alpha
            w_grid[1] = complex(np.real(w_grid[npol + 1]), d[1])
            for i in range(2 * r + 2, npol - 1):
                w_grid[npol + i] = w0[0] + ws * ((i - 1 - r) / 2.**(lp)
                                                 )**alpha
                w_grid[i] = complex(np.real(w_grid[npol + i]), d[1])
        return w_grid
    raise ValueError("Only '1l' or '2l' values are implemented")
