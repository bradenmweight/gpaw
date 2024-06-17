from __future__ import annotations
from typing import TYPE_CHECKING

from ase.parallel import parprint
from ase.utils.timing import Timer
from pathlib import Path
import numpy as np

from gpaw.new.ase_interface import ASECalculator
from gpaw.nlopt.basic import NLOData
from gpaw.utilities.progressbar import ProgressBar

if TYPE_CHECKING:
    from gpaw.nlopt.adapters import CollinearGSInfo, NoncollinearGSInfo
    from gpaw.typing import ArrayND


def get_mml(gs: CollinearGSInfo | NoncollinearGSInfo,
            spin: int,
            ni: int,
            nf: int,
            timer: Timer | None = None) -> ArrayND:
    """
    Compute momentum matrix elements.

    Parameters
    ----------
    gs
        Ground state adapter.
    spin
        Spin channel index (for spin-polarized systems 0 or 1).
    ni, nf
        First and last band to compute the mml (0 to nb).
    timer
        Timer to keep track of time.

    Returns
    -------
    p_kvnn
        Momentum matrix elements in atomic units gathered on master.
    """

    # Start the timer
    if timer is None:
        timer = Timer()
    parprint(f'Calculating momentum matrix elements for spin channel {spin}.')

    # Specify desired range and number of bands in calculation
    bands = slice(ni, nf)
    nb = nf - ni

    # Spin input
    assert spin < gs.ns, 'Wrong spin input'

    # Parallelisation and memory estimate
    ibzwfs = gs.ibzwfs
    kpt_comm = ibzwfs.kpt_comm
    rank = kpt_comm.rank
    master = (rank == 0)

    nk = len(ibzwfs.rank_k)  # Total number of k-points
    k_q = np.array(list(ibzwfs.q_k.keys()), int)  # k-index for each q-index
    nq = len(k_q)  # Number of k-points (q-indices) for each core
    est_mem = 2 * 3 * nk * nb**2 * 16 / 2**20
    parprint(f'At least {est_mem:.2f} MB of memory is required on master.')

    # Allocate the matrix elements
    p_qvnn = np.empty((nq, 3, nb, nb), dtype=complex)

    # Initial call to print 0 % progress
    if master:
        pb = ProgressBar()

    # Calculate matrix elements in loop over k-points
    for wfs_s in ibzwfs.wfs_qs:
        wfs = gs.get_wfs(wfs_s, spin)

        with timer('Contribution from pseudo wave functions'):
            G_plus_k_Gv, u_nG = gs.get_plane_wave_coefficients(
                wfs, bands=bands, spin=spin)
            p_vnn = np.einsum('Gv,mG,nG->vmn',
                              G_plus_k_Gv, u_nG.conj(), u_nG) * gs.ucvol

        with timer('Contribution from PAW corrections'):
            P_ani = gs.get_wave_function_projections(
                wfs, bands=bands, spin=spin)
            for P_ni, nabla_iiv in zip(P_ani.values(), gs.nabla_aiiv):
                p_vnn -= 1j * np.einsum('mi,nj,ijv->vmn',
                                        P_ni.conj(), P_ni, nabla_iiv)

        p_qvnn[wfs.q] = p_vnn

        if master:
            pb.update(wfs.q / nq)

    if master:
        pb.finish()

    with timer('Gather the data to master'):
        if not master:
            kpt_comm.send(np.array(nq, int), 0)
            kpt_comm.send(k_q, 0)
            kpt_comm.send(p_qvnn, 0)
        else:
            p_kvnn = np.empty((nk, 3, nb, nb), complex)
            p_kvnn[k_q] = p_qvnn
            for gather_rank in range(1, kpt_comm.size):
                _nq = np.empty(1, int)  # We can only communicate numpy arrays
                kpt_comm.receive(_nq, gather_rank)
                nq = _nq[0]

                k_q = np.empty(nq, int)
                kpt_comm.receive(k_q, gather_rank)

                p_qvnn = np.empty((nq, 3, nb, nb), complex)
                kpt_comm.receive(p_qvnn, gather_rank)
                p_kvnn[k_q] = p_qvnn

    # Print the timing
    if master:
        timer.write()

    if rank == 0:
        return p_kvnn
    else:
        return np.array([], dtype=complex)


def make_nlodata(calc: ASECalculator | str | Path,
                 spin_string: str = 'all',
                 ni: int | None = None,
                 nf: int | None = None) -> NLOData:
    """
    This function calculates and returns all required
    NLO data: w_sk, f_skn, E_skn, p_skvnn.

    Parameters
    ----------
    calc
        Calculator or string/path pointing to a .gpw file.
    spin_string
        String denoting which spin channels to include ('all', 's0' , 's1').
    ni
        First band to compute the mml.
    nf
        Last band to compute the mml (relative to number of bands for nf <= 0).

    Returns
    -------
    NLOData
        Data object carrying required matrix elements for NLO calculations.

    """

    if not isinstance(calc, ASECalculator):
        if not (isinstance(calc, str) or isinstance(calc, Path)):
            raise TypeError('Input must be a calculator or a string / path'
                            'pointing to a calculator.')
        from gpaw.new.ase_interface import GPAW
        calc = GPAW(calc, txt=None, parallel={'domain': 1, 'band': 1})
    assert not calc.symmetry.point_group, \
        'Point group symmetry should be off.'

    gs: CollinearGSInfo | NoncollinearGSInfo
    if calc.dft.state.density.collinear:
        from gpaw.nlopt.adapters import CollinearGSInfo
        gs = CollinearGSInfo(calc)
    else:
        from gpaw.nlopt.adapters import NoncollinearGSInfo
        gs = NoncollinearGSInfo(calc)

    # Parse spin string
    ns = gs.ns
    if spin_string == 'all':
        spins = list(range(ns))
    elif spin_string == 's0':
        spins = [0]
    elif spin_string == 's1':
        spins = [1]
        assert spins[0] < ns, 'Wrong spin input'
    else:
        raise NotImplementedError

    # Parse band input
    ibzwfs = gs.ibzwfs
    nb_full = ibzwfs.nbands
    ni = int(ni) if ni is not None else 0
    nf = int(nf) if nf is not None else nb_full
    nf = nb_full + nf if (nf <= 0) else nf

    # Start the timer
    timer = Timer()

    # Get the energy and Fermi-Dirac occupations (data is only in master)
    with timer('Get energies and fermi levels'):

        E_skn, f_skn = ibzwfs.get_all_eigs_and_occs()

        w_sk = np.array([ibzwfs.ibz.weight_k for _ in range(gs.ndensities)])
        w_sk *= gs.bzvol * ibzwfs.spin_degeneracy

    # Compute the momentum matrix elements
    with timer('Compute the momentum matrix elements'):
        p_skvnn = []
        for spin in spins:
            p_kvnn = get_mml(gs=gs, ni=ni, nf=nf,
                             spin=spin, timer=timer)
            p_skvnn.append(p_kvnn)
        if not gs.collinear:
            p_skvnn = [p_skvnn[0] + p_skvnn[1]]

    # Save the output to the file
    return NLOData(w_sk=w_sk,
                   f_skn=f_skn[:, :, ni:nf],
                   E_skn=E_skn[:, :, ni:nf],
                   p_skvnn=np.array(p_skvnn, complex),
                   comm=ibzwfs.kpt_comm)


def get_rml(E_n, p_vnn, pol_v, Etol=1e-6):
    """
    Compute the position matrix elements

    Parameters
    ----------
    E_n
        Band energies.
    p_vnn
        Momentum matrix elements.
    pol_v
        Tensor element.
    Etol
        Tolerance in energy to consider degeneracy.

    Returns
    -------
    r_vnn
        Position matrix elements.
    D_vnn
        Velocity difference matrix elements.

    """

    # Useful variables
    nb = len(E_n)
    r_vnn = np.zeros((3, nb, nb), complex)
    D_vnn = np.zeros((3, nb, nb), complex)
    E_nn = np.tile(E_n[:, None], (1, nb)) - \
        np.tile(E_n[None, :], (nb, 1))
    zeroind = np.abs(E_nn) < Etol
    E_nn[zeroind] = 1
    # Loop over components
    for v1 in set(pol_v):
        r_vnn[v1] = p_vnn[v1] / (1j * E_nn)
        r_vnn[v1, zeroind] = 0
        p_n = np.diag(p_vnn[v1])
        D_vnn[v1] = np.tile(p_n[:, None], (1, nb)) - \
            np.tile(p_n[None, :], (nb, 1))

    return r_vnn, D_vnn


def get_derivative(E_n, r_vnn, D_vnn, pol_v, Etol=1e-6):
    """
    Compute the generalised derivative of position matrix elements

    Parameters
    ----------
    E_n
        Band energies.
    r_vnn
        Momentum matrix elements.
    D_vnn
        Velocity difference matrix elements.
    pol_v
        Tensor element.
    Etol
        Tolerance in energy to consider degeneracy.

    Returns
    -------
    rd_vvnn
        Generalised derivative of position matrix elements.

    """

    # Useful variables
    nb = len(E_n)
    rd_vvnn = np.zeros((3, 3, nb, nb), complex)
    E_nn = np.tile(E_n[:, None], (1, nb)) - \
        np.tile(E_n[None, :], (nb, 1))
    zeroind = np.abs(E_nn) < Etol
    E_nn[zeroind] = 1
    for v1 in set(pol_v):
        for v2 in set(pol_v):
            tmp = (r_vnn[v1] * np.transpose(D_vnn[v2])
                   + r_vnn[v2] * np.transpose(D_vnn[v1])
                   + 1j * np.dot(r_vnn[v1], r_vnn[v2] * E_nn)
                   - 1j * np.dot(r_vnn[v2] * E_nn, r_vnn[v1])) / E_nn
            tmp[zeroind] = 0
            rd_vvnn[v1, v2] = tmp

    return rd_vvnn
