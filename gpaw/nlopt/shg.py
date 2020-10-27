
# Import the required modules: General
import numpy as np

# Import the required modules: GPAW/ASE
from ase.units import Bohr, _hbar, _e, _me, _eps0
from ase.utils.timing import Timer
from gpaw.mpi import world
from basic import load_data, parprint, print_progressbar
from nlobas import get_rml, calc_gender

# Compute the SHG spectrum and save it


def get_shg(
        freqs=[1.0],
        eta=0.05,
        pol='yyy',
        eshift=0.0,
        gauge='lg',
        ftol=1e-4, Etol=1e-6,
        band_n=None,
        out_name='shg.npy',
        mml_name='mml.npz'):
    """
    Calculate RPA SHG spectrum  for nonmagnetic semiconductors

    Input:
        freqs           Excitation frequency array (a numpy array or list)
        eta             Broadening, a number or an array (default 0.05 eV)
        pol             Tensor element (default 'yyy')
        gauge           Choose the gauge (lg or vg)
        Etol, ftol      Tol. in energy and fermi to consider degeneracy
        band_n          List of bands in the sum (default 0 to nb)
        out_name        Output filename (default 'shg.npy')
        mml_name        The momentum filename (default 'mml.npz')
    Output:
        shg.npy         Numpy array containing the spectrum and frequencies
    """

    # Start a timer
    timer = Timer()
    parprint('Calculating SHG spectrum (in {:d} cores).'.format(world.size))

    # Useful variables
    pol_v = ['xyz'.index(ii) for ii in pol]
    freqs = np.array(freqs)
    nw = len(freqs)
    w_lc = freqs + 1e-12 + 1j * eta  # Add small value to avoid 0
    # Use the TRS to reduce calculation time
    w_l = np.hstack((-w_lc[-1::-1], w_lc))
    nw = 2 * nw
    parprint('Calculation in the {} gauge for element {}.'.format(gauge, pol))

    # Load the required data
    with timer('Load and distribute the data'):
        k_info = load_data(mml_name=mml_name)
        k_ind, tmp = k_info.popitem()
        nb = len(tmp[1])
        nk = len(k_info) * world.size  # Approximately
        if band_n is None:
            band_n = list(range(nb))
        est_mem = 6 * 3 * nk * nb**2 * 16 / 2**20
        parprint('At least {:.2f} MB of memory is required.'.format(est_mem))
        # print(k_info.keys())

    # Initial call to print 0% progress
    count = 0
    ncount = len(k_info)
    print_progressbar(count, ncount)

    # Initialize the outputs
    sum2_l = np.zeros((nw), complex)
    sum3_l = np.zeros((nw), complex)

    # Do the calculations
    for k_ind, (we, f_n, E_n, p_vnn) in k_info.items():
        # Which gauge
        if gauge == 'vg':
            with timer('Sum over bands'):
                tmp = calc_shg_rvg(
                    w_l, f_n, E_n, p_vnn, band_n, pol_v,
                    ftol=ftol, Etol=Etol, eshift=eshift)
        elif gauge == 'lg':
            with timer('Position matrix elements calculation'):
                r_vnn, D_vnn = get_rml(E_n, p_vnn, pol_v, Etol=Etol)

            with timer('Compute generalized derivative'):
                rd_vvnn = calc_gender(E_n, r_vnn, D_vnn, pol_v, Etol=Etol)

            with timer('Sum over bands'):
                tmp = calc_shg_rlg(
                    w_l, f_n, E_n, r_vnn, rd_vvnn, D_vnn, band_n, pol_v,
                    ftol=ftol, Etol=Etol, eshift=eshift)
        else:
            parprint('Gauge ' + gauge + ' not implemented.')
            raise NotImplementedError

        # Add it to previous with a weight
        sum2_l += tmp[0] * we
        sum3_l += tmp[1] * we

        # Print the progress
        count += 1
        print_progressbar(count, ncount)

    with timer('Gather data from cores'):
        world.sum(sum2_l)
        world.sum(sum3_l)

    # Make the output in SI unit (2 is for spin)
    chi_l = 2 * make_output(gauge, sum2_l, sum3_l)

    # Save it to the file
    if world.rank == 0:
        # A multi-col output
        nw = len(freqs)
        chi_l = chi_l[nw:] + chi_l[nw - 1::-1]
        shg = np.vstack((freqs, chi_l))
        np.save(out_name, shg)

        # Print the timing
        timer.write()

# Implement the velocity gauge equation


def calc_shg_rvg(
        w_l, f_n, E_n, p_vnn, band_n, pol_v,
        ftol=1e-6, Etol=1e-9, eshift=0):
    """
    Loop over bands for computing in velocity gauge

    Input:
        w_l             Complex frequency array
        f_n             Fermi levels
        E_n             Energies
        p_vnn           Momentum matrix elements
        pol_v           Tensor element
        band_n          Band list
        Etol, ftol      Tol. in energy and fermi to consider degeneracy
        eshift          Bandgap correction
    Output:
        sum2_l, sum3_l  Output 2 and 3 bands terms
    """

    # Initialize variables
    sum2_l = np.zeros(w_l.size, complex)
    sum3_l = np.zeros(w_l.size, complex)

    # Loop over bands
    for nni in band_n:
        for mmi in band_n:
            # Remove non important term using TRS
            if mmi <= nni:
                continue

            # Useful variables
            fnm = f_n[nni] - f_n[mmi]
            Emn = E_n[mmi] - E_n[nni] + fnm * eshift

            # Comute the 2-band term
            if np.abs(Emn) > Etol and np.abs(fnm) > ftol:
                pnml = (p_vnn[pol_v[0], nni, mmi]
                        * (p_vnn[pol_v[1], mmi, nni]
                            * (p_vnn[pol_v[2], mmi, mmi]
                                - p_vnn[pol_v[2], nni, nni])
                            + p_vnn[pol_v[2], mmi, nni]
                            * (p_vnn[pol_v[1], mmi, mmi]
                                - p_vnn[pol_v[1], nni, nni])) / 2)
                sum2_l += 1j * fnm * np.imag(pnml) * \
                    (1 / (Emn**4 * (w_l - Emn)) -
                        16 / (Emn**4 * (2 * w_l - Emn)))

            # Loop over the last band index
            for lli in band_n:
                fnl = f_n[nni] - f_n[lli]
                fml = f_n[mmi] - f_n[lli]

                # Do not do zero calculations
                if np.abs(fnl) < ftol and np.abs(fml) < ftol:
                    continue

                # Compute the susceptibility with 1/w form
                Eln = E_n[lli] - E_n[nni] + fnl * eshift
                Eml = E_n[mmi] - E_n[lli] - fml * eshift
                pnml = (p_vnn[pol_v[0], nni, mmi]
                        * (p_vnn[pol_v[1], mmi, lli]
                            * p_vnn[pol_v[2], lli, nni]
                            + p_vnn[pol_v[2], mmi, lli]
                            * p_vnn[pol_v[1], lli, nni]))
                pnml = 1j * np.imag(pnml) / 2

                # Compute the divergence-free terms
                if np.abs(Emn) > Etol and np.abs(
                        Eml) > Etol and np.abs(Eln) > Etol:
                    ftermD = (16 / (Emn**3 * (2 * w_l - Emn))
                              * (fnl / (Emn - 2 * Eln)
                                 + fml / (Emn - 2 * Eml))) \
                        + fnl / (Eln**3 * (2 * Eln - Emn)
                                 * (w_l - Eln)) \
                        + fml / (Eml**3 * (2 * Eml - Emn)
                                 * (w_l - Eml))
                    sum3_l += pnml * ftermD

    # Return outputs
    return sum2_l, sum3_l

# Implement the length gauge equation


def calc_shg_rlg(
        w_l, f_n, E_n, r_vnn, rd_vvnn, D_vnn, band_n, pol_v,
        ftol=1e-6, Etol=1e-9, eshift=0):
    """
    Loop over bands for computing in length gauge

    Input:
        w_l             Complex frequency array
        f_n             Fermi levels
        E_n             Energies
        r_vnn           Momentum matrix elements
        rd_vvnn         Generalized derivative of position
        D_vnn           Velocity difference
        pol_v           Tensor element
        band_n          Band list
        Etol, ftol      Tol. in energy and fermi to consider degeneracy
        eshift          Bandgap correction
    Output:
        sum2_l, sum3_l  Output 2 and 3 bands terms
    """

    # Initialize variables
    sum2_l = np.zeros(w_l.size, complex)
    sum3_l = np.zeros(w_l.size, complex)

    # Loop over bands
    for nni in band_n:
        for mmi in band_n:
            # Remove the non important term using TRS
            if mmi <= nni:
                continue
            fnm = f_n[nni] - f_n[mmi]
            Emn = E_n[mmi] - E_n[nni] + fnm * eshift

            # Two band part
            if np.abs(fnm) > ftol:
                tmp = 2 * np.imag(
                    r_vnn[pol_v[0], nni, mmi]
                    * (rd_vvnn[pol_v[1], pol_v[2], mmi, nni]
                        + rd_vvnn[pol_v[2], pol_v[1], mmi, nni])) \
                    / (Emn * (2 * w_l - Emn))
                tmp += np.imag(
                    r_vnn[pol_v[1], mmi, nni]
                    * rd_vvnn[pol_v[2], pol_v[0], nni, mmi]
                    + r_vnn[pol_v[2], mmi, nni]
                    * rd_vvnn[pol_v[1], pol_v[0], nni, mmi]) \
                    / (Emn * (w_l - Emn))
                tmp += np.imag(
                    r_vnn[pol_v[0], nni, mmi]
                    * (r_vnn[pol_v[1], mmi, nni]
                        * D_vnn[pol_v[2], mmi, nni]
                        + r_vnn[pol_v[2], mmi, nni]
                        * D_vnn[pol_v[1], mmi, nni])) \
                    * (1 / (w_l - Emn)
                        - 4 / (2 * w_l - Emn)) / Emn**2
                tmp -= np.imag(
                    r_vnn[pol_v[1], mmi, nni]
                    * rd_vvnn[pol_v[0], pol_v[2], nni, mmi]
                    + r_vnn[pol_v[2], mmi, nni]
                    * rd_vvnn[pol_v[0], pol_v[1], nni, mmi]) \
                    / (2 * Emn * (w_l - Emn))
                sum2_l += 1j * fnm * tmp / 2  # 1j imag

            # Three band term
            for lli in band_n:
                fnl = f_n[nni] - f_n[lli]
                fml = f_n[mmi] - f_n[lli]
                Eml = E_n[mmi] - E_n[lli] - fml * eshift
                Eln = E_n[lli] - E_n[nni] + fnl * eshift
                # Do not do zero calculations
                if (np.abs(fnm) < ftol and np.abs(fnl) < ftol
                        and np.abs(fml) < ftol):
                    continue
                if np.abs(Eln - Eml) < Etol:
                    continue

                rnml = np.real(
                    r_vnn[pol_v[0], nni, mmi]
                    * (r_vnn[pol_v[1], mmi, lli]
                        * r_vnn[pol_v[2], lli, nni]
                        + r_vnn[pol_v[2], mmi, lli]
                        * r_vnn[pol_v[1], lli, nni])) / (2 * (Eln - Eml))
                if np.abs(fnm) > ftol:
                    sum3_l += 2 * fnm / (2 * w_l - Emn) * rnml
                if np.abs(fnl) > ftol:
                    sum3_l += -fnl / (w_l - Eln) * rnml
                if np.abs(fml) > ftol:
                    sum3_l += fml / (w_l - Eml) * rnml

    # Return outputs
    return sum2_l, sum3_l


# Make the output in SI unit

def make_output(gauge, sum2_l, sum3_l):
    """
    Make the output in SI unit and return chi

    Input:
        gauge       Chosen gauge
        sum2_l      2-bands term
        sum3_l      3-bands term
    Output:
        chi_l       Output chi as an array
    """
    # Make the output in SI unit
    if gauge == 'lg':
        dim_ee = _e**3 / (_eps0 * (2.0 * np.pi)**3)
        dim_sum = (_hbar / (Bohr * 1e-10))**3 / \
            (_e**5 * (Bohr * 1e-10)**3) * (_hbar / _me)**3
        dim_SI = dim_ee * dim_sum
        chi_l = dim_SI * (1j * sum2_l + sum3_l)
    elif gauge == 'vg':
        # Make the output in SI unit
        dim_vg = _e**3 * _hbar**2 / (_me**3 * (2.0 * np.pi)**3)
        dim_chi = 1j * _hbar / (_eps0 * 2.0 * _e)  # 2 beacuse of frequecny
        dim_sum = (_hbar / (Bohr * 1e-10))**3 / \
            (_e**4 * (Bohr * 1e-10)**3)
        dim_SI = dim_chi * dim_vg * dim_sum
        chi_l = dim_SI * (sum2_l + sum3_l)
    else:
        parprint('Gauge ' + gauge + ' not implemented.')
        raise NotImplementedError

    # Return output
    return chi_l
