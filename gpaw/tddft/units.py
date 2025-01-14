# https://www.answers.com/topic/hartree-units

from math import pi

from ase.units import _hbar, _eps0, _me, _e, _c, Bohr, Hartree, s

_a0 = Bohr * 1e-10
_autime = _hbar**3 * (4 * pi * _eps0)**2 / (_me * _e**4
                                            )  # 1 autime ~ 2.42e-17 s

# Conversion factors between ASE and GPAW units

attosec_to_autime = 1e-18 / _autime  # 1 as ~ 0.0413 autime
autime_to_attosec = _autime * 1e18  # 1 autime ~ 24.2 as

# ASE time unit is [Å sqrt(u/eV)]
# Hartree time unit is [hbar/Hartree]
asetime_to_autime = 1 / s / _autime  # 1 ASE time unit ~ 420.9 autime
autime_to_asetime = 1 / asetime_to_autime

eV_to_hartree = 1.0 / Hartree  # 1 eV ~ 0.0368 Eh
hartree_to_eV = Hartree  # 1 Eh ~ 27.2 eV

hartree_to_aufrequency = 1.0  # Eh ~ autime^(-1) since hbar=1
aufrequency_to_hartree = 1.0  # autime^(-1) ~ Eh since hbar=1

eV_to_aufrequency = 1.0 / Hartree  # 1 eV ~ 0.0368 autime^(-1)
aufrequency_to_eV = Hartree  # 1 autime^(-1) ~ 27.2 eV

eA_to_audipolemoment = 1 / Bohr  # 1 eÅ ~ 1.8897 aulength^(-1)
audipolemoment_to_eA = 1 / eA_to_audipolemoment

assert eV_to_aufrequency == eV_to_hartree * hartree_to_aufrequency
assert aufrequency_to_eV == aufrequency_to_hartree * hartree_to_eV

# Short-hand names
eV_to_au = eV_to_aufrequency
au_to_eV = aufrequency_to_eV
as_to_au = attosec_to_autime
au_to_as = autime_to_attosec
au_to_fs = au_to_as / 1e3
fs_to_au = as_to_au * 1e3
eA_to_au = eA_to_audipolemoment
au_to_eA = audipolemoment_to_eA

# Rotatory strength
# See https://doi.org/10.1016/0009-2614(95)01036-9
rot_au_to_SI = _e**2 * _hbar / _me * _a0  # 1 au = 1.6e-52 J C m / T
rot_au_to_cgs = rot_au_to_SI / (1e-6 / _c)  # 1 au = 4.7e-38 erg esu cm / gauss
