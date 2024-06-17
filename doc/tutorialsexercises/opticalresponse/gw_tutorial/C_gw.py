from gpaw.response.g0w0 import G0W0

gw = G0W0(calc='C_groundstate.gpw',
          nbands=30,  # number of bands for calculation of self-energy
          bands=(3, 5),          # VB and CB
          ecut=20.0,             # plane-wave cutoff for self-energy
          integrate_gamma='WS',  # Use supercell Wigner-Seitz truncation for W.
          filename='C-g0w0')
gw.calculate()
