from gpaw.xc.fxc import FXCCorrelation

fxc = FXCCorrelation('diam_kern.ralda.lda_wfcs.gpw', xc='rALDA',
                     txt='diam_kern.ralda_02_ralda_dens.txt',
                     ecut=[131.072])
E_i = fxc.calculate()
