from gpaw.xc.fxc import FXCCorrelation
from ase.parallel import paropen

fxc = FXCCorrelation('diam_kern.ralda.lda_wfcs.gpw', xc='rALDA',
                     ecut=[131.072],
                     txt='diam_kern.ralda_02_ralda_dens.txt',
                     avg_scheme='density')
E_i = fxc.calculate()

resultfile = paropen('diam_kern.ralda_kernel_comparison.dat', 'w')
resultfile.write(str(E_i[-1]) + '\n')
resultfile.close()
