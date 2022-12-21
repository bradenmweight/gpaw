from gpaw.xc.fxc import FXCCorrelation
from ase.parallel import paropen

fxc = FXCCorrelation('diam_kern.ralda.lda_wfcs.gpw', xc='JGMs',
                     ecut=[131.072],
                     txt='diam_kern.ralda_04_jgm.txt',
                     avg_scheme='wavevector',
                     Eg=7.3)
E_i = fxc.calculate()

resultfile = paropen('diam_kern.ralda_kernel_comparison.dat', 'a')
resultfile.write(str(E_i[-1]) + '\n')
resultfile.close()
