import numpy as np
from pylab import *

A = np.loadtxt('rpa_N2.dat').transpose()
plot(A[0]**(-1.5), A[1], 'o', label='Calculated points')

xs = np.array([A[0,0]+i*100000. for i in range(50000)])
plot(xs**(-1.5), -4.969+1993*xs**(-1.5), label='-4.969+1993*E^(-1.5)')

t = [int(A[0,i]) for i in range(len(A[0]))]
xticks(A[0]**(-1.5), t, fontsize=12)
axis([0.,150**(-1.5), None, -4.])
xlabel('Cutoff energy [eV]', fontsize=18)
ylabel('RPA correlation energy [eV]', fontsize=18)
legend(loc='lower right')
show()
savefig('extrapolate.png')
