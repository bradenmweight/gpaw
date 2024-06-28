# %%
"""
# Discussion
"""

# %%
"""
Now, we will finally visualize the absorption spectra we have calculated. The script below is setup for plotting and saving the absorption spectra for BN with 24x24x1 k-points. First we load the x, y, and z components and plot all three curves in different colors.
"""

# %%
#Plot the absorption spectrum:
import numpy as np
import matplotlib.pyplot as plt
plt.figure()

#Of course you need to input the name you gave the files earlier own yourself and only plot the components you calculated.
#Here is only showed hot to plot the x-component.

absox = np.loadtxt('CdTe_rpa_x.csv', delimiter=',') # student: absox = np.loadtxt('???_rpa_x.csv', delimiter=',')
plt.plot(absox[:, 0], absox[:, 4], label='RPA_x CdTe 12x12x4', lw=2, color='b')

plt.xlabel(r'$\hbar\omega\;[eV]$', size=20)
plt.ylabel(r'$Im(\epsilon)$', size=20)
plt.xticks(size=16)
plt.yticks(size=16)
plt.tight_layout()
plt.axis([0.0, 10.0, None, None])
plt.legend()

plt.show()


# %%
"""
Now open the saved figure and inspect the absorption spectra. Does it look as you expected? Can you guess the band gap from the absorption spectra and which of the band gaps calculated on day 2 does it match?

Now plot the absorption spectra for the other k-point meshes you calculated and compare. You could for instance modify the script above to plot, in different colors, all spectra in the same plot for comparison. Would you say that the calculation is converged?

Finally, talk with the other groups and compare your absorption spectra with their absorption spectra. Is the same number of k-points needed to obtain same degree of convergence for the different materials? Is there any difference for the 2D material (Boron-Nitride)?

For a more realistic description one would have to include excitonic effects (i.e. the electron-hole interaction). As an outlook you might consider reading:
https://wiki.fysik.dtu.dk/gpaw/documentation/bse/bse.html
"""
