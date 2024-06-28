# %%
# teacher
import ase.visualize as viz
viz.view = lambda atoms, repeat=None: None

# %%
"""

# Noncollinear magnetism - VI$_2$

For this part of the project, we will move on to a magnetic monolayer in which the dominant exchange coupling is antiferromagnetic, namely VI$_2$. We will still use the localized spin Hamiltonian

$$H = -\frac{1}{2}\sum_{i,j}J_{ij}\mathbf{S}_i\cdot \mathbf{S}_j+A\sum_i(S_i^z)^2$$

but here in the antiferromagnetic case, $J<0$.

## Optimizing the atomic structure

Since we will need to do LDA calculations later on, we will start of this part of the project by relaxing the atomic structure of VI$_2$ using the LDA functional. Usually, the difference in crystal structure between different magnetically ordered states is small, so we will perform the relaxation in the ferromagnetic state, which has a smaller unit cell.

1.   First you should download the relaxed PBE crystal structure. Either, browse the C2DB at https://cmrdb.fysik.dtu.dk/c2db and download the `.xyz` file for VI$_2$ or dowload it directly from the summer school tutorial website [here](https://wiki.fysik.dtu.dk/gpaw/summerschools/summerschool22/magnetism/magnetism.html)
2.   Fill in the expected ionic value for the V spins `S` below and run the cell to relax the crystal structure. The calculation takes about 17 minutes. (Hint: V has the electronic configuration [Ar]3d$^3$4s$^2$)

"""

# %%
from ase.io import read
from ase.visualize import view
from ase.optimize import BFGS
from ase.filters import UnitCellFilter
from gpaw import GPAW, PW

S = 3 / 2  # student: S = ???
m = S * 2
layer = read('VI2.xyz')
layer.set_initial_magnetic_moments([m, 0, 0])
view(layer)

calc = GPAW(mode=PW(400), xc='LDA', kpts=(4, 4, 1))
layer.calc = calc

uf = UnitCellFilter(layer, mask=[1, 1, 0, 0, 0, 1])
opt = BFGS(uf)
opt.run(fmax=0.1)

calc.set(symmetry='off')
calc.get_potential_energy()

calc.write('VI2_relaxed.gpw')

# %%
"""
## Magnetic anisotropy

In the cell above, we not only performed a structural optimization, but we also took the relaxed structure, switched off the $k$-point symmetries and did a final DFT calculation with all the $k$-points in the Brillouin zone. We did this in order to be able to evaluate the the magnetic anisotropy arising from the spin-orbit coupling.

1.   Adapt the code you used for CrI$_3$ to calculate the magnetic anisotropy in the cell below.
2.   Is the easy-axis in-plane or out-of-plane?
3.   Do you expect to the VI$_2$ to exhibit magnetic order at finite temperatures?

"""

# %%
from gpaw.spinorbit import soc_eigenstates
# teacher:

calc = GPAW('VI2_relaxed.gpw')
e_x, e_y, e_z = (
    soc_eigenstates(calc, theta=theta, phi=phi).calculate_band_energy()
    for theta, phi in [(90, 0), (90, 90), (0, 0)])
de_zx = e_z - e_x
de_zy = e_z - e_y
print(f'dE_zx = {de_zx * 1000:1.3f} meV')
print(f'dE_zy = {de_zy * 1000:1.3f} meV')
A = (de_zx + de_zy) / 2 / S**2
print(f'A = {A * 1000:1.3f} meV')

# %%
"""
## DFT calculations in a repeated cell

To realize that VI$_2$ is in fact an antiferromagnetic material, we need to do a calculation starting from an antiferromagnetic alignment of the spins. To do so, we need more than a single V atom in the unit cell. In the cell below, the atomic structure is repeated once to obtain two V atoms in the unit cell and it is shown how do a DFT calculation for the antiferromagnetic state.

1.   Fill in the `kpts` for the repeated unit cell.
2.   Replace the `...` with code to calculate the ferromagnetic state.
3.   When you have finalzed the code for the ferromagnetic state, run the cell. The calculation takes about 7 minutes.

"""

# %%
layer_afm = layer.repeat((2, 1, 1))
layer_afm.set_initial_magnetic_moments([m, 0, 0, -m, 0, 0])
view(layer_afm)

calc = GPAW(mode=PW(400),
            xc='LDA',
            kpts=(2, 4, 1))  # student: kpts=???))
layer_afm.calc = calc
layer_afm.get_potential_energy()
calc.write('V2I4_afm.gpw')

...
# teacher:
layer_fm = layer.repeat((2, 1, 1))
layer_fm.set_initial_magnetic_moments([m, 0, 0, m, 0, 0])
calc = GPAW(mode=PW(400), xc='LDA', kpts=(2, 4, 1))
layer_fm.calc = calc
layer_fm.get_potential_energy()
calc.write('V2I4_fm.gpw')

# %%
"""
## Calculating J

Finally, we are in a position to calculate the nearest neighbour Heisenberg exchange coupling $J$. Before doing so, please compare the output of the antiferromagnetic and ferromagnetic calculations from the cell above.

1.   Which state has the lowest energy?
2.   What will the sign of $J$ be?

You should find that the antiferromagnetic state has a lower energy than the ferromagnetic one, but does this also mean that the calculated configuration is the correct magnetic ground state? Rather, it implies that the spins prefer to be antialigned, i.e. that the exchange coupling $J$ is negative.

3.   Draw the structural arrangement of the V atoms on a piece of paper. Which type of magnetic lattice do they form?
4.   Fill in the spin configuration of the magnetic lattice and convince yourself that all spins cannot be antialigned to their nearest neighbours.

The latter finding means that the antiferromagnetic system is frustrated and the antiferromagnetic configuration we have computed will not be the true ground state of the system.

Leaving the magnetic frustration aside for the moment, we will first calculate $J$, which we can still do even though that we have not found the ground state yet.

5.   Use the classical Heisenberg model with nearest neighbor interaction only (and $A=0$) to derive the energy per magnetic site of the ferromagnetic and antiferromagnetic configurations respectively.

You should obtain the following:

$$E_{\mathrm{FM}} = E_0 - 3JS^2$$

and

$$E_{\mathrm{AFM}} = E_0 + JS^2$$

where $E_0$ is some reference energy.

6.   Use these expressions to eliminate $E_0$ and express $J$ in terms of the energy difference per magnetic site of the two configurations.
7.   Write code to extract `E_fm` and `E_afm` in the cell below.
8.   Fill in the formula for `J` and evaluate the cell to calculate it.

You should get a value for $J$ around -1.4 meV.
"""

# %%
E_fm = layer_fm.get_potential_energy() / 2  # student: E_fm = ???
E_afm = layer_afm.get_potential_energy() / 2  # student: E_afm = ???
dE = E_afm - E_fm  # student:
J = dE / 4 / S**2  # student: J = ???
print(f'J = {J * 1000:1.2f} meV')

# %%
"""
## Noncollinear configuration

As it turn out, the optimal spin structure on a trigonal lattice with antiferromagntice exchange coupling is to place all spins at 120$^\circ$ angles with respect its neighbors.

1.   Draw this structure and convince yourself that it is indeed possible to put every spin at a 120$^\circ$ angle with respect to its neighbors.
2.   What is the minimal number of magnetic atoms required in the magnetic unit cell to represent such a state?
3.   Verrify that the Heisenberg model with classical spins gives a lower energy with this configuration than in the antialigned structure calculated above.

You should obtain the following energy for the 120$^\circ$ noncolinnear configuration:

$$E_{\mathrm{NC}}=E_0+\frac{3}{2}JS^2.$$

We will now check if we can verify this prediction within the LSDA. To do that we need to perform a noncollinear DFT calculation, which is done in the cell below.

4.   Read and try to understand the code to perform a noncollinear LSDA calculation in the 120$^\circ$ noncolinnear configuration.
5.   Replace the `...` with code to make a noncolinnear LSDA calculation for the ferromagnetic state as well (we will need this for a late comparison).

Run the cell and verify that the energy per magnetic atom is lower in the 120$^\circ$ noncolinnear configuration compared to both of the previous calculated states.
"""

# %%
import numpy as np
from ase.io import read
from ase.visualize import view
from gpaw import GPAW, PW, MixerDif

m = 3
cell_cv = layer.get_cell()
layer_nc = layer.repeat((3, 1, 1))
new_cell_cv = [[3 * cell_cv[0, 0] / 2, 3**0.5 * cell_cv[0, 0] / 2, 0.0],
               [3 * cell_cv[0, 0] / 2, -3**0.5 * cell_cv[0, 0] / 2, 0.0],
               [0.0, 0.0, cell_cv[2, 2]]]
layer_nc.set_cell(new_cell_cv)
view(layer_nc)

magmoms = np.zeros((len(layer_nc), 3), float)
magmoms[0] = [m, 0, 0]
magmoms[3] = [m * np.cos(2 * np.pi / 3), m * np.sin(2 * np.pi / 3), 0]
magmoms[6] = [m * np.cos(2 * np.pi / 3), -m * np.sin(2 * np.pi / 3), 0]

calc = GPAW(mode=PW(400),
            xc='LDA',
            mixer=MixerDif(),
            symmetry='off',
            experimental={'magmoms': magmoms, 'soc': False},
            parallel={'domain': 1, 'band': 1},
            kpts=(2, 2, 1))
layer_nc.calc = calc
layer_nc.get_potential_energy()
calc.write('nc_nosoc.gpw')

...
# teacher:
magmoms = np.zeros((len(layer_nc), 3), float)
magmoms[0] = [m, 0, 0]
magmoms[3] = [m, 0, 0]
magmoms[6] = [m, 0, 0]
calc = GPAW(mode=PW(400),
            xc='LDA',
            mixer=MixerDif(),
            symmetry='off',
            experimental={'magmoms': magmoms, 'soc': False},
            parallel={'domain': 1, 'band': 1},
            kpts=(2, 2, 1))
layer_nc.calc = calc
layer_nc.get_potential_energy()
calc.write('fm_nosoc.gpw')

# %%
"""
## Anisotropy and exchange coupling from noncollinear DFT

In the cell above we could have included spin-orbit coupling in a self-consistent way by setting 'soc'=True'. However, it is more convenient for us to exclude it such that we can explicitly compute the single-ion anisotropy parameter $A$ afterwards. If the localized spin Hamiltonian with nearest neighbor exchange interactions is a good model, we should be able to obtain both $J$ and $A$ from the noncollinear calculation as well.

1.   Fill in the cell below with code to compute the nearest neighbour exchange coupling $J$ based on the noncollinear calculations performed above.
2.   Fill in code that computes the spin-orbit coupling corrected energies of the ferromagnetic state from the noncollinear calculation with spins directed along the $x$, $y$ and $z$ directions (Hint: the anisotropy is calculated by rotating the entire initial spin configuration first by $\theta$ and then by $\varphi$).
3.   Repeat point 2, but for the 120$^\circ$ noncolinnear configuration. In this case, you cannot align all spins to one direction, but takes as a reference the first V atom.

With the above code in place, please evaluate the cell. You should obtain a nearest neighbour exchange coupling of about -1.9 meV.

"""

# %%

# teacher:
calc = GPAW('nc_nosoc.gpw', txt=None)
E_nc = calc.get_potential_energy() / 3
e_x, e_y, e_z = (
    soc_eigenstates(calc, theta=theta, phi=phi).calculate_band_energy() / 3
    for theta, phi in [(0, 0), (0, 90), (90, 0)])
de_zx = e_z - e_x
de_zy = e_z - e_y
print(f'NC: A_zx = {de_zx * 1000:1.3f} meV')
print(f'NC: A_zy = {de_zy * 1000:1.3f} meV')
print()

calc = GPAW('fm_nosoc.gpw', txt=None)
E_fm = calc.get_potential_energy() / 3
e_x, e_y, e_z = (
    soc_eigenstates(calc, theta=theta, phi=phi).calculate_band_energy() / 3
    for theta, phi in [(0, 0), (0, 90), (90, 0)])
de_zx = e_z - e_x
de_zy = e_z - e_y
print(f'FM: A_zx = {de_zx * 1000:1.3f} meV')
print(f'FM: A_zy = {de_zy * 1000:1.3f} meV')
print()

dE = E_nc - E_fm
J = 2 * dE / 9 / S**2
print(f'J = {J * 1000:1.2f} meV')
print()
calc = GPAW('VI2_relaxed.gpw', txt=None)
E_fm = calc.get_potential_energy()
dE = E_nc - E_fm
J = 2 * dE / 9 / S**2
print(f'J = {J * 1000:1.2f} meV')


# %%
"""
## Critical temperature?

Based on the noncollinear calculations above:

1.   What is the easy axis of the system?
2.   Does the easy axis agree with your initial findings for the simple ferromagnetic state?

You should find that the energy of the 120$^\circ$ noncolinnear configuration does not depend strongly on whether the first V atom is aligned to the $x$ or the $y$ direction.

3.   If we assume full in-plane isotropy is there any rotational freedom left in the noncollinear ground state?
4.   What implications does this have for the critical temperature of the monolayer?

You might be able to convince yourself that some degree of in-plane anisotropy is required to obtain a finite critical temperature for the 120$^\circ$ noncolinnear magnetic order. Again, bear in mind that all the calculations in the present notebook ought to be properly converged with respect to $k$-points, plane wave cutoff etc. to achieve an accurate estimate of e.g. the in-plane anisotropy.

Clearly the noncollinear spin state of VI$_2$ is more difficult to describe than the ferromagnetic state in CrI$_3$ and we do not yet have a simple theoretical expression for the critical temperature as a function of anisotropy and exchange coupling constants. However, with the rapid development of experimental techniques to synthesize and characterize 2D materials it does seem plausible that such a noncollinear 2D material may be observed in the future.
"""
