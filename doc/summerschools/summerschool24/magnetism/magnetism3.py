# %%
"""
# Discover new ferromagnetic monolayers with high critical temperatures

Now that you know how to calculate exchange coupling constants and single-ion anisotropies, you are in a position to search for new ferromagnetic monolayers yourself.

First, we will try to find some suitable candidates based on the data, which is already in the C2DB.

1.   Go to the C2DB website at https://cmrdb.fysik.dtu.dk/c2db and search for magnetic monolayers with finite band gaps (the theory developed so far is not applicable to metals, since these typically have long range interactions). (Hint: Toggle `Magnetic` to `yes` and set a minimum value for the `Band gap range` in the search menu, then press the search icon.)
2.   Refine the query to sort out the magnetic monolayers with a ferromagnetic ground state. (Hint: Key values such as the nearest neighbor exchange coupling can be queried by entering `J>0` into the search field.)
3.   Show the spin state, exchange coupling and magnetic anisotropy in the property overview. (Hint: You can add key values to the overview through the `Add Column` botton.)
4.   Decide on a material you find promising to have a high Curie temperature, based on what you learned in previous notebooks. You may also want to choose a material with a limited number of atoms in the unit cell in order for the calculations to run fast.

Secondly, you will write your own code below to

5.   Calculate the energy difference per magnetic site of a ferromagnetic and an antiferromagnetic configuration.
6.   Calculate $J$.
7.   Calculate the single-ion anisotropy $A$.
8.   Estimate the Curie temperature $T_c$.

You are welcome to repeat the process for as many monolayers as you like.
"""

# %%
