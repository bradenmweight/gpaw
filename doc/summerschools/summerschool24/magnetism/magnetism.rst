.. _magnetism:

===============
Magnetism in 2D
===============

This exercise investigates magnetic order in two dimensions. While a collinear
spin-density functional theory calculation might reveal that it is
energetically favorable for a given 2D material to order magnetically, two
dimensional magnetic order at finite temperatures also requires presence of
magnetic anisotropy, usually arising from the spin-orbit coupling.

This exercise will teach you how to extract magnetic exchange and anisotropy
parameters for a localized spin model based on first principles calculations.
It will also touch upon the Mermin-Wagner theorem and show why anisotropy is
crucial in order to sustain magnetic order in two dimensions.

In the first part of the project, you will calculate the Curie temperature of
a |CrI3| monolayer. Afterwards, you will investigate the magnetic order in
|VI2|, which has antiferromagnetic coupling and noncollinear order. Finally,
you will finish the project by performing a search for new magnetic 2D materials
with high critical temperatures based on a database of hypothetical monolayers.


Part 1: Critical temperature of |CrI3|
======================================

Following the instructions in the Jupyter notebook
:download:`magnetism1.ipynb`
you will in this first part of the project set up and relax a monolayer of |CrI3|
after which you will calculate its critical temperature using a prerelaxed
structure file (:download:`CrI3.xyz`).

The procedure will be as follows:

1) Set up the atomic structure and optimize the geometry of |CrI3|
2) Calculate the nearest neighbor Heisenberg exchange coupling based on a total
   energy mapping analysis
3) Show that the magnetic ground state is thermodynamically unstable when
   anisotropy is neglected (The Mermin-Wagner theorem)
4) Calculate the single-ion magnetic anisotropy and estimate the critical
   temperature


Part 2: Noncollinear magnetism in |VI2|
=======================================

In materials where the dominant magnetic exchange coupling is antiferromagnetic,
or in cases where different exchange couplings compete, the ground state may
have a complicated noncollinear magnetic order. Completing the notebook
:download:`magnetism2.ipynb`,
you will examine a prototypical monolayer with a noncollinear ground state,
namely |VI2|. Starting from the structure file
:download:`VI2.xyz`,
you will:

1) Relax the atomic structure using LDA
2) Compare a collinear antiferromagnetic structure with the ferromagnetic state
3) Obtain the noncollinear ground state
4) Calculate the magnetic anisotropy and discuss whether or not the material
   will exhibit magnetic order at low temperatures


Part 3: Find new ferromagnetic monolayers with high critical temperatures
=========================================================================

In this last part of the project, you will try to find new ferromagnetic
monolayers that can preserve their magnetic ordering at elevated temperatures.
Using the notebook
:download:`magnetism3.ipynb`,
you will:

1) Search through a database of monolayers to pick a material you might expect
   to have a high critical temperature
2) Carry out a total energy mapping analysis to obtain exchange coupling and
   anisotropy parameters
3) Calculate a first principles estimate of the critical temperature

You are welcome to repeat this procedure for as many monolayers as you like.

.. |CrI3| replace:: CrI\ :sub:`3`

.. |VI2| replace:: VI\ :sub:`2`
