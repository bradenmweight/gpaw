.. _inducedfield:

===============================================================
Induced density oscillations and electric near field from TDDFT
===============================================================

Calculating Fourier transform of the density oscillation at
the resonant frequencies of a systems is a great way of analyzing
the excitations.
Induced density oscillations and the resulting electric near field
can be calculated for finite systems from
:ref:`inducedfield_timepropagation` or :ref:`inducedfield_casida`.

For details of the implementation, see Ref. [#Rossi2013]_.

.. _inducedfield_timepropagation:

Time-propagation TDDFT
===========================================

See :ref:`timepropagation` for instructions how to use time-propagation TDDFT.

Within time-propagation TDDFT, the induced density is obtained as an on-the-fly
Fourier transform.
A restriction of this iterative approach is that the frequencies of interest
must be given at initialization, that is, before time propagation.

Example code for time-propagation calculation

.. literalinclude:: timepropagation_calculate.py

Example code for continuing time-propagation

.. literalinclude:: timepropagation_continue.py

Induced electric potential and near field are calculated after time-propagation as follows:

.. literalinclude:: timepropagation_postprocess.py

Plotting example (run in serial mode, i.e., with one process)

.. literalinclude:: timepropagation_plot.py

|na2_td_Frho| |na2_td_Fphi| |na2_td_Ffe|

.. |na2_td_Frho| image:: na2_td_Frho.png
.. |na2_td_Fphi| image:: na2_td_Fphi.png
.. |na2_td_Ffe| image:: na2_td_Ffe.png


TODO: notes about AE/comp corrections, extending grid

.. _inducedfield_casida:

Linear response TDDFT (Casida's equation)
===========================================

See :ref:`lrtddft` for instructions how to use linear response TDDFT.

Example code for linear response calculation

.. literalinclude:: casida_calculate.py

Induced electric potential and near field are calculated as post-processing step as follows:

.. literalinclude:: casida_postprocess.py

Plotting example (same as in time propagation)

.. literalinclude:: casida_plot.py

|na2_casida_Frho| |na2_casida_Fphi| |na2_casida_Ffe|

.. |na2_casida_Frho| image:: na2_casida_Frho.png
.. |na2_casida_Fphi| image:: na2_casida_Fphi.png
.. |na2_casida_Ffe| image:: na2_casida_Ffe.png


References
==========

.. [#Rossi2013]
   T. P. Rossi,
   Simulating electric field enhancement in plasmonic nanomaterials,
   Master's thesis, Aalto University (2013).
   `<https://urn.fi/URN:NBN:fi:aalto-201309137662>`_
