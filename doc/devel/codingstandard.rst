.. _codingstandard:

==================
Coding Conventions
==================

Python Coding Conventions
=========================

Follow ASE's :ref:`ase:coding conventions`.


C-code
======

Code C in the C99 style::

  for (int i = 0; i < 3; i++) {
      double f = 0.5;
      a[i] = 0.0;
      b[i + 1] = f * i;
  }

and try to follow PEP7_.

Use **M-x c++-mode** in emacs.

.. _PEP7: https://www.python.org/dev/peps/pep-0007
