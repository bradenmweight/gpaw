.. _writing documentation:

=====================
Writing documentation
=====================

.. highlight:: bash

We use the Sphinx_ tool to generate the GPAW documentation.

First, you should take a look at the documentation for Sphinx_ and
reStructuredText_.  Also, read carefully the
:ref:`Writing documentation for ASE <ase:writing_documentation_ase>`
page.

.. _reStructuredText: https://docutils.sf.net/rst.html
.. _Sphinx: https://www.sphinx-doc.org

**Structure**

When writing documentation easy accessibility and readability is key. To that
end the documentation is split into several parts:

:ref:`Documentation/Basic usage <basic>`: This part contains basic usage
instructions for GPAW, including references to parameters for the GPAW
calculator object. This part should not contain extended examples, theory or
code references.

:ref:`Documentation/Advanced topics <advanced>`: This part contains
explanations of the various features of GPAW. The focus here is on
implementation specific information, not theory, as well as code references.

:ref:`Documentation/Theory <theory>`: This is the place for theoretical
descriptions of methods used in the code. Reference to literature should be
given.

:ref:`Tutorials and Exercises <tutorials>`: As the name suggests, this is the
heading for any worked out examples, tutorials and exercises. Entries are
further sorted into fields of physics or application.

One should always the different pages relating to one topic for easy navigation
between theory, implementation and example sections.

**Getting started**

If you don't already have your own copy of the GPAW package, then
see :ref:`development workflow` for how to set up your invironment.
You will also need to ``pip install sphinx-rtd-theme``.

Then :command:`cd` to the :file:`doc` directory and build the html-pages::

  $ cd ~/gpaw/doc
  $ make

Make your changes to the ``.rst`` files, run the
:command:`make` command again, check the results and if things
looks ok, commit::

    $ emacs index.rst
    $ make
    $ make browse
    $ git add index.rst
    $ git commit -m "..."


**Adding figures and tables**

We don't want to have png and csv files committed to Git.  Instead, you should
add the Python scripts that generate the figures and table data so that we can
always generate them again if needed.

For quick scripts (no more than 5 seconds), see :ref:`ase:generated`.  For
more expensive scripts you can use :ref:`AGTS <agts>` for running long jobs
that create figures or table data for this web-page.  For an example, look at
the source code :git:`here <doc/tutorialsexercises/electronic/stm>` which will
produce this: :ref:`stm tutorial`.


.. automodule:: gpaw.doctools.aamath

.. autoclass:: gpaw.typing.Vector
