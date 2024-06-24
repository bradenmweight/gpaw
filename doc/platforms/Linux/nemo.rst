====
nemo
====

Information about `nemo <https://www.hpc.uni-freiburg.de/nemo>`__.

Building GPAW
=============

We assume that the installation will be located in ``$HOME/source``.

Setups
------

The setups of your choice must be installed
(see also :ref:`installation of paw datasets`)::

  cd
  GPAW_SETUP_SOURCE=$PWD/source/gpaw-setups
  mkdir -p $GPAW_SETUP_SOURCE
  cd $GPAW_SETUP_SOURCE
  wget https://wiki.fysik.dtu.dk/gpaw-files/gpaw-setups-0.9.11271.tar.gz
  tar xzf gpaw-setups-0.9.11271.tar.gz

Let gpaw know about the setups::

  export GPAW_SETUP_PATH=$GPAW_SETUP_SOURCE/gpaw-setups-0.9.11271

Using the module environment
----------------------------

It is very handy to add our installation to the module environment::

  cd
  mkdir -p modulefiles/gpaw-setups
  cd modulefiles/gpaw-setups
  echo -e "#%Module1.0\nprepend-path       GPAW_SETUP_PATH    $GPAW_SETUP_SOURCE/gpaw-setups-0.9.11271" > 0.9.11271

We need to let the system know about our modules
(add this command to ``~/.profile`` or ``~/.bashrc`` to execute automatically)::

  module use $HOME/modulefiles

such that we also see them with::

  module avail

libxc
-----

GPAW relies on libxc (see the `libxc web site <https://www.tddft.org/programs/octopus/wiki/index.php/Libxc:download>`__).
To install libxc we assume that ``MYLIBXCDIR`` is set to
the directory where you want to install
(e.g. ``MYLIBXCDIR=$HOME/source/libxc``)::

 mkdir -p $MYLIBXCDIR
 cd $MYLIBXCDIR
 wget https://www.tddft.org/programs/libxc/down.php?file=5.2.2/libxc-5.2.2.tar.gz -O libxc-5.2.2.tar.gz
 tar xvzf libxc-5.2.2.tar.gz
 cd libxc-5.2.2
 mkdir install
 module purge
 module load compiler/gnu
 ./configure CFLAGS="-fPIC" --prefix=$PWD/install -enable-shared
 make |tee make.log
 make install

This will have installed the libs ``$MYLIBXCDIR/libxc-5.2.2/install/lib``
and the C header
files to ``$MYLIBXCDIR/libxc-5.2.2/install/include``.
We create a module for libxc::

  cd
  mkdir modulefiles/libxc
  cd modulefiles/libxc

and edit the module file  :file:`5.2.2` that should read::

  #%Module1.0

  #                                    change this to your path
  set             source               $::env(HOME)/source
  set             libxchome            $source/libxc/libxc-5.2.2/install
  prepend-path    C_INCLUDE_PATH       $libxchome/include
  prepend-path    LIBRARY_PATH         $libxchome/lib
  prepend-path    LD_LIBRARY_PATH      $libxchome/lib

ASE release
-----------

You might want to install a stable version of ASE::

  cd
  ASE_SOURCE=$PWD/source/ase
  mkdir -p $ASE_SOURCE
  cd $ASE_SOURCE
  git clone -b 3.18.1 https://gitlab.com/ase/ase.git 3.18.1

We add our installation to the module environment::

  cd
  mkdir -p modulefiles/ase
  cd modulefiles/ase

Edit the module file  :file:`3.18.1` that should read::

  #%Module1.0

  devel/python/3.7.10

  #           change this to your path
  set source $::env(HOME)/source
  set asehome $source/ase/3.18.1
  prepend-path       PYTHONPATH    $asehome
  prepend-path       PATH          $asehome/tools

ASE origin
----------

We get ASE origin::

  cd
  ASE_SOURCE=$PWD/source/ase
  mkdir -p $ASE_SOURCE
  cd $ASE_SOURCE
  git clone https://gitlab.com/ase/ase.git origin

which can be updated using::

  cd $ASE_SOURCE/origin
  git pull

We add our installation to the module environment::

  cd
  mkdir -p modulefiles/ase
  cd modulefiles/ase

and edit the module file  :file:`origin` that should read::

  #%Module1.0

  devel/python/3.7.10

  #           change this to your path
  set source $::env(HOME)/source
  set asehome $source/ase/origin
  prepend-path       PYTHONPATH    $asehome
  prepend-path       PATH          $asehome/tools

Building GPAW
-------------

We create a place for gpaw and get the origin version::

 cd
 GPAW_SOURCE=$PWD/source/gpaw
 mkdir -p $GPAW_SOURCE
 cd $GPAW_SOURCE
 git clone https://gitlab.com/gpaw/gpaw.git origin

The current origin version can then be updated by::

 cd $GPAW_SOURCE/origin
 git pull

A specific tag can be loaded by::

 cd $GPAW_SOURCE/origin
 # list tags
 git tag
 # load version 1.2.0
 git checkout 1.2.0

To build the current origin version of GPAW we need to create
a file :file:`siteconfig.py` that reads

.. literalinclude:: nemo_siteconfig.py

Then we build the executable::

 module purge
 module load libxc
 module load compiler/intel
 module load mpi/impi
 module load numlib/mkl
 module load ase

 cd $GPAW_SOURCE/origin
 unset CC
 python3 setup.py build

which builds GPAW to ``$GPAW_SOURCE/origin/build``.
We create a module that creates the necessary definitions::

  cd
  mkdir -p modulefiles/gpaw
  cd modulefiles/gpaw

The file  :file:`origin` that should read::

 #%Module1.0

 if {![is-loaded ase]} {module load ase}
 if {![is-loaded libxc]} {module load libxc}
 if {![is-loaded mpi]}  {module load mpi/impi}
 if {![is-loaded gpaw-setups]}  {module load gpaw-setups}

 # change the following directory definition to your needs
 set source $::env(HOME)/source
 set gpawhome $source/gpaw/origin
  # this can stay as is
 prepend-path    PATH                 $gpawhome/tools:$gpawhome/build/scripts-3.7
 prepend-path    PYTHONPATH           $gpawhome:$gpawhome/build/lib.linux-x86_64-3.7

Running GPAW
------------

A gpaw script :file:`test.py` can be submitted with the help
of :file:`gpaw-runscript` to run on 20 cpus like this::

  > module load gpaw
  > gpaw-runscript test.py 20
  using nemo
  run.nemo written
  > msub run.nemo

See options of :file:`gpaw-runscript` with::

  > gpaw-runscript -h

