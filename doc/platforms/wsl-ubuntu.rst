.. _wsl:
.. index:: WSL

=================
WSL Ubuntu 20.04+
=================

This guide assumes that you have a running Ubuntu_ installation in Windows Subsystem for Linux (WSL). There are two WSL types: WSL1 and WSL2. WSL1 is older but gives more performance in the calculations. However, WSL2 is more supported by Microsoft, and the future of WSL will be WSL2 with WSLg (graphical interface). You must be running Windows 10 version 2004 and higher (Build 19041 and higher) or Windows 11. More information about installation of WSL can be found from Microsoft WSL Installation_ website. After installing WSL and the Ubuntu distribution, make sure that things are up-to-date in Ubuntu::

    $ sudo apt update
    $ sudo apt upgrade

There are two ways to use GPAW on WSL: Classical installation and Conda installation.


Classical Installation
======================

You need Tk library for GUI, unzip for file unzipping, and for further package installations, we need PIP installer::

    $ sudo apt install python3-tk python3-pip unzip python-is-python3

Install ASE and other math, parallel, dev libraries::

    $ pip3 install --upgrade --user base


At this point, PIP can give some warnings as:

    WARNING: The scripts f2py, f2py3, and f2py3.8 are installed in '/home/YOURUSERNAME/.local/bin' which is not on PATH.
    Consider adding this directory to PATH or, if you prefer to suppress this warning, use --no-warn-script-location.
    WARNING: The scripts ase, ase-build, ase-db, ase-gui, ase-info and ase-run are installed in
    '/home/YOURUSERNAME/.local/bin'
    which is not on PATH.
    Consider adding this directory to PATH or, if you prefer to suppress this warning, use --no-warn-script-location.

Add the following line at the end of your ~/.bashrc file::

    export PATH=/home/YOURUSERNAME/.local/bin:$PATH

After editing ~/.bashrc file quit the current shell session and start a new one (or you can use source ~/.bashrc command). Then continue::

    $ sudo apt install python3-dev libopenblas-dev libxc-dev libscalapack-mpi-dev libfftw3-dev

Create a siteconfig.py file::

    $ mkdir -p ~/.gpaw
    $ cat > ~/.gpaw/siteconfig.py
    fftw = True
    scalapack = True
    libraries = ['xc', 'blas', 'fftw3', 'scalapack-openmpi']
    ^D

**NOTE:** If the user wants to use exchange correlations listed in libxc library, ‘xc’ must be listed in the libraries line as shown above. Otherwise, only built in exchange correlation functionals can be used.

Then install gpaw (it will install with its dependencies: ASE_, Numpy, SciPy)::

    $ pip3 install --upgrade --user gpaw

Use gpaw info to see information about installation. However, PAW-datasets are not installed yet. To install it, firstly create a directory under ~/.gpaw then install PAW datasets::

    $ mkdir ~/.gpaw/gpaw-setups
    $ gpaw install-data ~/.gpaw/gpaw-setups/


Conda Installation (easier than the Classical installation)
===========================================================

Download and install the miniconda. Say 'yes' or 'no' to initialization after installing it::

    $ wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    $ bash Miniconda3-latest-Linux-x86_64.sh

or::

    $ ./Miniconda3-latest-Linux-x86_64.sh

then you can update miniconda::

    $ eval "$(/home/$USER/miniconda3/bin/conda shell.bash hook)"
    $ conda update conda

Now, we can create an environment (here 'gpaw-env' name is used. You can use any name) and activate it::

    $ conda create --name gpaw-env
    $ conda activate gpaw-env

Lastly, you can install GPAW::

    $ conda install -c conda-forge gpaw

.. _Ubuntu: https://www.ubuntu.com/
.. _ASE: https://wiki.fysik.dtu.dk/ase/
.. _Installation: https://docs.microsoft.com/en-us/windows/wsl/install
