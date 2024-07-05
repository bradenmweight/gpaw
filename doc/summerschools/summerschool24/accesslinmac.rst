.. _accesslinmac:

=========================================================
Starting and accessing a Jupyter Notebook (Linux / macOS)
=========================================================

To run a Jupyter Notebook in the DTU databar while displaying it output in your browser requires three steps.

* Starting the notebook on an interactive compute node.

* Make a connection to the relevant compute node, bypassing the firewall.

* Connecting your browser to the Jupyter Notebook process.


Logging into the databar
========================


If you are not already logged into the databar, do so by starting
a Terminal window.  Log in to the databar front-end with the command::

  ssh -XY USERNAME@login.gbar.dtu.dk

Replace ``USERNAME`` with your DTU username.

Once you are logged in on the front-end, get a session on an interactive compute node by typing the command::

  linuxsh -X


Starting a Jupyter Notebook
===========================

Change to the folder where you keep your notebooks (most likely ``CAMD2024``) and start the Jupyter Notebook server::

  cd CAMD2024
  camdnotebook

The command ``camdnotebook`` is a local script.  It checks that you
are on a compute server (and not on the front-end) and that X11
forwarding is enabled.  Then it starts a jupyter notebook by running
the command ``jupyter notebook --no-browser --port=40000 --ip=$HOSTNAME``
(you can also use this command yourself if you prefer).

The Notebook server replies by printing a few status lines, as seen here

.. image:: JupyterRunningMac.png
   :width: 66%

The important line is the second from the bottom, it shows on which
computer and port number the notebook is running (here ``n-62-30-6``
and 40000, respectively).


Create an SSH Tunnel to the notebook
====================================

You now need to create an SSH tunnel to the server directly from your laptop.  It is also done with an ``ssh`` command, which unfortunately is a bit cryptic.  *Open a new Terminal window on your laptop,*  and write the following command::

  ssh USERNAME@login.gbar.dtu.dk -g -LPORT:HOSTNAME:PORT -N

I this line, you should replace ``USERNAME`` with your DTU username,
``HOSTNAME`` with the servername you see in the other terminal window
(it has the form ``n-XX-YY-ZZ`` or ``n-XX-YY-Z``) and ``PORT`` with the port number you
see in that line (typically 40000 or close, note that it should be
specified in two places!).  There is no feedback in form of stars when
you type the password.  If you type the password correctly (and press
enter) then the command gives *no feedback indicating that it is
running!*


Starting a browser.
===================

Start a browser (Chrome and Firefox are known to work well).  Find the
link that begins with ``http://127.0.0.1`` and paste it into the
address bar.  It will look similar to this::

  http://127.0.0.1:40000/?token=401bb4a3a4faeafd2fd948a137b0f6857ba4aa6e6fc47b7d

You are now ready to open one of the notebooks, and run the exercises.

Logging out
===========

When you are done for the day, please

* Save all notebooks, then select ``Close and Halt`` on the file menu.

* Stop the SSH tunnel.

* Stop the Jupyter Notebook server by pressing Control-C twice in the
  window where it is running.

* Log out of the databar by typing ``exit`` twice in the window(s).
