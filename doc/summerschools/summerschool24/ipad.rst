.. _ipad:

==============
 Using an iPad
==============

Using an iPad for the computer exercises is **not recommended,** mainly because you cannot run the ASE graphical user interface.  But it is a possibility if there are no other options.

Terminal app
============

A terminal app that is known to work is Termius.  Create a login by choosing ``Hosts`` in the sidebar, then select the Plus sign at the top.  Select "New Host"

* Alias: DTU gbar (or whatever you want)
* Hostname: login.gbar.dtu.dk
* Click on the person icon where it says Username, then choose New Identity.  Set Username and Password to your DTU username and password.  Save and click the back arrow.  Save the profile

Click on the new profile, you should be logged in.

Now follow the instructions for :ref:`logging in the first time with a Mac <setuplinmac>`, with the important change that the command ``linuxsh -X`` should be replaced with ``linuxsh``.  The omitted ``-X`` makes it possible for the graphical user interface to show atomic structures on your screen, that does not work on an iPad and a workaround is described below.

Connecting to Jupyter Notebooks
===============================

Then follow the guide for :ref:`Starting and accessing a Jupyter Notebook <accesslinmac>`  When you start the Jupyter server with the ``camdnotebook`` command, you will be warned that X11 forwarding is not set up.  Ignore the warning and continue by pressing ENTER.

You now need to set up port forwarding.  Make a note of the computer name and host number, i.e. n-62-27-19 and 40042.  Then click on the button with a < sign in the lower right corner.

* Click on Port Forwarding.
* Select type "local" and click on Continue 
* Set local port to 8080, leave the other field blank.  Press Continue
* Select a host.  Choose the same setup that you made when you logged in.
* Destination address is the host name you noted, i.e. n-62-27-19 (it will be different for you!)
* Destination port is the port number you noted, i.e. 40042 (it will also be different for you)
* Click Done, then Save

Then click on the new forwarding button to start it.

**You cannot edit the forwarding rule.  When you log out and in again and get a new host and portnumber, you have to define a new rule.  Make sure the old one is not running (maybe even delete it).**

Now click on your terminal connection to go back to the terminal window.

Using the notebook
==================

You need to have Termius and your browser running simultaneously.  Click on the three dots on the top of the window, and choose the split view.  Then open your browser.  You should now have Termius and the browser running side by side.  In the address field of the browser, type ``localhost:8080`` to connect to the notebook.  You should see a line being printed in the terminal, and the browser asks for your *jupyter*  password (not the DTU password!).

Viewing atomic structures
=========================

The ``view(atoms)`` command in ASE will not work for you, as you do not have an  X11 server on your ipad.  Instead use the ``plot_atoms`` command::

  from ase.visualize.plot import plot_atoms
  plot_atoms(atoms)

You rotate the atoms like this (45 degrees along the x-axis)::

  plot_atoms(atoms, rotation='45x')


  
  

  
