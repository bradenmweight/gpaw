.. _setupwin:

===================================
Setting up the first time (Windows)
===================================

You will need the program MobaXterm which will help you with the
following tasks

  * Logging in to the databar.

  * Displaying the ASE graphical user interface on your laptop.

  * Forward a connection to the Jupyter server (SSH tunnel).

  * Accessing files on the remote server from your laptop.


Installing MobaXterm
====================

You download the program from the website
https://mobaxterm.mobatek.net.  Choose Download, select then "Free
Home Edition".  You are now given a choice between a "Portable
Edition" and an "Installer Edition".  The Installer Edition is
installed like any other Windows program; the Portable Edition comes
as a ZIP file that you need to unpack.  The program is then in a
folder together with a data file, and you can run it from this
folder.  We have tested the installer version, but assume that both
work.

You must unpack the downloaded ZIP file before you run the installer,
running it from inside the ZIP does not work!



Configuring your secret key
===========================

For security reasons, all access to DTU computers require two-factor
authentification.  In the case of the HPC installation, the two
factors are your password and an encryption key (*SSH key*).  You get
your secret SSH key as described in the document with your username
and password.  Save it on your laptop.

Start MobaXterm.  You will see a window with a row of buttons at the
top.  Click on the Settings button (one of the last ones).  On the
configuration window that opens (shown below) you select the SSH tab.
The last box is labelled "SSH agents", select "Use internal SSH agent
MobAgent".

Click the Plus sign next to the box labelled "Load the following keys
at MobAgent startup" and then select the key with your private SSH key
(i.e. **not** the file with file type .pub).  Then press OK; MobaXterm
will then need to restart.

.. image:: MobaXtermKeyConfiguration.png
   :width: 66%

Connecting the first time
=========================

Start MobaXterm.  You will see a window with a row of buttons at the
top.  Click on the *Session* button, you will now see a window as
shown below.

.. image:: Moba_ssh.png
   :width: 66%

In the tab *Basic SSH settings* you should choose *Remote host* to be
``login.gbar.dtu.dk``.  The user name is your DTU user name (external
participants got it in the registration material).  The port number
must remain 22.  Click *OK*  and give your DTU password in the text
window when prompted.  **NOTE** Nothing is written when you type the
password, not even stars.

**We do not recommend allowing MobaXterm to remember your password!**

You now have a command-line window on the DTU login-computer, as shown
below.

.. image:: Logged_in_Win.png
   :width: 66%

The two last lines are the command prompt.  The first line indicates
your current working directory, in this case your home folder symbolized by
the ~ (tilde) character.  The lower line gives the name of the
computer (``gbarlogin``) and the user name (``jasc`` in the figure)
followed by a dollar sign.

This computer (``gbarlogin``) may not be used for calculations, as it
would be overloaded.  You therefore need to log in to the least loaded
interactive computer by writing the command::

  linuxsh -X

(the last X is a capital X, you get no error message if you type it
wrong, but the ASE graphical user interface will not work).


Get access to the software
==========================

To get access to the software you need for this course, please run
the command::

  source ~jasc/setup2024

Note the tilde in the beginning of the second word.

The script will install ASE, GPAW and related software.
The script will also copy a selection of draft notebooks to a folder
called CAMD2024 in your DTU databar account.


Carrying on
===========

Now read the guide for :ref:`Starting and accessing a Jupyter Notebook
<accesswin>`

