.. _setuplinmac:

=========================================
Setting up the first time (Linux / macOS)
=========================================

You will be using the command line in a Terminal to access the DTU
central computers, and to set up the connections you need to display
the output on your laptop.


Mac users: Install XQuartz
==========================

Linux users should skip this step as an X11 server is part of all
normal Linux distributions.

As a Mac user, you should install an X11 server, it is needed to
display the ASE Graphical User Interface on your laptop.  If you do
not install it, you can still run Jupyter Notebooks, but the command
to view structures with the GUI will not work.

Go to https://www.xquartz.org/ and download the newest version of
XQuartz.  It is an open-source port of the X11 windows system for
macOS, it used to be part of macOS until version 10.7, but was then
removed.  Although it is no longer an official part of macOS, it is
still Apple developers that maintain the project.

After installing it, you will have to log out of your Mac and log in
again.


Open two Terminal windows
=========================

You will need two open terminal windows, where you can write Unix
commands directly to the operating system.

macOS
   The program is called ``Terminal`` and is in the folder ``Other``
   in Launchpad.  You can also find it in Spotlight by searching for
   Terminal.

Linux
   The name and placement of the program depends on the distribution.
   It is typically called ``Terminal``,  ``LXterminal``, ``xterm`` or
   something similar.  Search for Terminal in spotlight, or look in
   menus named "System tools" or similar.

Once you have opened a window you can typically get a new one by
pressing Cmd-N (Mac) or Crtl-Shift-N (Linux).


Configuring your secret key
===========================

For security reasons, all access to DTU computers require two-factor
authentification.  In the case of the HPC installation, the two
factors are your password and an encryption key (*SSH key*).  You get
your secret SSH key as described in the document with your username
and password.  Save it on your laptop.

You now need to install the key so it is used automatically.  This unfortunately
involves editing a configuration file in a hidden folder, which
possibly does not exist yet.  This can be done in a terminal window.
In a terminal window, run the following three commands exactly as
written here.  Respect the spaces, and be aware that the folder name
``.ssh`` starts with a period.::

  cd
  mkdir -p .ssh
  touch .ssh/config
  open .ssh

The first line switches to your home folder in the unlikely case you
are not there already.  The second line creates the hidden ``.ssh``
folder, if it does not already exist.  The third line creates an empty
configuration file for the SSH program, if a file already exists it
merely updates the time stamp.  The fourth line opens a Finder window
in the hidden folder (you cannot find it otherwise with Finder), it is
MacOS specific.  On a Linux machine you should instead open the file
browser and find the .ssh folder with it, normally you can access
hidden folders without difficulty on a Linux machine.

Now, use the Finder to move the SSH key file that you just downloaded
into this hidden folder.  Then, right-click on the file named
``config`` and open it with TextEdit or another text editor.  The file may already contain
stuff or it may be empty.  Add the following to the end of the file
(after a blank line, if there is already stuff in the file)::

  Host *.gbar.dtu.dk
      User XXXXXXXX
      IdentityFile ~/.ssh/YYYYYYY
      ForwardX11 yes
      ForwardX11Timeout 36000

where XXXXXXX should be replaced with your DTU user name, and YYYYYYY
with the name of the file containing your secret SSH key (there are
two files, you should give the name of the one that does **NOT**  end
in .pub).  Save the config file.

Finally, you need to make sure your key file is not readable by other
users on your computer - even if you are the only user, the SSH
program will only use the key file if the permissions are sufficiently
restrictive.  In the Terminal window, write::

  chmod og-rwx .ssh/YYYYYYYY

where YYYYYYY is the name of the key file.  The strange string
``og-rwx`` contains a minus in the middle, and removes read, write and
eXecute permissions for *other users* and *group users*.  Check that
permissions are set correctly::

  ls -l .ssh

The output should look something like::

  total 384
  -rw-r--r--@ 1 myname  staff    735 Jul  5 11:27 config
  -rw-r-----@ 1 myname  staff   1315 Jul  2 14:11 config.old
  -rw-r--r--@ 1 myname  staff   1315 Jul  2 14:11 config~
  -rw-------  1 myname  staff   1843 Jul  4 15:32 gbartest_rsa
  -rw-r--r--  1 myname  staff    413 Jul  4 15:32 gbartest_rsa.pub
  -rw-------  1 myname  staff   1675 Sep  8  2010 id_rsa
  -rw-r--r--  1 myname  staff    401 Sep  8  2010 id_rsa.pub
  -rw-------  1 myname  staff  80395 Jul  1 15:55 known_hosts

Check that the line for the sectet key file file starts with ``-rw-------``
(there may be a @ or another character at the end, that is OK).

Your key file should now automatically be used when logging into the
DTU data bar.



Log into the databar
====================

You use the ``ssh`` (Secure SHell) command to create a secure
(encrypted) connection to the databar computers.  In the terminal,
write::

  ssh -XY USERNAME@login.gbar.dtu.dk

where ``USERNAME`` is your DTU user name (external participants got it
in their registration material).  You can leave it out if you have
written it in your ``.ssh/config`` file.  Note the ``-XY`` option, it is a minus
followed by a capital X and a capital Y, it tells ssh to let the
remote computer open windows on your screen.

Note that when you write your DTU password, you cannot see what you
type (not even stars or similar!).

You now have a command-line window on the DTU login-computer.  This
computer (``gbarlogin`` a.k.a. ``login.gbar.dtu.dk``) may not be used
to calculations, as it would be overloaded.  You therefore need to log
in to the least loaded interactive computer by writing the command::

  linuxsh -X

You now have a command-line window on an interactive compute node, as shown
below.

.. image:: Logged_in_Mac.png
   :width: 66%

The two last lines are the command prompt.  The first line indicates
your current working directory, here your home folder symbolized by
the ~ (tilde) character.  The lower line gives the name of the
computer (here ``n-62-27-23``) and the user name (``jasc`` in the figure)
followed by a dollar sign.



Get access to the software
==========================

To give access to the software you need for this course, please run
the command::

  source ~jasc/setup2024

Note the tilde in the beginning of the second word.
The script installs ASE, GPAW and related software.
The script will also copy a selection of draft notebooks to a folder
called CAMD2024 in your DTU databar account.


Carrying on
===========

Now read the guide for :ref:`Starting and accessing a Jupyter Notebook
<accesslinmac>`

