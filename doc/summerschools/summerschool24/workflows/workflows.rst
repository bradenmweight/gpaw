========================================
Workflows with Atomic Simulation Recipes
========================================

In this exercise we will write and run computational workflows
using the atomic simulation recipes, ASR.

The basic unit of computation in ASR is a *task*.  A task is a Python
function along with a specification of the inputs to that function.
The inputs can be either concrete values like lists, strings, or numbers,
or references to the outputs of other tasks.
Tasks that depend on one another and form a graph.  An
ASR workflow is a Python class which defines such a graph, along with
metadata about how the tasks should run and be stored.  Workflows can
then be parametrized so they run on a collection of materials,
for example.

When using ASR, we define workflows and tasks using Python.
However, the tools used to *run* workflows and tasks are command-line
tools.  Therefore, for this exercise we will be using the terminal
rather than the usual notebooks.  Basic knowledge of shell commands
is an advantage.


This exercise consists of three parts.  Specifically, we will:

 * Write a workflow which defines a structure optimization task
 * Extend the workflow with ground state and band structure tasks
 * Parametrize the workflow to apply it to multiple materials

When actually using ASR, many tasks and workflows are already written.
Thus, we would be able to import and use those features directly.
But in this tutorial we write everything from scratch.


Part 1: Create a repository and define a workflow
=================================================

First, go to a clean directory and create an ASR repository::

 human@computer:~$ mkdir myworkflow
 human@computer:~$ cd myworkflow
 human@computer:~/myworkflow$ asr init
 Created repository in /home/askhl/myworkflow

The repository will store calculations under the newly created,
currently empty folder named ``tree/``.  The ``asr info`` command
will tell us a few basic things about the repository::

 human@computer:~/myworkflow$ asr info
 Root:     /home/askhl/myworkflow
 Tree:     /home/askhl/myworkflow/tree
 db-file:  /home/askhl/myworkflow/registry.dat (0 entries)
 Tasks:    /home/askhl/myworkflow/tasks.py (not created)

Let's perform a structure optimization of bulk Si.
We write a function which performs such an optimization:

.. literalinclude:: tasks.py
   :end-before: end-snippet-1

This function uses a cell filter to expose the cell degrees of freedom
for the standard BFGS optimizer (see the ASE documentation on optimizers
and cell filters if interested).

Since workflows run on the local computer whereas computational tasks
(generally) run on compute nodes, we separate *workflow* code
and *computational*
code in different files.  ASR can load user-defined functions from the
special file ``tasks.py`` mentioned by info command.
Create that file and save the above function to it.

Next, we write a workflow with a task that will call the function:

.. literalinclude:: workflow.py
   :end-before: end-snippet-1

Explanation:

* The ``@asr.workflow`` decorator tells ASR to regard the class as a
  workflow.  In particular, it equips the class with a constructor
  with appropriate input arguments.

* ``asr.var()`` is used to declare input variables.  The names
  ``atoms`` and ``calculator`` imply
  that we want this workflow to take atoms and calculator parameters
  as input.

* The method ``relax()`` defines our task.
  By naming the method ``relax()``, we choose that the task will run
  in a directory called ``tree/relax``.

* The method returns ``asr.node(...)``, which is a specification of
  the *actual* calculation: The name of the task
  (``'relax'``, which must exist in ``tasks.py``) is given as a string.
  The inputs are then assigned, and will be forwarded to the
  ``relax()`` function in ``tasks.py``.
  The attributes ``self.atoms`` and ``self.calculator``
  will refer to the input variables.

  When defining a node, ASR calculates a hash (i.e. checksum) of the inputs;
  the hash will become different if any inputs are changed.

* The decorator ``@asr.task`` can be used to attach information
  about *how* the task runs, such as computational
  resources.


The workflow class serves as a *static declaration* of information, not as
statements or commands to be executed (yet).
To actually run it, we must at least choose a material and then tell
the computer to run the workflow on it.
We do this by adding a standalone function called ``workflow``
for ASR to call:

.. literalinclude:: workflow.py
   :pyobject: workflow


ASR will take care of creating a "runner" and passing it to the function.
(Note: In a future version of the code, this syntax will be simplified.)


Save the code (both the class and the ``workflow()`` function)
to a file, e.g. named ``workflow.py``.  Then execute the
workflow by issuing the command::

  asr workflow workflow.py

The command executes the workflow and creates a folder under the ``tree/``
directory for each task.
We can run ``asr ls`` to see a list of the tasks we generated::

  541d427c new      tree/relax                     relax(atoms=…, calculator=…)

The task is identified to the computer as a hash value (541d427c....),
whereas to a human user, the location in the directory tree,
``tree/relax``, will be more descriptive.

Feel free to look at the contents of the ``tree/relax`` directory.
The task is listed as "new" because we did not run it yet
— we only *created* it, so far.  While developing workflows,
we will often want to
create and inspect the tasks before we submit anything expensive.
If we made a mistake, we can remove the task with
``asr remove tree/relax``, then fix the mistake and run the workflow again.

Once we're happy with the task, let's run the task on the local computer::

  asr run tree/relax

If everything worked as intended, the task will now be "done",
which we can see by running ``asr ls`` again::

  541d427c done     tree/relax                     relax(atoms=…, calculator=…)


We can use the very handy ``tree`` command to see the whole
directory tree::

  human@computer:~/myworkflow$ tree tree/
  tree/
  └── relax
      ├── gpaw.txt
      ├── input.json
      ├── input.resolved.json
      ├── opt.log
      ├── opt.traj
      ├── output.json
      └── state.dat

  1 directory, 7 files

Be sure to open the trajectory file in (e.g. in ASE GUI) to check
that the optimization ran as expected.  Also the logfiles
``gpaw.txt`` and ``opt.log`` are there.


Part 2: Add ground state and band structure tasks
=================================================

After the relaxation, we want to run a ground state
calculation to save a ``.gpw`` file, which we subsequently want
to pass to a non self-consistent calculation to get the band structure.

Add a ``groundstate()`` function to ``tasks.py``:

.. literalinclude:: tasks.py
   :pyobject: groundstate

In order to "return" the gpw file, we actually return a ``Path`` object
pointing to it.  When passing the path to another task, ASR
resolves it with respect to the task's own directory such
that the human will not need to remember or care about the actual directories
where the tasks run.

Let's add a corresponding groundstate method to the workflow:

.. literalinclude:: workflow.py
   :pyobject: MyWorkflow.groundstate


By calling ``asr.node(..., atoms=self.relax)``, we are specifying
that the atoms should be taken as the *output* of the ``relax`` task,
creating a dependency.

We can now run the workflow again.  The old task still exists and
will remain unchanged, whereas the new task should now appear
in the ``tree/groundstate`` directory.

Run the ground state task and check that the ``.gpw`` file was created as
expected.

Finally, we write a band structure task in ``tasks.py``:

.. literalinclude:: tasks.py
   :pyobject: bandstructure

A corresponding method should be added on the workflow:

.. literalinclude:: workflow.py
   :pyobject: MyWorkflow.bandstructure

Now run the workflow and the resulting tasks.
The code saves the Brillouin zone path and band structure separately to
ASE JSON files.  Once it runs, we can go to the directory and check
that it looks correct::

  ase reciprocal tree/bandstructure/bandpath.json

::

   ase band-structure tree/bandstructure/bs.json

Note that here we are using the ``ase`` tool, *not* the ``asr`` tool.

You can delete all the tasks with ``asr remove tree/`` and run them from
scratch by ``asr run tree/``, ``asr run tree/*``, or simply ``asr run
tree/bandstructure``.
The run command always executes tasks in
topological order, i.e., each task runs only when its dependencies
are done.

The ``asr ls`` command can also be used to list tasks in topological
order following the dependency graph::

  human@computer:~/myworkflow$ asr ls --parents tree/bandstructure/
  541d427c done     tree/relax                     relax(atoms=…, calculator=…)
  5ca14caa done     tree/groundstate               groundstate(atoms=<541d427c>, calculator=…)
  b5875ebd done     tree/bandstructure             bandstructure(gpw=<5ca14caa>)


This way, we can comfortably work with larger numbers of tasks.
Note how the hash values are consistent:
The band structure's input includes the hash value of the
ground state, and the ground state's input includes the hash value of the
relaxation.

If we edit the workflow such that tasks receive different inputs,
then the hash values will change, and ASR will raise an error
because the new hash is inconsistent with the old one in that directory.
Such a conflict can be solved by removing the old calculations.


Part 3: Run workflow on multiple materials
==========================================

The current workflow creates directories right under the repository root.
For a proper materials workflow, it will be helpful to work
with a structure that nests the tasks by material.

ASR contains a feature called ``totree`` which deploys a dataset
to the tree, such as defining initial structures for materials.
One then parametrizes a workflow (such as the one we just wrote)
on the materials.

The following workflow defines a function which returns a set of materials,
then specifies to ASR that those must be added to the tree.

.. literalinclude:: totree.py

Add this to a new file, named e.g. ``totree.py``, and execute the workflow::

 human@computer:~/myworkflow$ asr workflow totree.py
       Add: 889575c5 new      tree/Al/material               define(obj=…)
       Add: 5e39fb8e new      tree/Si/material               define(obj=…)
       Add: 9612a07a new      tree/Ti/material               define(obj=…)
       Add: 7153df81 new      tree/Cu/material               define(obj=…)
       Add: 155d59ee new      tree/Ag/material               define(obj=…)
       Add: e9b41657 new      tree/Au/material               define(obj=…)

The totree command created some tasks for us.
Actually they are not really tasks — they are just static pieces of data.
But now that they exist, we can run other tasks that depend on them.

In the old workflow file (``workflow.py``),
replace the ``workflow()`` function with the following function which
tells ASR to parametrize the workflow by "globbing" over the materials:


.. literalinclude:: materials.py
   :pyobject: workflow

The workflow will now be called once for each material.
Run the workflow and it will create our three well-known tasks
for each material, now nested by material.

As before, we can inspect the newly created tasks, e.g.::

 human@computer:~/myworkflow$ asr ls tree/Au/bandstructure/ --parents
 e9b41657 new      tree/Au/material               define(obj=…)
 5306d226 new      tree/Au/relax                  relax(atoms=<e9b41657>, calculator=…)
 a54f98a7 new      tree/Au/groundstate            groundstate(atoms=<5306d226>, calculator=…)
 7fbfa099 new      tree/Au/bandstructure          bandstructure(gpw=<a54f98a7>)


Since it may take a while to run on the front-end node,
we can tell ASR to submit one or more tasks using MyQueue_::

  asr submit tree/Au

The submit command works much like the run command, only it calls
myqueue which will then talk to the scheduler (slurm, torque, ...).
After submitting, we can use standard myqueue commands to monitor
the jobs, such as ``mq ls`` or ``mq rm``.  See the `myqueue documentation
<https://myqueue.readthedocs.io/en/latest/cli.html>`_.

If everything works well, we can submit the whole tree::

  asr submit tree/

Note: In the current version, myqueue and ASR do not perfectly
share the state of a task.  This can lead to
mild misbehaviours if using both ``asr run`` and ``asr submit``,
such as a job executing twice.


.. _MyQueue: https://myqueue.readthedocs.io/
