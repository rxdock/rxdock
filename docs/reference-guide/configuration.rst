Configuration
=============

Before launching |Dock|, make sure the following environment variables are
defined. Precise details are likely to be site-specific.

* ``RBT_ROOT`` environment variable: should be defined to point to the |Dock|
  installation directory.
* ``RBT_HOME`` environment variable: is optional, but can be defined to to point
  to a user project directory containing |Dock| input files and customised data
  files.
* ``PATH`` environment variable: ``$RBT_ROOT/bin`` should be added to the
  ``$PATH`` environment variable.
* ``LD_LIBRARY_PATH``: ``$RBT_ROOT/lib`` should be added to the
  ``$LD_LIBRARY_PATH`` environment variable.

Input file locations
--------------------

The search path for the majority of input files for |Dock| is:

* Current working directory
* ``$RBT_HOME``, if defined
* The appropriate subdirectory of ``$RBT_ROOT/data/``. For example, the default
  location for scoring function files is ``$RBT_ROOT/data/sf/``.

The exception is that input ligand SD files are always specified as an absolute
path. If you wish to customise a scoring function or docking protocol, it is
sufficient to copy the relevant file to the current working directory or to
``$RBT_HOME``, and to modify the copied file.

Launching executables
---------------------

For small scale experimentation, the |Dock| executables can be launched directly
from the command line. However, serious virtual screening campaigns will likely
need access to a compute farm. In common with other docking tools, |Dock| uses
the embarrassingly parallel approach to distributed computing. Large ligand
libraries are split into smaller chunks, each of which is docked independently
on a single machine. Docking jobs are controlled by a distributed resource
manager (DRM) such as Condor or SGE.
