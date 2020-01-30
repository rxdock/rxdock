Building
========

|Dock| is written in C++ (with a small amount of C code from Numerical Recipes)
and makes heavy use of the C++ Standard Template Library (STL). The majority of
the source code is compiled into a single shared library (``libRbt.so``). The
executable programs themselves are relatively light-weight command-line
applications linked with ``libRbt.so``.

The tmake build systems (from Trolltech) is used to generate makefiles
automatically for a particular build target (i.e. combination of operating
system and compiler). The source distribution comes with tmake templates
defining the compiler options and flags for three Linux build targets
(``linux-g++`` and ``linux-g++-64``). The build targets have been tested under
openSUSE 11.3 (2.6.34.10-0.2 kernel) with GNU g++ (versions 3.3.4, 4.5.0,
and 4.7.2).

.. table:: Standard tmake build targets provided

   +------------------+--------------+----------+--------------------------------+
   | Target Name      | Architecture | Compiler | Compiler flags (release build) |
   +==================+==============+==========+================================+
   | ``linux-g++``    | 32-bit       | ``g++``  | ``-m32 -O3 -ffast-math``       |
   +------------------+--------------+----------+--------------------------------+
   | ``linux-g++-64`` | 64-bit       | ``g++``  | ``-m64 -O3 -ffast-math``       |
   +------------------+--------------+----------+--------------------------------+

Customising the tmake template for a build target
-------------------------------------------------

If none of the tmake templates are suitable for your machine, or if you wish to
customise the compiler options, you should first customise one of the existing
templates. The tmake template files are stored under
``[RBT_ROOT]/build/tmakelib/``. Locate and edit the ``tmake.conf`` file for the
build target you wish to customise. For example, to customise the ``linux-g++``
build target, edit ``[RBT_ROOT]/build/tmakelib/linux-g++/tmake.conf`` and
localise the values to suit your compiler.

Build procedure
---------------

To build |Dock|, first go to the ``[RBT_ROOT]/build/`` directory.

.. code-block:: bash

  $ cd [RBT_ROOT]/build

Compile
^^^^^^^

Make one of the build targets listed below.

.. code-block:: bash

   $ make linux-g++
   $ make linux-g++-64

Test
^^^^

Run the |Dock| unit tests to check build integrity. If no failed tests are reported you should be all set.

.. code-block:: bash

   $ make test

Cleanup (optional)
^^^^^^^^^^^^^^^^^^

To remove all intermediate build files from ``[RBT_ROOT]/build/``, leaving just
the final executables (in ``[RBT_ROOT]/bin/``) and shared libraries (in
``[RBT_ROOT]/lib/``):

.. code-block:: bash

   $ make clean

To remove the final executables and shared libraries as well, returning to a
source-only distribution:

.. code-block:: bash

   $ make distclean
