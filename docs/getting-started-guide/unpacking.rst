Unpacking the distribution files
================================

The rDock source files and test suite files are provided as independent gzipped
tar (.tar.gz) distributions. Depending on your requirements, the two
distributions can be unpacked to entirely separate locations, or can be unpacked
under the same location. In this example they are unpacked under the same
location.

.. table:: rDock distribution files

   +---------------------------------+-----------------------------------+
   | File                            | Description                       |
   +=================================+===================================+
   | ``rDock_[CODELINE]_src.tar.gz`` | rDock source distribution         |
   +---------------------------------+-----------------------------------+
   | ``[TEST]_rDock_TestSet.tar.gz`` | Test suite data files and scripts |
   +---------------------------------+-----------------------------------+

Here ``[CODELINE]``, and ``[TEST]`` will vary depending on the release and test
set. ``[CODELINE]`` represents the major version string (for example, 2013.1)
and ``[TEST]`` represents the given dataset (ASTEX, RNA or DUD).

Example unpacking procedure
---------------------------

Create a new directory for building rDock.

.. code-block:: bash

   $ mkdir ~/dev

The directory you created is referred to as ``[BUILDDIR]`` in the subsequent
steps.

Copy or download the distribution files to ``[BUILDDIR]``.

.. code-block:: bash

   $ cp ~/mydownloads/rDock_2013.1_src.tar.gz ~/dev/

Extract the distributions.

.. code-block:: bash

   $ cd ~/dev/
   $ tar -xvzf rdock_2013.1_src.tar.gz

The distributions contain files with relative path names, and you should find
the following subdirectories created under ``rDock_[CODELINE]_src``. Note that
the ``./rDock_2013.1_src`` subdirectory may have a different name depending on
the major version string (see above).

.. code-block:: bash

   $ find . -type d
   .
   ./fw
   ./src
   ./src/daylight
   ./src/lib
   ./src/exe
   ./src/GP
   ./build
   ./build/test
   ./build/test/RBT_HOME
   ./build/tmakelib
   ./build/tmakelib/linux-pathCC-64
   ./build/tmakelib/linux-g++-64
   ./build/tmakelib/linux-g++
   ./build/tmakelib/unix
   ./data
   ./data/filters
   ./data/sf
   ./data/pmf
   ./data/pmf/smoothed
   ./data/scripts
   ./lib
   ./import
   ./import/tnt
   ./import/tnt/include
   ./import/simplex
   ./import/simplex/src
   ./import/simplex/include
   ./docs
   ./docs/images
   ./docs/newDocs
   ./include
   ./include/GP
   ./bin

Make a note of the following locations for later use.

The rDock root directory is ``[BUILDDIR]/rDock_[CODELINE]_src`` and will be
referred to as ``[RBT_ROOT]`` in later instructions. In this example,
``[RBT_ROOT]`` is ``~/dev/rDock_2013.1_src/``.
