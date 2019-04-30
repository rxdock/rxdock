Prerequisites
=============

Compilers
---------

rDock is supplied as source code, which means that you will have to compile the
binary files (run-time libraries and executable programs) before you use them.
rDock has been developed largely on the Linux operating systems, most recently
with GNU g++ compiler (tested under openSUSE 11.3). The code will almost
certainly compile and run under other Linux distributions with little or no
modification. For the moment, it has been tested in the latest Ubuntu and
openSUSE releases for both 32 and 64 bits system architectures (by
November 2013) and compilation was possible without any code modification.
However, no other distributions or compilers have been tested extensively to
date.

For full production use, you would typically compile rDock on a separate build
machine and run the docking calculations on a cluster of compute machines.
However, for the purposes of getting started, these instructions assume that you
will be compiling rDock and running the initial validation experiments on the
same machine.

Required packages
-----------------

Make sure you have the following packages installed on your machine before you
continue. The versions listed are appropriate for openSUSE 11.3; other versions
may be required for your particular Linux distribution.

.. table:: Required packages for building and running rDock

   +-------------------+----------------------------+--------------+-----------+
   | Package           | Description                | Required at  | Version   |
   +===================+============================+==============+===========+
   | ``gcc``           | GNU C compiler             | Compile-time | >=3.3.4   |
   +-------------------+----------------------------+--------------+-----------+
   | ``g++``           | GNU C++ compiler           | Compile-time | >=3.3.4   |
   +-------------------+----------------------------+--------------+-----------+
   | ``cppunit``       | C++ unit test framework    | Compile-time | >=1.10.2  |
   +-------------------+----------------------------+--------------+-----------+
   | ``cppunit-devel`` | Development files for      | Compile-time | >=1.10.2  |
   |                   | cppunit                    |              |           |
   +-------------------+----------------------------+--------------+-----------+
