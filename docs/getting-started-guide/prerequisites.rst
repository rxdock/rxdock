Prerequisites
=============

Compilers
---------

|Dock| is supplied as source code, which means that you will have to compile the
binary files (run-time libraries and executable programs) before you use them.
|Dock| has been developed largely on the Linux operating systems, most recently
with GNU g++ compiler (tested under openSUSE 11.3). The code will almost
certainly compile and run under other Linux distributions with little or no
modification. For the moment, it has been tested in the latest Ubuntu and
openSUSE releases for both 32 and 64 bits system architectures (by
November 2013) and compilation was possible without any code modification.
However, no other distributions or compilers have been tested extensively to
date.

For full production use, you would typically compile |Dock| on a separate build
machine and run the docking calculations on a cluster of compute machines.
However, for the purposes of getting started, these instructions assume that you
will be compiling |Dock| and running the initial validation experiments on the
same machine.

Required packages
-----------------

Make sure you have the following packages installed on your machine before you
continue. The versions listed are appropriate for openSUSE 11.3; other versions
may be required for your particular Linux distribution.

.. table:: Required packages for building and running |Dock|

   +-------------------+----------------------------+--------------+-----------+
   | Package           | Description                | Required at  | Version   |
   +===================+============================+==============+===========+
   | ``gcc``           | GNU C compiler             | Compile-time | >=3.3.4   |
   +-------------------+----------------------------+--------------+-----------+
   | ``g++``           | GNU C++ compiler           | Compile-time | >=3.3.4   |
   +-------------------+----------------------------+--------------+-----------+
   | ``gtest``         | Google's C++ unit testing  | Compile-time | >=1.8.0   |
   |                   | framework                  |              |           |
   +-------------------+----------------------------+--------------+-----------+
   | ``gtest-devel``   | Development files for      | Compile-time | >=1.8.0   |
   |                   | gtest                      |              |           |
   +-------------------+----------------------------+--------------+-----------+
