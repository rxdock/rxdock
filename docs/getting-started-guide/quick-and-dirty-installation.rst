.. _quick-and-dirty-installation:

Quick and dirty installation
============================

In this section you will have short instructions to make a typical installation
of |Dock|.

To get the full documentation of all |Dock| software package and methods, please
go to the :ref:`Reference guide <reference-guide>`.

Moreover, you can also check the following information:

* :ref:`Getting started <getting-started-guide>`: installation and validation
  instructions for first-time users.
* :ref:`Validation experiments <validation-experiments>`: instructions and
  examples for re-running the validation sets we have carried out.
* :ref:`Calculating ROC curves <calculating-roc-curves>`: tutorial for
  generating ROC curves and other statistics after running |Dock| docking jobs.

Installation in 3 steps
-----------------------

We have been able to compile |Dock| in the following Linux systems:

* CentOS 5.5 64 bits
* openSUSE 11.3 32 and 64 bits
* openSUSE 12.3 32 and 64 bits
* openSUSE 13.1 32 and 64 bits
* Ubuntu 12.04 32 and 64 bits

Step 1
^^^^^^

First of all, you will need to install several packages before compiling and
running |Dock|:

* gcc and g++ compilers version > 3.3
* make
* cppunit and cppunit-devel

.. note::

   **For Ubuntu users:**

   If you are trying to use |Dock| in Ubuntu, please note that csh shell is not
   included in a default installation. We recommend to install csh in case some
   error arises (``sudo apt-get install csh``), even with all the above-stated
   dependencies installed.

Afterwards, download the source code compressed file or get it by SVN in
:ref:`Downloads <download>` section.

Step 2
^^^^^^

Then, run the following commands:

.. code-block:: bash

   $ tar -xvzf rxdock-0.1.0.tar.gz
   $ cd rxdock-0.1.0/build/

and, for 32 bits computers:

.. code-block:: bash

   $ make linux-g++

for 64 bits computers:

.. code-block:: bash

   $ make linux-g++-64

Step 3
^^^^^^

After compiling successfully, type the following command to make a test and
check that your compiled version works good and the results are correct.

.. code-block:: bash

   $ make test

If the test has succeed, you are done, enjoy using |Dock|!

Otherwise, please check your dependencies and all the previous commands or go to
:ref:`Support Section <support>` to ask for help.

Just as a concluding remark, don't forget to set the necessary environmental
variables for running |Dock| in the command line (for example, in Bash shell):

.. code-block:: bash

   $ export RBT_ROOT=/path/to/rxdock/installation/
   $ export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$RBT_ROOT/lib
   $ export PATH=$PATH:$RBT_ROOT/bin
