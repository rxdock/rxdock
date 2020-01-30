Running docking jobs in parallel
================================

In this short tutorial we will try to explain how to run |Dock| on a computer with multiple CPUs or a cluster with different calculation nodes.

.. note::

   |Dock| has not an MPI version to be run in parallel on a computation cluster. The approach |Dock| uses to parallelize the jobs is rather simple: as each molecule can be run in an independent way, the input structure file is splitted in multiple files and each of them is run independently.

For this example, we have a set of **200 molecules** (``input.sdf``) and we want to run it in **10 CPUs**.

Step 1: Split molecules input file
----------------------------------

To split an SDF file (|Dock| needs the input in SDF format), there is a script in |Dock| package called ``sdsplit`` that does this.

.. code-block:: bash

   $ sdsplit
   Splits SD records into multiple files of equal size

   Usage:  sdsplit [-<RecSize>] [-o<OutputRoot>] [sdFiles]

           -<RecSize>      record size to split into (default = 1000 records)
           -o<OutputRoot>  Root name for output files (default = tmp)

           If SD file list not given, reads from standard input

In our case, to split 200 molecules in 10 files (with 20 molecules each), we will have to run the following command that will generate 10 files called ``split[1-10].sd``:

.. code-block:: bash

   sdsplit -20 -osplit input.sdf

Moreover, you can use the following code which allows you to specify the number of files you want instead of the number of molecules in each file (e.g., save it in a file named ``splitMols.sh``):

.. code-block:: bash

   #!/bin/bash
   #Usage: splitMols.sh <input> #Nfiles <outputRoot>
   fname=$1
   nfiles=$2
   output=$3
   molnum=$(grep -c '$$$$' $fname)
   echo "$molnum molecules found"
   echo "Dividing '$fname' into $nfiles files"
   rawstep=`echo $molnum/$nfiles | bc -l`
   let step=$molnum/$nfiles
   if [ ! `echo $rawstep%1 | bc` == 0 ]; then
           let step=$step+1;
   fi;
   sdsplit -$step -o$output $1

To get the same as in the first case, run:

.. code-block:: bash

   splitMols.sh input.sdf 10 split

Step 2: Run docking jobs with splitted files
--------------------------------------------

We have two options:

* Run docking jobs locally: send it over 10 CPUs.
* Run docking jobs using a job scheduler.

Option 1: Run docking jobs locally
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To run |Dock| (standard mode, 50 runs per ligand) in 10 CPUs, be sure that all the necessary files are located in the working directory: receptor mol2 file, prm file, cavity as file, and reference ligand for cavity definition (if used) and run the following command:

.. code-block:: bash

   for file in split*sd; do rbdock -i $file -o ${file%%.*}_out -r <PRMFILE> -p dock.prm -n 50 &; done

This will send 10 independent docking jobs and will eventually generate 10 output files ``split[1-10]_out.sd``.

So that's it, you are done!

Option 2: Run docking jobs with job scheduler
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Same as in Option 1, but instead of running the command above, you have to create a queueing submission file for each of the files and submit them to the queue.

There are several options to use as a job scheduler. In our particular case, we use SGE and a typical submission file looks like this:

.. code-block:: bash

   #!/bin/sh
   #$ -N rdock_job1
   #$ -S /bin/sh
   #$ -q serial
   #$ -o out.log
   #$ -e err.log
   #-cwd
   export RBT_ROOT=/data/soft/rdock/2006.1
   export LD_LIBRARY_PATH=$RBT_ROOT/lib
   #next is optional
   export RBT_HOME=/path/to/job/files 

   # These are the comands to be executed.
   cd /path/to/job/files
   $RBT_ROOT/bin/rbdock -i <INPUT>.sd -o <OUTPUT> -r <PRMFILE> -p dock.prm -n 50

This is highly recommended for running docking jobs of big molecule libraries.

For example, to run a Virtual Screening Campaign of a million compounds, you can split the molecules in 10000 files in order to have individual files with 100 molecules each and use a job scheduler to control their execution.
