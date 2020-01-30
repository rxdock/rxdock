Multi-step protocol for HTVS
============================

For high-throughput virtual screening (HTVS) applications, where computing
performance is important, the recommended |Dock| protocol is to limit the search
space (i.e. rigid receptor), apply the grid-based scoring function and/or to use
a multi-step protocol to stop sampling of poor scorers as soon as possible.

Using a multi-step protocol for the DUD system COMT, the computational time can
be reduced by 7.5-fold without affecting performance by:

1. Running 5 docking runs for all ligands;
2. ligands achieving a score of -22 or lower run 10 further runs;
3. for those ligands achieving a score of -25 or lower, continue up to 50 runs.

The optimal protocol is specific for each particular system and parameter-set,
but can be identified with a purpose-built script (see the :ref:`Reference guide
<reference-guide>`, section ``rbhtfinder``).

Here you will find a tutorial to show you how to create and run a multi-step
protocol for a HTVS campaign.

Step 1: Create the multi-step protocol
--------------------------------------

These are the instructions for running rbhtfinder:

.. code-block:: bash

   1st) exhaustive docking of a small representative part of the
        whole library.
   2nd) Store the result of sdreport -t over that exhaustive dock.
        in file  that will be the input of this
        script.
   3rd) rbhtfinder <sdreport_file> <output_file> <thr1max> <thr1min> <ns1> <ns2>
        <ns1> and <ns2> are the number of steps in stage 1 and in
        stage 2. If not present, the default values are 5 and 15
        <thrmax> and <thrmin> setup the range of thresholds that will
        be simulated in stage 1. The threshold of stage 2 depends
        on the value of the threshold of stage 1.
        An input of -22 -24 will try protocols:
                5   -22     15      -27
                5   -22     15      -28
                5   -22     15      -29
                5   -23     15      -28
                5   -23     15      -29
                5   -23     15      -30
                5   -24     15      -29
                5   -24     15      -30
                5   -24     15      -31
        Output of the program is a 7 column values. First column
        represents the time. This is a percentage of the time it
        would take to do the docking in exhaustive mode, i.e. 
        docking each ligand 100 times. Anything
        above 12 is too long.
        Second column is the first percentage. Percentage of
        ligands that pass the first stage.
        Third column is the second percentage. Percentage of
        ligands that pass the second stage.
        The four last columns represent the protocol.
        All the protocols tried are written at the end.
        The ones for which time is less than 12%, perc1 is
        less than 30% and perc2 is less than 5% but bigger than 1%
        will have a series of *** after, to indicate they are good choices
        WARNING! This is a simulation based in a small set.
        The numbers are an indication, not factual values.

Step 1, substep 1: Exhaustive docking
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Hence, as stated, the first step is to run an **exhaustive docking** of a
representative part of the whole desired library to dock.

For |Dock|, exhaustive docking means doing **100 runs** for each ligand, whereas
standard docking means 50 runs for each ligand:

.. code-block:: bash

   $ rbdock -i INPUT.sd -o OUTPUT -r PRMFILE.prm -p dock.prm -n 100

Step 1, substep 2: ``sdreport`` summary
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Once the exhaustive docking has finished, the results have to be saved in a
**single file** and the output of the script ``sdreport -t`` will be used as
**input for** ``rbhtfinder``:

.. code-block:: bash

   $ sdreport -t OUTPUT.sd > sdreport_results.txt

Step 1, substep 3: ``rbhtfinder`` script
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The **last step** is to run the ``rbhtfinder`` script (download
:download:`sdreport_results.txt <_downloads/sdreport_results.txt.zip>` for
testing):

.. code-block:: bash

   $ rbhtfinder sdreport_results.txt htvs_protocol.txt -10 -20 7 25

Which will result in a file called ``htvs_protocol.txt``.

The parameters are explained in the script instructions. They are not always the
same and as they depend on the system, you will probably have to play a little
with different values in order to **obtain good parameters sets** (marked with
``***`` in the output).

This will happen when **time** is less than 12%, **perc1** (number of ligands
that pass the first filter) is less than 30% and **perc2** (number of ligands
that pass the second filter) is less than 5% but bigger than 1%.

Step 2: Run docking with the multi-step protocol
------------------------------------------------

The script finished with two good parameters sets:
::

   TIME PERC1 PERC2 N1 THR1 N2 THR2
   [...]
   11.928, 27.461, 3.207, 7, -12, 25, -17 ***
   [...]
   10.508, 18.773, 1.511, 7, -13, 25, -18 ***
   [...]

These parameters have to be adapted to a **file** with the HTVS **protocol
format** that |Dock| understands.

A **template file** looks as follows (``THR1``, ``THR2``, ``N1`` and ``N2``
are the parameters found above):
::

   3
   if - <THR1> SCORE.INTER 1.0 if - SCORE.NRUNS <N1-1> 0.0 -1.0,
   if - <THR2> SCORE.INTER 1.0 if - SCORE.NRUNS <N2-1> 0.0 -1.0,
   if - SCORE.NRUNS 49 0.0 -1.0,
   1
   - SCORE.INTER -10,

It is divided in 2 sections, **Running Filters** and **Writing Filters**
(defined by the lines with one number).

The first line (the number 3) indicates the number of lines in the Running
Filters:

* The first filter is defined as follows: if the number of runs reaches ``N1``
  and the score is lower than ``THR1``, continue to filter 2, else stop with
  that ligand and go to the next one.
* The second filter is defined similar to the first one: if the number of runs
  reaches ``N2`` and the score is lower than ``THR2``, continue to filter 3,
  else stop and go to the next ligand.
* If a ligand has passed the first two filters, continue up to 50 runs.

The fifth line (the number 1 after the three Running Filters) indicates the
number of lines in the Writing Filters:

* Only print out all those poses where ``SCORE.INTER`` is lower than -10 (for
  avoiding excessive printing).

For the parameters obtained in the first Section of this tutorial (first line
with ``***``), we will have to generate a file as follows:
::

   3
   if - -12 SCORE.INTER 1.0 if - SCORE.NRUNS 6 0.0 -1.0,
   if - -17 SCORE.INTER 1.0 if - SCORE.NRUNS 24 0.0 -1.0,
   if - SCORE.NRUNS 49 0.0 -1.0,
   1
   - SCORE.INTER -10,

Please note that the parameters ``N1`` and ``N2`` are 7 and 25 but we write 6
and 24, respectively, as stated in the template.

Finally, **run** |Dock| changing the flag ``-n XX`` for ``-t PROTOCOLFILE.txt``:

.. code-block:: bash

   $ rbdock -i INPUT.sd -o OUTPUT -r PRMFILE.prm -p dock.prm -t PROTOCOLFILE.txt
