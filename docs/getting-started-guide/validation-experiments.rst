.. _validation-experiments:

Validation experiments
======================

In this section you will find the instructions about how to reproduce our
validation experiments using different test sets. Three different sets were
analyzed for three different purposes:

* ASTEX set for binding mode prediction in proteins.
* RNA set for assess RNA-ligand docking.
* DUD set for database enrichment.

Binding mode prediction in proteins
-----------------------------------

First of all, please go to `ASTEX set SourceForge download page
<https://sourceforge.net/projects/rdock/files/Validation_Sets/ASTEX_rDock_TestSet.tar.gz/download>`__
to download a compressed file with the necessary data.

After downloading the file ``ASTEX_rDock_TestSet.tar.gz``, uncompress the file
with the following command, which will create a folder called
``ASTEX_rDock_TestSet``:

.. code-block:: bash

   $ tar -xvzf ASTEX_rDock_TestSet.tar.gz
   $ cd ASTEX_rDock_TestSet/

Here you will have the instructions for one of the systems (``1sj0``), to run
with the rest of the systems, just change the pdb code with the one desired.
Then, make sure that the necessary environmental variables for running |Dock| are
well defined and run the following commands for entering to the folder and
running |Dock| with the same settings that we have used:

.. code-block:: bash

   $ cd 1sj0/

   #first create the cavity using rbcavity
   $ rbcavity -r 1sj0_rdock.prm -was > 1sj0_cavity.log

   #then use rbdock to run docking
   $ rbdock -r 1sj0_rdock.prm -p dock.prm -n 100 -i 1sj0_ligand.sd \
   -o 1sj0_docking_out > 1sj0_docking_out.log

   #sdsort for sorting the results according to their score
   $ sdsort -n -f'SCORE' 1sj0_docking_out.sd > 1sj0_docking_out_sorted.sd

   #calculate rmsd from the output comparing with the crystal structure of the ligand
   $ sdrmsd 1sj0_ligand.sd 1sj0_docking_out_sorted.sd

Binding mode prediction in RNA
------------------------------

In a similar way of the section above, here you will find a brief tutorial on
how to run |Dock| with the RNA TestSet used in the validation. As in the first
section, please go to `RNA set SourceForge download page
<https://sourceforge.net/projects/rdock/files/Validation_Sets/RNA_rDock_TestSet.tar.gz/download>`__
to download a compressed file with the necessary data.

After downloading the file ``RNA_rDock_TestSet.tar.gz``, uncompress the file
with the following command, which will create a folder called
``RNA_rDock_TestSet``:

.. code-block:: bash

   $ tar -xvzf RNA_rDock_TestSet.tar.gz
   $ cd RNA_rDock_TestSet/

Here you will have the instructions for one of the systems (``1nem``), to run
with the rest of the systems, just change the pdb code with the one desired.
Then, make sure that the necessary environmental variables for running |Dock| are
well defined and run the following commands for entering to the folder and
running |Dock| with the same settings that we have used (if you have run the
previous set, the variables should already be correctly defined):

.. code-block:: bash

   $ cd 1nem/

   #first create the cavity using rbcavity
   $ rbcavity -r 1nem_rdock.prm -was > 1nem_cavity.log

   #then use rbdock to run docking
   $ rbdock -r 1nem_rdock.prm -p dock.prm -n 100 -i 1nem_lig.sd \
   -o 1nem_docking_out > 1nem_docking_out.log

   #sdsort for sorting the results according to their score
   $ sdsort -n -f'SCORE' 1nem_docking_out.sd > 1nem_docking_out_sorted.sd

   #calculate rmsd from the output comparing with the crystal structure of the ligand
   $ sdrmsd 1nem_lig.sd 1nem_docking_out_sorted.sd

Database enrichment (actives vs. decoys -- for HTVS)
----------------------------------------------------

In this section you will find a brief tutorial on how to run |Dock| with the DUD
test set used in the validation and how to perform different analysis of the
results. As in the sections above, please go to `DUD set SourceForge download
page <https://sourceforge.net/projects/rdock/files/Validation_Sets/DUD_rDock_TestSet.tar.gz/download>`__
to download a compressed file with the necessary data.

After downloading the file ``DUD_rDock_TestSet.tar.gz``, uncompress the file
with the following command, which will create a folder called
``DUD_rDock_TestSet``:

.. code-block:: bash

   $ tar -xvzf DUD_rDock_TestSet.tar.gz
   $ cd DUD_rDock_TestSet/

Here you will have the instructions for one of the systems (``hivpr``), to run
with the rest of the systems, just change the DUD system code with the one
desired. Then, make sure that the necessary environmental variables for running
|Dock| are well defined and run the following commands for entering to the folder
and running |Dock| with the same settings that we have used (if you have run the
previous sets, the variables should already be correctly defined):

.. code-block:: bash

   $ cd hivpr/

   #first create the cavity using rbcavity
   $ rbcavity -r hivpr_rdock.prm -was > hivpr_cavity.log

As the number of ligands to dock is very high, we suggest you to use any
distributed computing environments, such as SGE or Condor, and configure |Dock|
to run in multiple CPUs. Namely, split the input ligands file in as many parts
as desired (very easy using ``sdsplit`` tool) and run independent |Dock| docking
jobs for each "splitted" input file. However, for this example purpose, you will
have the instructions for running all set of actives and decoys in one docking
job:

.. code-block:: bash

   #uncompress ligand file
   $ gunzip hivpr_ligprep.sdf.gz

   #use rbdock to run docking
   $ rbdock -r hivpr_rdock.prm -p dock.prm -n 100 -i hivpr_ligprep.sdf \
   -o hivpr_docking_out > hivpr_docking_out.log

   #sdsort with -n and -s flags will sort internally each ligand by increasing
   #score and sdfilter will get only the first entry of each ligand
   $ sdsort -n -s -fSCORE hivpr_docking_out.sd | sdfilter \
   -f'$_COUNT == 1' > hivpr_1poseperlig.sd

   #sdreport will print all the scores of the output in a tabular format and,
   #with command awk, we will format the results
   $ sdreport -t hivpr_1poseperlig.sd | awk '{print $2,$3,$4,$5,$6,$7}' > dataforR_uq.txt

At this point, you should have a file called ``hivpr_docking_out.sd`` with all
docking poses written by |Dock| (100 * number of input ligands), a file called
``hivpr_1poseperlig.sd`` with the best scored docking pose for each ligand and a
file called ``dataforR_uq.txt`` that will be used for calculating ROC curves
using R. The next step is to calculate ROC curves and other statistics. To do
so, please visit section :ref:`Calculating ROC curves <calculating-roc-curves>`
and jump to the subsection "R commands for generating ROC curves".
