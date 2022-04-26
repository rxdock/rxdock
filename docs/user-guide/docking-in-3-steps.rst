.. _docking-in-3-steps:

Docking in 3 steps
==================

You will find in this page a short tutorial for running |Dock|.

It has been divided in 3 steps:

1. System definition
2. Cavity generation
3. Docking

Step 1: System definition
-------------------------

First of all, we need to define the system.

Below these lines you have an example for a DUD system of a typical prm file (See :ref:`Documentation <reference-guide>` for more information):

.. code-block:: python

   RBT_PARAMETER_FILE_V1.00
   TITLE gart_DUD

   RECEPTOR_FILE gart_rdock.mol2
   RECEPTOR_FLEX 3.0

   ##############################################
   ## CAVITY DEFINITION: REFERENCE LIGAND METHOD
   ##############################################
   SECTION MAPPER
       SITE_MAPPER RbtLigandSiteMapper
       REF_MOL xtal-lig.sd
       RADIUS 6.0
       SMALL_SPHERE 1.0
       MIN_VOLUME 100
       MAX_CAVITIES 1
       VOL_INCR 0.0
   GRIDSTEP 0.5
   END_SECTION

   ############################
   ## CAVITY RESTRAINT PENALTY
   ############################
   SECTION CAVITY
       SCORING_FUNCTION RbtCavityGridSF
       WEIGHT 1.0
   END_SECTION

You will need this generated ``.prm`` file, a receptor structure mol2 file (``gart_rdock.mol2``) and a ligand file in the cavity (``xtal-lig.sd``) for going to next stage.

.. note::

   The receptor ``.mol2`` file must be preparated (protonated, charged, etc.) prior to this stage. The program chosen to do so is up to the user. As a suggestion, we usually work with MOE and/or Maestro.

Step 2: Cavity generation
-------------------------

Once the files are ready, a simple command will generate the cavity:

.. code-block:: bash

   rbcavity -was -d -r <PRMFILE>

With the ``-d`` flag a grid ``.grd`` file is generated. This file can be visualized in a molecular viewer to check the generated cavity.

For example, in PyMOL (after loading by: ``pymol <RECEPTOR>.mol2 <LIGAND>.sd <GRID>.grd``), write the following command in the console:

.. code-block:: python

   isomesh cavity, <GRID>.grd, 0.99

Step 3: Docking
---------------

Once the cavity is defined and generated, a 50 runs-per-ligand |Dock| job can be run straightforwardly with the following command:

.. note::

   The ``.prm`` file, receptor, reference ligand and ``.as`` cavity file must be in the working directory or pointed by the environmental variable ``RBT_HOME``.

.. code-block:: bash

   rbdock -i <INPUT>.sd -o <OUTPUT> -r <PRMFILE> -p dock.prm -n 50
