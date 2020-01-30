Pharmacophoric restraints
=========================

In this short tutorial you will find how to prepare and run Docking with
pharmacophoric restraints.

.. note::

   |Dock| assumes the user knows how to compute and find pharmacophores. The user
   will need the coordinates, tolerance and type of restraint, which will be the
   input for |Dock|.

Step 1: Pharmacophoric restraints file
--------------------------------------

The first step is to create the input file for |Dock| with the necessary info.

As you can find in the :ref:`Reference guide <reference-guide>`, this file needs
one line per pharmacophore with the next structure (each element separated with
a space):

.. code-block:: bash

   x y z (coords of restraint centre), tolerance radius (in Angstroms), restraint type (string)

The pharmacophore types accepted by |Dock| are the following:

+--------+-------------------------+-------------------------------------------+
| String | Description             | Matches                                   |
+========+=========================+===========================================+
| Any    | Any atom                | Any non-hydrogen atom                     |
+--------+-------------------------+-------------------------------------------+
| Don    | H-Bond donor            | Any neutral donor hydrogen                |
+--------+-------------------------+-------------------------------------------+
| Acc    | H-Bond acceptor         | Any neutral acceptor                      |
+--------+-------------------------+-------------------------------------------+
| Aro    | Aromatic                | Any aromatic ring centre (pseudo atom)    |
+--------+-------------------------+-------------------------------------------+
| Hyd    | Hydrophobic             | Any C or S sp3, any C od S not bonded     |
|        |                         | to O sp2, any Cl, Br, I                   |
+--------+-------------------------+-------------------------------------------+
| Hal    | Hydrophobic, aliphatic  | Subset of Hyd, sp3 atoms only             |
+--------+-------------------------+-------------------------------------------+
| Har    | Hydrophobic, aromatic   | Subset of Hyd, aromatic atoms only        |
+--------+-------------------------+-------------------------------------------+
| Ani    | Anionic                 | Any atom with negative distributed        |
|        |                         | formal charge                             |
+--------+-------------------------+-------------------------------------------+
| Cat    | Cationic                | Any atom with positive distributed        |
|        |                         | formal charge                             |
+--------+-------------------------+-------------------------------------------+

A sample file (``pharma.restr``) has been created containing two restraints
(``Acc`` and ``Hyd``) with a tolerance radius of 2 and located at points (-1.75,
1.25, 0.25) and (-2, 2, -3) respectively.

.. code-block:: python

   -1.75 1.25 0.25 2.0 Acc
   -2.00 2.00 -3.0 2.0 Hyd

Step 2: |Dock| system definition file
-------------------------------------

The second and final step is to modify the system definition file (``FILE.prm``)
to take into account the defined restraints.

Just add a ``SECTION PHARMA`` (see the :ref:`Reference guide <reference-guide>`
for more info) with the following lines:

.. code-block:: python

   SECTION PHARMA
       SCORING_FUNCTION RbtPharmaSF
       WEIGHT 1.0
       CONSTRAINTS_FILE pharma.restr
   END_SECTION

With the ``FILE.prm`` finally being:

.. code-block:: python
   :emphasize-lines: 32-36

   RBT_PARAMETER_FILE_V1.00
   TITLE title

   RECEPTOR_FILE receptor_file.mol2
   RECEPTOR_FLEX 3.0

   ##############################################
   ## CAVITY DEFINITION: REFERENCE LIGAND METHOD
   ##############################################
   SECTION MAPPER
       SITE_MAPPER RbtLigandSiteMapper
       REF_MOL reference.sdf
       RADIUS 5.0
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

   #############################
   ## PHARMACOPHORIC RESTRAINTS
   #############################
   SECTION PHARMA
       SCORING_FUNCTION RbtPharmaSF
       WEIGHT 1.0
       CONSTRAINTS_FILE pharma.restr
   END_SECTION

.. note::

   This ``FILE.prm`` is an example file for the sake of the tutorial. The point
   here is to clarify how to define pharmacophoric restraints and how to
   configure |Dock| to take them into account.

Finally, when running |Dock|, the user can check if the program writes similar
lines as the following to be sure that the restraints have been correctly
read.

.. code-block:: bash

   [...]
   RbtPharmaSF: Reading mandatory ph4 constraints from /path/to/pharma.restr
   (-1.75, 1.25, 0.25)  2.0     Acc
   (-2.0, 2.0, -3.0)    2.0     Hyd
   RbtPharmaSF: No optional ph4 constraints file found
   [...]

Step 3: Optional constraints
----------------------------

This tutorial is an example for mandatory constraints. Optional constraints can
also be configured in a different file (same format as pharma.restr created
above). The ``SECTION PHARMA`` in the |Dock| System Configuration File should be
modified as follows (NOPT=1 means that only one of the optional restraints has
to be met):

.. code-block:: python

   SECTION PHARMA
       SCORING_FUNCTION RbtPharmaSF
       WEIGHT 1.0
       CONSTRAINTS_FILE pharma.restr
       OPTIONAL_FILE optional_pharma.restr
       NOPT 1
   END_SECTION

.. tip::

   For more information about the pharmacophoric restraints and the parameters
   in ``SECTION PHARMA``, please go to the :ref:`Reference guide
   <reference-guide>`.
