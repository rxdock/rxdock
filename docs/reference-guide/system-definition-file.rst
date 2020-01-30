.. _system-definition-file:

System definition file
======================

Although known previously as the receptor ``.prm`` file, the system definition
file has evolved to contain much more than the receptor information. The system
definition file is used to define:

* Receptor input files and flexibility parameters (the section called
  :ref:`Receptor definition <sdf-receptor>`)
* Explicit solvent input file and flexibility parameters (the section called
  :ref:`Solvent definition <sdf-solvent>`)
* Ligand flexibility parameters (the section called :ref:`Ligand definition
  <sdf-ligand>`).
* External restraint terms to be added to the total scoring function (e.g.
  :ref:`cavity restraint <sdf-cavity-restraint>`, :ref:`pharmacophoric restraint
  <sdf-pharmacophore-restraints>`)

.. _sdf-receptor:

Receptor definition
-------------------

The receptor can be loaded from a single MOL2 file, or from a combination of
Charmm PSF and CRD files. In the former case the MOL2 file provides the topology
and reference coordinates simultaneously, whereas in the latter case the
topology is loaded from the PSF file and the reference coordinates from the CRD
file. For historical compatibility reasons, receptor definition parameters are
all defined in the top-level namespace and should not be placed between
``SECTION`` / ``END_SECTION`` pairs.

.. caution::

   If MOL2 and PSF/CRD parameters are defined together, the MOL2 parameters take
   precedence and are used to load the receptor model.

.. table:: Receptor definition parameters

   +---------------------------+---------------------------------------------------------------------------+-------------------------+------------------+----------------------------------------------+
   | Parameter                 | Description                                                               | Type                    | Default          | Range of values                              |
   +===========================+===========================================================================+=========================+==================+==============================================+
   | **Parameters specific to loading receptor in MOL2 file format**                                                                                                                                   |
   +---------------------------+---------------------------------------------------------------------------+-------------------------+------------------+----------------------------------------------+
   | ``RECEPTOR_FILE``         | Name of receptor MOL2 file                                                | Filename string         | No default value | Valid MOL2 filename                          |
   +---------------------------+---------------------------------------------------------------------------+-------------------------+------------------+----------------------------------------------+
   | **Parameters specific to loading receptor in Charmm PSF/CRD file format**                                                                                                                         |
   +---------------------------+---------------------------------------------------------------------------+-------------------------+------------------+----------------------------------------------+
   | ``RECEPTOR_TOPOL_FILE``   | Name of receptor Charmm PSF file                                          | Filename string         | No default value | Valid Charmm PSF filename                    |
   +---------------------------+---------------------------------------------------------------------------+-------------------------+------------------+----------------------------------------------+
   | ``RECEPTOR_COORD_FILE``   | Name of receptor Charmm CRD file                                          | Filename string         | No default value | Valid Charmm CRD filename                    |
   +---------------------------+---------------------------------------------------------------------------+-------------------------+------------------+----------------------------------------------+
   | ``RECEPTOR_MASSES_FILE``  | Name of |Dock|-annotated Charmm masses file                               | Filename string         | No default value | ``masses.rtf``, ``top_all2_prot_na.inp``     |
   +---------------------------+---------------------------------------------------------------------------+-------------------------+------------------+----------------------------------------------+
   | **General receptor parameters, applicable to either file format**                                                                                                                                 |
   +---------------------------+---------------------------------------------------------------------------+-------------------------+------------------+----------------------------------------------+
   | ``RECEPTOR_SEGMENT_NAME`` | List of molecular segment names to read from either MOL2 or PSF/CRD file. | Comma separated list    | Empty            | Comma separated list of segment name strings |
   |                           | If this parameter is defined, then any segment/chains not listed          | of segment name strings | (i.e. all        |                                              |
   |                           | are not loaded. This provides a convenient way to remove cofactors,       | (without any spaces)    | segments read    |                                              |
   |                           | counterions and solvent without modifying the original file.              |                         | from file)       |                                              |
   +---------------------------+---------------------------------------------------------------------------+-------------------------+------------------+----------------------------------------------+
   | ``RECEPTOR_FLEX``         | Defines terminal OH and NH3+ groups withing this distance of docking      | float (Angstroms)       | Undefined        | >0.0 (3.0 is a reasonable value)             |
   |                           | volume as flexible.                                                       |                         | (rigid receptor) |                                              |
   +---------------------------+---------------------------------------------------------------------------+-------------------------+------------------+----------------------------------------------+
   | **Advanced parameters (should not need to be changed by the user)**                                                                                                                               |
   +---------------------------+---------------------------------------------------------------------------+-------------------------+------------------+----------------------------------------------+
   | ``RECEPTOR_ALL_H``        | Disable the removal of explicit non-polar hydrogens from the receptor     | boolean                 | FALSE            | TRUE or FALSE                                |
   |                           | model. *Not recommended*                                                  |                         |                  |                                              |
   +---------------------------+---------------------------------------------------------------------------+-------------------------+------------------+----------------------------------------------+
   | ``DIHEDRAL_STEP``         | Maximum mutation step size for receptor dihedral degrees of freedom.      | float (degrees)         | 30.0             | >0.0                                         |
   +---------------------------+---------------------------------------------------------------------------+-------------------------+------------------+----------------------------------------------+

.. _sdf-ligand:

Ligand definition
-----------------

Ligand definition parameters need only be defined if you wish to introduce
tethering of some or all of the ligand degrees of freedom. If you are running
conventional free docking then this section is not required. All ligand
definition parameters should be defined in ``SECTION LIGAND``. Note that the
ligand input SD file continues to be specified directly on the ``rbdock``
command-line and not in the system definition file.

.. table:: Ligand definition parameters

   +--------------------+--------------------------------------------------------------------+---------------------------+----------+--------------------------+
   | Parameter          | Description                                                        | Type                      | Default  | Range of values          |
   +====================+====================================================================+===========================+==========+==========================+
   | **Main user parameters**                                                                                                                                  |
   +--------------------+--------------------------------------------------------------------+---------------------------+----------+--------------------------+
   | ``TRANS_MODE``     | Sampling mode for ligand translational degrees of freedom          | enumerated string literal | ``FREE`` | ``FIXED``, ``TETHERED``, |
   |                    |                                                                    |                           |          | ``FREE``                 |
   +--------------------+--------------------------------------------------------------------+---------------------------+----------+--------------------------+
   | ``ROT_MODE``       | Sampling mode for ligand whole-body rotational degrees of freedom  | enumerated string literal | ``FREE`` | ``FIXED``, ``TETHERED``, |
   |                    |                                                                    |                           |          | ``FREE``                 |
   +--------------------+--------------------------------------------------------------------+---------------------------+----------+--------------------------+
   | ``DIHEDRAL_MODE``  | Sampling mode for ligand internal dihedral degrees of freedom      | enumerated string literal | ``FREE`` | ``FIXED``, ``TETHERED``, |
   |                    |                                                                    |                           |          | ``FREE``                 |
   +--------------------+--------------------------------------------------------------------+---------------------------+----------+--------------------------+
   | ``MAX_TRANS``      | (for ``TRANS_MODE = TETHERED`` only) Maximum deviation allowed     | float (Angstroms)         | 1.0      | >0.0                     |
   |                    | from reference centre of mass                                      |                           |          |                          |
   +--------------------+--------------------------------------------------------------------+---------------------------+----------+--------------------------+
   | ``MAX_ROT``        | (for ``ROT_MODE = TETHERED`` only) Maximum deviation allowed from  | float (degrees)           | 30.0     | >0.0--180.0              |
   |                    | orientation for reference principle axes                           |                           |          |                          |
   +--------------------+--------------------------------------------------------------------+---------------------------+----------+--------------------------+
   | ``MAX_DIHEDRAL``   | (for ``DIHEDRAL_MODE = TETHERED`` only) Maximum deviation allowed  | float (degrees)           | 30.0     | >0.0--180.0              |
   |                    | from reference dihedral angles for any rotatable bond              |                           |          |                          |
   +--------------------+--------------------------------------------------------------------+---------------------------+----------+--------------------------+
   | **Advanced parameters (should not need to be changed by the user)**                                                                                       |
   +--------------------+--------------------------------------------------------------------+---------------------------+----------+--------------------------+
   | ``TRANS_STEP``     | Maximum mutation step size for ligand translational degrees of     | float (Angstroms)         | 2.0      | >0.0                     |
   |                    | freedom                                                            |                           |          |                          |
   +--------------------+--------------------------------------------------------------------+---------------------------+----------+--------------------------+
   | ``ROT_STEP``       | Maximum mutation step size for ligand whole-body rotational        | float (degrees)           | 30.0     | >0.0                     |
   |                    | degrees of freedom                                                 |                           |          |                          |
   +--------------------+--------------------------------------------------------------------+---------------------------+----------+--------------------------+
   | ``DIHEDRAL_STEP``  | Maxium mutation step size for ligand internal dihedral degrees of  | float (degrees)           | 30.0     | >0.0                     |
   |                    | freedom                                                            |                           |          |                          |
   +--------------------+--------------------------------------------------------------------+---------------------------+----------+--------------------------+

.. _sdf-solvent:

Solvent definition
------------------

Solvent definition parameters need only be defined if you wish to introduce
explicit structural waters into the docking calculation, otherwise this section
is not required. All solvent definition parameters should be defined in
``SECTION SOLVENT``.

.. _table-solvent-definition-parameters:

.. table:: Solvent definition parameters

   +--------------------+--------------------------------------------------------------------+---------------------------+------------------+------------------+
   | Parameter          | Description                                                        | Type                      | Default          | Range of values  |
   +====================+====================================================================+===========================+==================+==================+
   | **Main user parameters**                                                                                                                                  |
   +--------------------+--------------------------------------------------------------------+---------------------------+------------------+------------------+
   | ``FILE``           | Name of explicit solvent PDB file                                  | File name string          | No default value | Valid PDB        |
   |                    |                                                                    |                           | (mandatory       | filename         |
   |                    |                                                                    |                           | parameter)       |                  |
   +--------------------+--------------------------------------------------------------------+---------------------------+------------------+------------------+
   | ``TRANS_MODE``     | Sampling mode for solvent translational degrees of freedom. If     | enumerated string literal | ``FREE``         | ``FIXED``,       |
   |                    | defined here, the value overrides the per-solvent translational    |                           |                  | ``TETHERED``,    |
   |                    | sampling modes defined in the solvent PDB file                     |                           |                  | ``FREE``         |
   +--------------------+--------------------------------------------------------------------+---------------------------+------------------+------------------+
   | ``ROT_MODE``       | Sampling mode for solvent whole-body rotational degrees of         | enumerated string literal | ``FREE``         | ``FIXED``,       |
   |                    | freedom. If defined here, the value overrides the per-solvent      |                           |                  | ``TETHERED``,    |
   |                    | rotational sampling modes defined in the solvent PDB file          |                           |                  | ``FREE``         |
   +--------------------+--------------------------------------------------------------------+---------------------------+------------------+------------------+
   | ``MAX_TRANS``      | (for ``TRANS_MODE = TETHERED`` waters only) Maximum deviation      | float (Angstroms)         | 1.0              | >0.0             |
   |                    | allowed from reference water oxygen positions. The same value is   |                           |                  |                  |
   |                    | applied to all waters with ``TRANS_MODE = TETHERED``; it is no     |                           |                  |                  |
   |                    | possible currently to define per-solvent ``MAX_TRANS`` values      |                           |                  |                  |
   +--------------------+--------------------------------------------------------------------+---------------------------+------------------+------------------+
   | ``MAX_ROT``        | (for ``ROT_MODE = TETHERED`` waters only) Maximum deviation        | float (degrees)           | 30.0             | >0.0--180.0      |
   |                    | allowed from orientation of reference principal axes. The same     |                           |                  |                  |
   |                    | value is applied to all waters with ``ROT_MODE = TETHERED``; it is |                           |                  |                  |
   |                    | no possible currently to define per-solvent ``MAX_ROT`` values     |                           |                  |                  |
   +--------------------+--------------------------------------------------------------------+---------------------------+------------------+------------------+
   | ``OCCUPANCY``      | Controls occupancy state sampling for all explicit solvent. If     | float                     | 1.0              | 0.0--1.0         |
   |                    | defined here, the values overrides the per-solvent occupancy       |                           |                  |                  |
   |                    | states defined in the solvent PDB file                             |                           |                  |                  |
   +--------------------+--------------------------------------------------------------------+---------------------------+------------------+------------------+
   | **Advanced parameters (should not need to be changed by the user)**                                                                                       |
   +--------------------+--------------------------------------------------------------------+---------------------------+------------------+------------------+
   | ``TRANS_STEP``     | Maximum mutation step size for solvent translational degrees of    | float (Angstroms)         | 2.0              | >0.0             |
   |                    | freedom                                                            |                           |                  |                  |
   +--------------------+--------------------------------------------------------------------+---------------------------+------------------+------------------+
   | ``ROT_STEP``       | Maximum mutation step size for solvent wholebody rotational        | float (degrees)           | 30.0             | >0.0             |
   |                    | degrees of freedom                                                 |                           |                  |                  |
   +--------------------+--------------------------------------------------------------------+---------------------------+------------------+------------------+
   | ``OCCUPANCY_STEP`` | Maximum mutation step size for solvent occupancy state degrees     | float (degrees)           | 1.0              | 0.0--1.0         |
   |                    | of freedom                                                         |                           |                  |                  |
   +--------------------+--------------------------------------------------------------------+---------------------------+------------------+------------------+

**Solvent occupancy state sampling** ``OCCUPANCY = 0`` permanently disables all
solvent; ``OCCUPANCY = 1.0`` permanently enables all solvent; ``OCCUPANCY``
between 0 and 1 activates variable occupancy state sampling, where the value
represents the initial probability that the solvent molecule is enabled. For
example, ``OCCUPANCY = 0.5`` means that the solvent is enabled in 50 % of the
initial GA population. However, the probability that the solvent is actually
enabled in the final docking solution will depend on the particular ligand, the
scoring function terms, and on the penalty for solvent binding. The occupancy
state chromosome value is managed as a continuous variable between 0.0 and 1.0,
with a nominal mutation step size of 1.0. Chromosome values lower than the
occupancy threshold (defined as ``1.0 - OCCUPANCY``) result in the solvent being
disabled; values higher than the threshold result in the solvent being enabled.

.. _sdf-cavity-mapping:

Cavity mapping
--------------

The cavity mapping section is mandatory. You should choose one of the mapping
algorithms shown below. All mapping parameters should be defined in
``SECTION MAPPER``.

.. table:: Two sphere site mapping parameters

   +------------------+-----------------------------------------+--------------------------------+-------------------------+-----------------------------------+
   | Parameter        | Description                             | Type                           | Default                 | Range of values                   |
   +==================+=========================================+================================+=========================+===================================+
   | **Main user parameters**                                                                                                                                  |
   +------------------+-----------------------------------------+--------------------------------+-------------------------+-----------------------------------+
   | ``SITE_MAPPER``  | Mapping algorithm specifier             | string literal                 | ``RbtSphereSiteMapper`` | fixed                             |
   +------------------+-----------------------------------------+--------------------------------+-------------------------+-----------------------------------+
   | ``CENTER``       | (x,y,z) center of cavity mapping region | Bracketed cartesian coordinate | None                    | None                              |
   |                  |                                         | string (x,y,z)                 |                         |                                   |
   +------------------+-----------------------------------------+--------------------------------+-------------------------+-----------------------------------+
   | ``RADIUS``       | Radius of cavity mapping region         | float (Angstroms)              | 10.0                    | >0.0 (10.0--20.0 suggested range) |
   +------------------+-----------------------------------------+--------------------------------+-------------------------+-----------------------------------+
   | ``SMALL_SPHERE`` | Radius of small probe                   | float (Angstroms)              | 1.5                     | >0.0 (1.0--2.0 suggested range)   |
   +------------------+-----------------------------------------+--------------------------------+-------------------------+-----------------------------------+
   | ``LARGE_SPHERE`` | Radius of large probe                   | float (Angstroms)              | 4.0                     | >\ ``SMALL_SPHERE`` (3.5--6.0     |
   |                  |                                         |                                |                         | suggested range)                  |
   +------------------+-----------------------------------------+--------------------------------+-------------------------+-----------------------------------+
   | ``MAX_CAVITIES`` | Maximum number of cavities to accept    | integer                        | 99                      | >0                                |
   |                  | (in descending order of size)           |                                |                         |                                   |
   +------------------+-----------------------------------------+--------------------------------+-------------------------+-----------------------------------+
   | **Advanced parameters (less frequently changed by the user)**                                                                                             |
   +------------------+-----------------------------------------+--------------------------------+-------------------------+-----------------------------------+
   | ``VOL_INCR``     | Receptor atom radius increment for      | float (Angstroms)              | 0.0                     | >=0.0                             |
   |                  | excluded volume                         |                                |                         |                                   |
   +------------------+-----------------------------------------+--------------------------------+-------------------------+-----------------------------------+
   | ``GRID_STEP``    | Grid resolution for mapping             | float (Angstroms)              | 0.5                     | >0.0 (0.3--0.8 suggested range)   |
   +------------------+-----------------------------------------+--------------------------------+-------------------------+-----------------------------------+
   | ``MIN_VOLUME``   | Minimum cavity volume to accept         | float (Angstroms\ :sup:`3`)    | 100                     | >0 (100--300 suggested range)     |
   |                  | (in Å\ :sup:`3`, not grid points)       |                                |                         |                                   |
   +------------------+-----------------------------------------+--------------------------------+-------------------------+-----------------------------------+

.. table:: Reference ligand site mapping parameters

   +------------------+-----------------------------------------+--------------------------------+-------------------------+-----------------------------------+
   | Parameter        | Description                             | Type                           | Default                 | Range of values                   |
   +==================+=========================================+================================+=========================+===================================+
   | **Main user parameters**                                                                                                                                  |
   +------------------+-----------------------------------------+--------------------------------+-------------------------+-----------------------------------+
   | ``SITE_MAPPER``  | Mapping algorithm specifier             | string literal                 | ``RbtLigandSiteMapper`` | fixed                             |
   +------------------+-----------------------------------------+--------------------------------+-------------------------+-----------------------------------+
   | ``REF_MOL``      | Reference ligand SD file name           | string                         | ``ref.sd``              | None                              |
   +------------------+-----------------------------------------+--------------------------------+-------------------------+-----------------------------------+
   | ``RADIUS``       | Radius of cavity mapping region         | float (Angstroms)              | 10.0                    | >0.0 (10.0--20.0 suggested range) |
   +------------------+-----------------------------------------+--------------------------------+-------------------------+-----------------------------------+
   | ``SMALL_SPHERE`` | Radius of small probe                   | float (Angstroms)              | 1.5                     | >0.0 (1.0--2.0 suggested range)   |
   +------------------+-----------------------------------------+--------------------------------+-------------------------+-----------------------------------+
   | ``LARGE_SPHERE`` | Radius of large probe                   | float (Angstroms)              | 4.0                     | >\ ``SMALL_SPHERE`` (3.5--6.0     |
   |                  |                                         |                                |                         | suggested range)                  |
   +------------------+-----------------------------------------+--------------------------------+-------------------------+-----------------------------------+
   | ``MAX_CAVITIES`` | Maximum number of cavities to accept    | integer                        | 99                      | >0                                |
   |                  | (in descending order of size)           |                                |                         |                                   |
   +------------------+-----------------------------------------+--------------------------------+-------------------------+-----------------------------------+
   | **Advanced parameters (less frequently changed by the user)**                                                                                             |
   +------------------+-----------------------------------------+--------------------------------+-------------------------+-----------------------------------+
   | ``VOL_INCR``     | Receptor atom radius increment for      | float (Angstroms)              | 0.0                     | >=0.0                             |
   |                  | excluded volume                         |                                |                         |                                   |
   +------------------+-----------------------------------------+--------------------------------+-------------------------+-----------------------------------+
   | ``GRID_STEP``    | Grid resolution for mapping             | float (Angstroms)              | 0.5                     | >0.0 (0.3--0.8 suggested range)   |
   +------------------+-----------------------------------------+--------------------------------+-------------------------+-----------------------------------+
   | ``MIN_VOLUME``   | Minimum cavity volume to accept         | float (Å\ :sup:`3`)            | 100                     | >0 (100--300 suggested range)     |
   |                  | (in Å\ :sup:`3`, not grid points)       |                                |                         |                                   |
   +------------------+-----------------------------------------+--------------------------------+-------------------------+-----------------------------------+

.. _sdf-cavity-restraint:

Cavity restraint
----------------

The cavity restraint penalty function is mandatory and is designed to prevent
the ligand from exiting the docking site. The function is calculated over all
non-hydrogen atoms in the ligand (and over all explicit water oxygens that can
translate). The distance from each atom to the nearest cavity grid point is
calculated. If the distance exceeds the value of ``RMAX``, a penalty is imposed
based on the value of (``distance - RMAX``). The penalty can be either linear or
quadratic depending on the value of the ``QUADRATIC`` parameter. It should not
be necessary to change any the parameters in this section. Note that the docking
protocol itself will manipulate the ``WEIGHT`` parameter, so any changes made to
``WEIGHT`` will have no effect.

.. code-block:: python

   SECTION CAVITY
      SCORING_FUNCTION RbtCavityGridSF
      WEIGHT 1.0
      RMAX 0.1
      QUADRATIC FALSE
   END_SECTION

.. _sdf-pharmacophore-restraints:

Pharmacophore restraints
------------------------

This section need only be defined if you wish to dock with pharmacophore
restraints. If you are running conventional free docking then this section is
not required. All pharmacophore definition parameters should be defined in
``SECTION PHARMA``.

.. table:: Pharmacophore restraint parameters

   +----------------------+---------------------------------------------------------+------------------+------------+------------------------------------------+
   | Parameter            | Description                                             | Type             | Default    | Range of values                          |
   +======================+=========================================================+==================+============+==========================================+
   | ``CONSTRAINTS_FILE`` | Mandatory pharmacophore restraints file                 | File name string | None       | Valid file name                          |
   |                      |                                                         |                  | (mandatory |                                          |
   |                      |                                                         |                  | parameter) |                                          |
   +----------------------+---------------------------------------------------------+------------------+------------+------------------------------------------+
   | ``OPTIONAL_FILE``    | Optional pharmacophore restraints file                  | File name string | None       | Valid file name, or empty                |
   |                      |                                                         |                  | (optional  |                                          |
   |                      |                                                         |                  | parameter) |                                          |
   +----------------------+---------------------------------------------------------+------------------+------------+------------------------------------------+
   | ``NOPT``             | Number of optional restraints that should be met        | Integer          | 0          | Between 0 and number of restraints in    |
   |                      |                                                         |                  |            | ``OPTIONAL_FILE``                        |
   +----------------------+---------------------------------------------------------+------------------+------------+------------------------------------------+
   | ``WRITE_ERRORS``     | Ligands with insufficient pharmacophore features to     | Boolean          | FALSE      | TRUE of FALSE                            |
   |                      | match the mandatory restraints are always removed prior |                  |            |                                          |
   |                      | to docking. If this parameter is true, the pre-filtered |                  |            |                                          |
   |                      | ligands are written to an error SD file with the same   |                  |            |                                          |
   |                      | root name as the docked pose output SD file, but with   |                  |            |                                          |
   |                      | an ``_errors.sd`` suffix. If false, the pre-filtered    |                  |            |                                          |
   |                      | ligands are not written.                                |                  |            |                                          |
   +----------------------+---------------------------------------------------------+------------------+------------+------------------------------------------+
   | ``WEIGHT``           | Overall weight for the pharmacophore penalty function   | Float            | 1.0        | >=0.0                                    |
   +----------------------+---------------------------------------------------------+------------------+------------+------------------------------------------+

**Calculation of mandatory restraint penalty** The list of ligand atoms that
matches each restraint type in the mandatory restraints file is precalculated
for each ligand as it is loaded. If the ligand contains insufficient features
to satisfy all of the mandatory restraints the ligand is rejected and is not
docked. Note that the rejection is based purely on feature counts and does not
take into account the possible geometric arrangements of the features. Rejected
ligands are optionally written to an error SD file. The penalty for each
restraint is based on the distance from the nearest matching ligand atom to the
pharmacophore restraint centre. If the distance is less than the defined
tolerance (restraint sphere radius), the penalty is zero. If the distance is
greater than the defined tolerance a quadratic penalty is applied, equal to
(nearest distance - tolerance)\ :sup:`2`.

**Calculation of optional restraint penalty** The individual restraint penalties
for each restraint in the optional restraints file are calculated in the same
way as for the mandatory penalties. However, only the ``NOPT`` lowest scoring
(least penalised) restraints are summed for any given docking pose. Any
remaining higher scoring optional restraints are ignored and do not contribute
to the total pharmacophore restraint penalty.

**Calculation of overall restraint penalty** The overall pharmacophore restraint
penalty is the sum of the mandatory restraint penalties and the ``NOPT`` lowest
scoring optional restraint penalties, multiplied by the ``WEIGHT`` parameter
value.

NMR restraints
--------------

To be completed. However, this feature has rarely been used.

Example system definition files
-------------------------------

Full system definition file with all sections and common parameters enumerated explicitly:

.. code-block:: python

   RBT PARAMETER_FILE_V1.00
   TITLE HSP90-PU3-lig-cavity, solvent flex=5
   RECEPTOR_FILE PROT_W3_flex.mol2
   RECEPTOR_SEGMENT_NAME PROT
   RECEPTOR_FLEX 3.0
   SECTION SOLVENT
      FILE PROT_W3_flex_5.pdb
      TRANS_MODE TETHERED
      ROT_MODE TETHERED
      MAX_TRANS 1.0
      MAX_ROT 30.0
      OCCUPANCY 0.5
   END_SECTION
   SECTION_LIGAND
      TRANS_MODE FREE
      ROT_MODE FREE
      DIHEDRAL_MODE FREE
      MAX_TRANS 1.0
      MAX_ROT 30.0
      MAX_DIHEDRAL 30.0
   END_SECTION
   SECTION MAPPER
      SITE_MAPPER RbtLigandSiteMapper
      REF_MOL ref.sd
      RADIUS 5.0
      SMALL_SPHERE 1.0
      MIN_VOLUME 100
      MAX_CAVITIES 1
      VOL_INCR 0.0
      GRIDSTEP 0.5
   END_SECTION
   SECTION CAVITY
      SCORING_FUNCTION RbtCavityGridSF
      WEIGHT 1.0
   END_SECTION
   SECTION PHARMA
      SCORING_FUNCTION RbtPharmaSF
      WEIGHT 1.0
      CONSTRAINTS_FILE mandatory.const
      OPTIONAL FILE optional.const
      NOPT 3
      WRITE_ERRORS TRUE
   END_SECTION
