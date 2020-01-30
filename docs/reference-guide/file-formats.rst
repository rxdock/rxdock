File formats
============

.prm file format
----------------

The .prm file format is an |Dock|-specific text format and is used for:

* system definition files (known previously as receptor ``.prm`` files)
* scoring function definition files
* search protocol definition files

The format is simple and allows for an arbitrary number of named parameter/value
pairs to be defined, optionally divided into named sections. Sections provide a
namespace for parameter names, to allow parameter names to be duplicated within
different sections. The key features of the format are:

* The first line of the file must be ``RBT_PARAMETER_FILE_V1.00`` with no
  preceeding whitespace.
* Subsequent lines may contain either:

  1. comment lines
  2. reserved keywords ``TITLE``, ``SECTION``, or ``END_SECTION``
  3. parameter name/value pairs

* Comment lines should start with a ``#`` character in the first column with no
  preceeding whitespace, and are ignored.
* The reserved words must start in the first column with no preceeding
  whitespace.
* The ``TITLE`` record should occur only once in the file and is used to provide
  a title string for display by various scripts such as run ``rbscreen.pl``. The
  keyword should be followed by a single space character and then the title
  string, which may contain spaces. If the ``TITLE`` line occurs more than once,
  the last occurence is used.
* ``SECTION`` records can occur more than once, and should always be paired with
  a closing ``END_SECTION`` record. The keyword should be followed by a single
  space character and then the section name, which may NOT itself contain
  spaces. All section names must be unique within a .prm file. All parameter
  name/value pairs within the ``SECTION`` / ``END_SECTION`` block belong to that
  section.
* Parameter name/value pairs are read as free-format tokenised text and can have
  preceeding, trailing, and be separated by arbitrary whitespace. This implies
  that the parameter name and value strings themselves are not allowed to
  contain any spaces. The value strings are interpreted as numeric, string, or
  boolean values as appropriate for that parameter. Boolean values should be
  entered as ``TRUE`` or ``FALSE`` uppercase strings.

.. caution::

   The current implementation of the .prm file reader does not tolerate a tab
   character immediately following the ``TITLE`` and ``SECTION`` keywords. It is
   very important that the first character after the ``SECTION`` keyword in
   particular is a true space character, otherwise the reserved word will not be
   detected and the parameters for that section will be ignored.

**Example .prm file** In the following example, ``RECEPTOR_FILE`` is defined in
the top level namespace. The remaining parameters are defined in the ``MAPPER``
and ``CAVITY`` namespaces. The indentation is for readability, and has no
significance in the format.

.. code-block:: python

   RBT_PARAMETER_FILE_V1.00
   TITLE 4dfr oxido-reductase

   RECEPTOR_FILE 4dfr.mol2

   SECTION MAPPER
      SITE_MAPPER RbtLigandSiteMapper
      REF_MOL 4dfr_c.sd
      RADIUS 6.0
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

.. _water-pdb-file-format:

Water PDB file format
---------------------

|Dock| requires explicit water PDB files to be in the style as output by the
`Dowser <https://www.ks.uiuc.edu/Research/vmd/plugins/dowser/>`__ program. In
particular:

* Records can be ``HETATM`` or ``ATOM``
* The atom names must be ``OW``, ``H1`` and ``H2``
* The atom records for each water molecule must belong to the same subunit ID
* The subunit IDs for different waters must be distinct, but do not have to be
  consecutive
* The atom IDs are not used and do not have to be consecutive (they can even be
  duplicated)
* The order of the atom records within a subunit is unimportant
* The temperature factor field of the water oxygens can be used to define the
  per-solvent flexibility modes. The temperature factors of the water hydrogens
  are not used.

.. _table-conversion-temperature-solvent-flexibility:

.. table:: Conversion of temperature values to solvent flexibility modes

   +-------------+---------------+--------------+
   | PDB         | Solvent       | Solvent      |
   | temperature | translational | rotational   |
   | factor      | flexibility   | flexibility  |
   +=============+===============+==============+
   | 0           | ``FIXED``     | ``FIXED``    |
   +-------------+---------------+--------------+
   | 1           | ``FIXED``     | ``TETHERED`` |
   +-------------+---------------+--------------+
   | 2           | ``FIXED``     | ``FREE``     |
   +-------------+---------------+--------------+
   | 3           | ``TETHERED``  | ``FIXED``    |
   +-------------+---------------+--------------+
   | 4           | ``TETHERED``  | ``TETHERED`` |
   +-------------+---------------+--------------+
   | 5           | ``TETHERED``  | ``FREE``     |
   +-------------+---------------+--------------+
   | 6           | ``FREE``      | ``FIXED``    |
   +-------------+---------------+--------------+
   | 7           | ``FREE``      | ``TETHERED`` |
   +-------------+---------------+--------------+
   | 8           | ``FREE``      | ``FREE``     |
   +-------------+---------------+--------------+

**Example** Valid |Dock| PDB file for explicit, flexible waters:

.. code-block:: python

   REMARK tmp 1YET.pdb xtal_hoh.pdb
   HETATM 3540 OW HOH W 106 28.929 12.684 20.864 1.00 1.0
   HETATM 3540 H1 HOH W 106 28.034 12.390 21.200 1.00
   HETATM 3540 H2 HOH W 106 29.139 12.204 20.012 1.00
   HETATM 3542 OW HOH W 108 27.127 14.068 22.571 1.00 2.0
   HETATM 3542 H1 HOH W 108 26.632 13.344 23.052 1.00
   HETATM 3542 H2 HOH W 108 27.636 13.673 21.806 1.00
   HETATM 3679 OW HOH W 245 27.208 10.345 27.250 1.00 3.0
   HETATM 3679 H1 HOH W 245 27.657 10.045 26.409 1.00
   HETATM 3679 H2 HOH W 245 26.296 10.693 27.036 1.00
   HETATM 3680 OW HOH W 246 31.737 12.425 21.110 1.00 4.0
   HETATM 3680 H1 HOH W 246 31.831 12.448 22.106 1.00
   HETATM 3680 H2 HOH W 246 30.775 12.535 20.863 1.00

Pharmacophore restraints file format
------------------------------------

Pharmacophore restraints are defined in a simple text file, with one restraint
per line. Each line should contain the following values, separated by commas or
whitespace:

.. code-block:: bash

   x y z coords of restraint centre, tolerance (in Angstroms), restraint type string

The supported restraint types are:

.. table:: Pharmacophore restraint types

   +---------+-----------------+-----------------------------------------------+
   | String  | Description     | Matches                                       |
   +=========+=================+===============================================+
   | ``Any`` | Any atom        | Any non-hydrogen atom                         |
   +---------+-----------------+-----------------------------------------------+
   | ``Don`` | H-bond donor    | Any neutral donor hydrogen                    |
   +---------+-----------------+-----------------------------------------------+
   | ``Acc`` | H-bond acceptor | Any neutral acceptor                          |
   +---------+-----------------+-----------------------------------------------+
   | ``Aro`` | Aromatic        | Any aromatic ring centre (pseudo atom)        |
   +---------+-----------------+-----------------------------------------------+
   | ``Hyd`` | Hydrophobic     | Any non-polar hydrogens (if present), any C   |
   |         |                 | sp3 or S sp3, any C or S not bonded to O sp2, |
   |         |                 | any Cl, Br, I                                 |
   +---------+-----------------+-----------------------------------------------+
   | ``Hal`` | Hydrophobic,    | Subset of Hyd, sp3 atoms only                 |
   |         | aliphatic       |                                               |
   +---------+-----------------+-----------------------------------------------+
   | ``Har`` | Hydrophobic,    | Subset of Hyd, aromatic atoms only            |
   |         | aromatic        |                                               |
   +---------+-----------------+-----------------------------------------------+
   | ``Ani`` | Anionic         | Any atom with negative distributed formal     |
   |         |                 | charge                                        |
   +---------+-----------------+-----------------------------------------------+
   | ``Cat`` | Cationic        | Any atom with positive distributed formal     |
   |         |                 | charge                                        |
   +---------+-----------------+-----------------------------------------------+
