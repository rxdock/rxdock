Molecular files and atom typing
===============================

Macromolecular targets (protein or RNA) are input from Tripos MOL2 files
(``RbtMOL2FileSource`` class) or from pairs of Charmm PSF (``RbtPsfFileSource``
class) and CRD (``RbtCrdFileSource`` class) files. Ligands are input from MDL
Information Systems (MDL) structure data (SD) files (``RbtMdlFileSource``
class). Explicit structural waters are input optionally from PDB files
(``RbtPdbFileSource`` class). Ligand docking poses are output to MDL SD files.

The |Dock| scoring functions have been defined and validated for implicit
non-polar hydrogen (extended carbon) models only. If you provide all-atom
models, be aware that the non-polar hydrogens will be removed automatically.
Polar hydrogens must be defined explicitly in the molecular files, and are not
added by |Dock|. Positive ionisable and negative ionisable groups can be
automatically protonated and deprotonated respectively to create common charged
groups such as guanidinium and carboxylic acid groups.

MOL2 is now the preferred file format for |Dock| as it eliminates many of the
atom typing issues inherent in preparing and loading PSF files. The use of
PSF/CRD files is strongly discouraged. The recommendation is to prepare an
all-atom MOL2 file with correct Tripos atom types assigned, and allow |Dock| to
remove non-polar hydrogens on-the-fly.

Atomic properties
-----------------

|Dock| requires the following properties to be defined per atom. Depending on the
file format, these properties may be loaded directly from the molecular input
file, or derived internally once the model is loaded:

* Cartesian (x,y,z) coordinates
* Element (atomic number)
* Formal hybridisation state (sp, sp2, sp3, aromatic, trigonal planar)
* Formal charge
* Distributed formal charge (known informally as group charge)
* Tripos force field type (|Dock| uses a modified version of the Sybyl 5.2 types,
  extended to include carbon types with implicit non-polar hydrogens)
* Atom name
* Substructure (residue) name
* Atomic radius (assigned per element from ``$RBT_ROOT/data/RbtElements.dat``)

.. note::

   The |Dock| scoring functions do not use partial charges and therefore partial
   charges do not have to be defined. The atomic radii are simplified radii
   defined per element, and are used for cavity mapping and in the polar scoring
   function term, but are not used in the vdW scoring function term. The latter
   has its own indepedent parameterisation based on the Tripos force field
   types.

Difference between formal charge and distributed formal charge
--------------------------------------------------------------

The formal charge on an atom is always an integer. For example, a charged
carboxylic acid group (COO-) can be defined formally as a formal double bond to
a neutral oxygen sp2, and a formal single bond to a formally charged oxygen sp3.
In reality of course, both oxygens are equivalent. |Dock| distributes the integer
formal charge across all equivalent atoms in the charged group that are
topologically equivalent. In negatively charged acid groups, the formal charge
is distributed equally between the acid oxygens. In positively charged amines,
the formal charge is distributed equally between the hydrogens. In charged
guanidinium, amidinium, and imidazole groups, the central carbon also receives
an equal portion of the formal charge (in addition to the hydrogens). The
distributed formal charge is also known as the group charge. The polar scoring
functions in |Dock| use the distributed formal charge to scale the polar
interaction strength of the polar interactions.

Parsing a MOL2 file
-------------------

``MOLECULE``, ``ATOM``, ``BOND``, and ``SUBSTRUCTURE`` records are parsed. The
atom name, substructure name, Cartesian coordinates and Tripos atom type are
read directly for each atom. The element type (atomic number) and formal
hybridisation state are derived from the Tripos type using an internal lookup
table. Formal charges are not read from the MOL2 file and do not have to be
assigned correctly in the file. Distributed formal charges are assigned directly
by |Dock| based on standard substructure and atom names as described below.

Parsing an SD file
------------------

Cartesian coordinates, element and formal charge are read directly for each
atom. Formal bond orders are read for each bond. Atom names are derived from
element name and atom ID (e.g. C1, N2, C3). The substructure name is ``MOL``.
Formal hybridisation states are derived internally for each atom based on
connectivity patterns and formal bond orders. The Tripos types are asssigned
using internal rules based on atomic number, formal hybridisation state and
formal charges. The integer formal charges are distributed automatically
across all topologically equivalent atoms in the charged group.

Assigning distributed formal charges to the receptor
----------------------------------------------------

|Dock| provides a file format independent method for assigning distributed formal
charges directly to the receptor atoms, which is used by the MOL2 and PSF/CRD
file readers. The method uses a lookup table based on standard substructure and
atom names, and does not require the integer formal charges to be assigned to
operate correctly.

The lookup table file is ``$RBT_ROOT/data/sf/RbtIonicAtoms.prm``. Each section
name represents a substructure name that contains formally charged atoms. The
entries within the section represent the atom names and distributed formal
charges for that substructure name. The file provided with |Dock| contains
entries for all standard amino acids and nucleic acids, common metals, and
specific entries required for processing the GOLD CCDC/Astex validation sets.

.. important::

   You may have to extend ``RbtIonicAtoms.prm`` if you are working with
   non-standard receptor substructure names and/or atom names, in order for the
   distributed formal charges to be assigned correctly.
