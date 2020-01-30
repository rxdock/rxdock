.. _docking-strategies:

Docking strategies
==================

This section does not pretend to be a comprehensive user guide. It does,
however, highlight the key steps the user must take for different docking
strategies, and may serve as a useful checklist in writing such a guide
in the future.

Standard docking
----------------

By standard docking, we refer to docking of a flexible, untethered ligand to a
receptor in the absence of explicit structural waters or any experimental
restraints.

Standard docking workflow
^^^^^^^^^^^^^^^^^^^^^^^^^

1. Prepare a MOL2 file for the protein or nucleic acid target, taking into
   account the atom typing issues described above for MOL2 file parsing. The
   recommendation is to prepare an all-atom MOL2 file and allow |Dock| to remove
   the non-polar hydrogens on-the-fly.

   .. important::

      Make sure that any non-standard atom names and substructure names are
      defined in ``$RBT_ROOT/data/sf/RbtIonicAtoms.prm`` in order for the
      assignment of distributed formal charges to work correctly. Make sure that
      the Tripos atom types are set correctly. |Dock| uses the Tripos types to
      derive other critical atomic properties such as atomic number and
      hybridisation state.

   .. note::

      The |Dock| MOL2 parser was developed to read the CCDC/Astex protein
      ``.mol2`` files, therefore this validation set is the de facto standard
      reference. You should compare against the format of the CCDC/Astex MOL2
      files if you are in doubt as to whether a particular MOL2 file is suitable
      for |Dock|.

2. Prepare a system definition file. At a minimum, you need to define the
   receptor parameters, the cavity mapping parameters (``SECTION MAPPER``) and
   the cavity restraint penalty (``SECTION CAVITY``). Make sure you define the
   ``RECEPTOR_FLEX`` parameter if you wish to activate sampling of terminal OH
   and NH3+ groups in the vicinity of the docking site.

3. Generate the docking site (.as) file using ``rbcavity``. You will require a
   reference bound ligand structure in the coordinate space of the receptor if
   you wish to use the reference ligand cavity mapping method.

4. Prepare the ligand SD files you wish to dock, taking into account the atom
   typing issues described above for SD file parsing. In particular, make sure
   that formal charges and formal bond order are defined coherently so that
   there are no formal valence errors in the file. |Dock| will report any
   perceived valence errors but will dock the structures anyway. Note that |Dock|
   never samples bond lengths, bond angles, ring conformations, or non-rotatable
   bonds during docking so initial conformations should be reasonable.

5. Run a small test calculation to check that the system is defined correctly.
   For example, run ``rbdock`` from the command line with a small ligand SD
   file, with the score-only protocol (``-p score.prm``) and with the ``-T 2``
   option to generate verbose output. The output will include receptor atom
   properties, ligand atom properties, flexibility parameters, scoring function
   parameters and docking protocol parameters.

6. When satisfied, launch the full-scale calculations. A description of the
   various means of launching |Dock| is beyond the scope of this guide.

Tethered scaffold docking
-------------------------

In tethered scaffold docking, the ligand poses are restricted and forced to
overlay the substructure coordinates of a reference ligand. The procedure is
largely as for standard docking, except that:

* Ligand SD files must be prepared with the ``rbtether`` utility to annotate
  each record with the matching substructure atom indices, and to transform
  the coordinates of each ligand so that the matching substructure coordinates
  are overlaid with the reference substructure coordinates. This requires a
  Daylight SMARTS toolkit license.

* The system definition file should contain a ``SECTION LIGAND`` to define which
  of the the ligand degrees of freedom should be tethering to their reference
  values. Tethering can be applied to position, orientation and dihedral degrees
  of freedom independently. Note that the tethers are applied directly within
  the chromosome representation used by the search engine (where they affect the
  randomisation and mutation operators), and therefore external restraint
  penalty functions to enforce the tethers are not required.

.. important::

   The reference state values for each tethered degree of freedom are defined
   directly from the initial conformation of each ligand as read from the input
   SD file, and not from the reference SD file used by ``rbtether``. This is
   why the ligand coordinates are transformed by ``rbtether``, such that each
   ligand record can act as its own reference state. The reference SD file used
   by ``rbtether`` is not referred to by the docking calculation itself.

It follows from the above that tethered ligand docking is inappropriate for
input ligand SD files that have not already been transformed to the coordinate
space of the docking site, either by ``rbtether`` or by some other means.

Example ligand definition for tethered scaffold
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This definition will tether the position and orientation of the tethered
substructure, but will allow free sampling of ligand dihedrals.

.. code-block:: python

   SECTION LIGAND
      TRANS_MODE TETHERED
      ROT_MODE TETHERED
      DIHEDRAL_MODE FREE
      MAX_TRANS 1.0
      MAX_ROT 30.0
   END_SECTION

Docking with pharmacophore restraints
-------------------------------------

In pharmacophore restrained docking, ligand poses are biased to fit user-defined
pharmacophore points. The bias is introduced through the use of an external
penalty restraint, which penalises docking poses that do not match the
pharmacophore restraints. Unlike tethered scaffold docking, there is no
modification to the chromosome operators themselves, hence the search can be
inefficient, particularly for large numbers of restraints and/or for ligands
with large numbers of matching features. Pre-screening of ligands is based
purely on feature counts, and not on geometric match considerations.

The implementation supports both mandatory and optional pharmacophore
restraints. The penalty function is calculated over all mandatory restraints,
and over (any ``NOPT`` from N) of the optional restraints. For example, you may
wish to ensure that any 4 from 7 optional restraints are satisfied in the
generated poses.

The procedure is largely as for standard docking, except that:

* You should prepare separate pharmacophore restraint files for the mandatory
  and optional restraints. Note that optional restraints do not have to be
  defined, it is sufficient to only define at least one mandatory restraint.
* The system definition file should contain a ``SECTION PHARMA`` to add the
  pharmacophore restraint penalty to the scoring function.

Docking with explicit waters
----------------------------

Explicit structural waters can be loaded from an external PDB file,
independently from the main receptor model, by adding a ``SECTION SOLVENT``
to the system definition file. The user has fine control over the flexibility of
each water molecule. A total of 9 flexibility modes are possible, in which the
translational and rotational degrees of freedom of each water can be set
independently to ``FIXED``, ``TETHERED``, or ``FREE``. Thus, for example, it is
possible to define a water with a fixed oxygen coordinate (presumably at a
crystallographically observed position), but freely rotating such that the
orientation of the water hydrogens can be optimised by the search engine (and
can be ligand-dependent).

.. note::

   In the current implementation, solvent refers strictly to water molecules,
   and the format of the water PDB file is very strictly defined. In future
   implementations it is anticipated that other, larger (and possibly flexible)
   molecules will be loadable as solvent, and that other file formats will be
   supported.

Explicit waters workflow
^^^^^^^^^^^^^^^^^^^^^^^^

1. Prepare a separate PDB file for the explicit waters according to the format
   prescribed (the section called :ref:`Water PDB file format
   <water-pdb-file-format>`).

2. Add a ``SECTION SOLVENT`` to the system definition file and define the
   relevant flexibility parameters (:numref:`Table %s
   <table-solvent-definition-parameters>`). The minimal requirement is to define
   the ``FILE`` parameter.

3. Decide whether you wish to have different per-solvent flexibility modes
   (defined via the occupancy values and temperature factor values in the PDB
   file (:numref:`Table %s
   <table-conversion-temperature-solvent-flexibility>`)), or whether you wish to
   have a single flexibility mode applied to all waters (defined via the
   ``TRANS_MODE`` and ``ROT_MODE`` values in the ``SECTION SOLVENT`` of the
   receptor .prm file).

   .. important::
   
      If you wish to use per-solvent flexibility modes (that is, you wish to set
      different modes for different waters) make sure that you do not define
      ``TRANS_MODE`` or ``ROT_MODE`` entries in the ``SECTION SOLVENT`` as these
      values will override the per-solvent values derived from the temperature
      factors in the PDB file.

4. If you have defined any waters with ``TETHERED`` translational or rotational
   degrees of freedom, define ``MAX_TRANS`` and/or ``MAX_ROT`` values as
   appropriate (or accept the default values). The tethered ranges are applied
   to all tethered waters and can not be defined on a per-solvent basis at
   present.
