.. _programs:

Programs
========

Programs summary tables:

.. table:: Core |Dock| C++ executables

   +-------------------+-------------+---------------------------------------------------------------------------------+
   | Executable        | Used for    | Description                                                                     |
   +===================+=============+=================================================================================+
   | ``rbcavity``      | Preparation | Cavity mapping and preparation of docking site (.as) file.                      |
   +-------------------+-------------+---------------------------------------------------------------------------------+
   | ``rbcalcgrid``    | Preparation | Calculation of vdW grid files (usually called by ``make_grid.csh`` wrapper      |
   |                   |             | script).                                                                        |
   +-------------------+-------------+---------------------------------------------------------------------------------+
   | ``rbdock``        | Docking     | The main |Dock| docking engine itself.                                          |
   +-------------------+-------------+---------------------------------------------------------------------------------+

.. table:: Auxiliary |Dock| programs

   +-------------------+-------------+---------------------------------------------------------------------------------+
   | Executable        | Used for    | Description                                                                     |
   +===================+=============+=================================================================================+
   | ``sdtether``      | Preparation | Prepares a ligand SD file for tethered scaffold docking. Annotates ligand SD    |
   |                   |             | file with tethered substructure atom indices. Requires Open Babel Python        |
   |                   |             | bindings.                                                                       |
   +-------------------+-------------+---------------------------------------------------------------------------------+
   | ``rbhtfinder``    | Preparation | Used to optimise a high-throughput docking protocol from an initial exhaustive  |
   |                   |             | docking of a small representative ligand library. Parametrize a multi-step      |
   |                   |             | protocol for your system.                                                       |
   +-------------------+-------------+---------------------------------------------------------------------------------+
   | ``make_grid.csh`` | Preparation | Creates the vdW grid files required for grid-based docking protocols            |
   |                   |             | (``dock_grid.prm`` and ``dock_solv_grid.prm``). Simple front-end to             |
   |                   |             | ``rbcalcgrid``.                                                                 |
   +-------------------+-------------+---------------------------------------------------------------------------------+
   | ``rbconvgrid``    | Analysis    | Converts |Dock| vdW grids to InsightII grid format for visualisation.           |
   +-------------------+-------------+---------------------------------------------------------------------------------+
   | ``rbmoegrid``     | Analysis    | Converts |Dock| vdW grids to MOE grid format for visualisation.                 |
   +-------------------+-------------+---------------------------------------------------------------------------------+
   | ``rblist``        | Analysis    | Outputs miscellaneous information for ligand SD file records.                   |
   +-------------------+-------------+---------------------------------------------------------------------------------+
   | ``sdrmsd``        | Analysis    | Calculation of ligand Root Mean Squared Displacement (RMSD) between reference   |
   |                   |             | and docked poses, taking into account ligand topological symmetry. Requires     |
   |                   |             | Open Babel Python bindings.                                                     |
   +-------------------+-------------+---------------------------------------------------------------------------------+
   | ``sdfilter``      | Analysis    | Utility for filtering SD files by arbitrary data field expressions. Useful for  |
   |                   |             | simple post-docking filtering by score components.                              |
   +-------------------+-------------+---------------------------------------------------------------------------------+
   | ``sdsort``        | Analysis    | Utility for sorting SD files by arbitrary data field. Useful for simple         |
   |                   |             | post-docking filtering by score components.                                     |
   +-------------------+-------------+---------------------------------------------------------------------------------+
   | ``sdreport``      | Analysis    | Utility for reporting SD file data field values. Output in tab-delimited or CSV |
   |                   |             | format.                                                                         |
   +-------------------+-------------+---------------------------------------------------------------------------------+
   | ``sdsplit``       | Utility     | Splits an SD file into multiple smaller SD files of fixed number of records.    |
   +-------------------+-------------+---------------------------------------------------------------------------------+
   | ``sdmodify``      | Utility     | Sets the molecule title line of each SD record equal to a given SD data field.  |
   +-------------------+-------------+---------------------------------------------------------------------------------+

Programs reference
------------------

rbdock
^^^^^^

``rbdock`` -- The |Dock| docking engine itself.

.. code-block:: bash

   $RBT_ROOT/bin/rbdock
   {-i input ligand MDL SD file}
   {-o output MDL SD file}
   {-r system definition.prm file}
   {-p docking protocol.prm file}
   [-n number of docking runs/ligand]
   [-s random seed]
   [-T debug trace level]
   [[-t SCORE.INTER threshold] | [-t filter definition file]]
   [-ap -an -allH -cont]

Simple exhaustive docking
"""""""""""""""""""""""""

The minimum requirement for rbdock is to specify the input (``-i``) and output
(``-o``) ligand SD file names, the system definition .prm file (``-r``) and the
docking protocol .prm file (``-p``). This will perform one docking run per
ligand record in the input SD file and output all docked ligand poses to the
output SD file. Use ``-n`` to increase the number of docking runs per ligand
record.

High-throughput docking, option 1
"""""""""""""""""""""""""""""""""

The ``-t`` and ``-cont`` options can be used to construct high-throughput
protocols. If the argument following ``-t`` is numeric it is interpreted as a
threshold value for ``SCORE.INTER``, the total intermolecular score between
ligand and receptor/solvent. In the absence of -cont, the threshold acts as an
early termination filter, and the docking runs for each ligand will be
terminated early once the threshold value has been exceeded. Note that the
threshold is applied only at the end of each individual docking run, not during
the runs themselves. If the ``-cont`` (continue) option is specified as well,
the threshold acts as an output pose filter instead of a termination filter.
The docking runs for each ligand run to completion as in the exhaustive case,
but only the docking poses that exceed the threshold value of ``SCORE.INTER``
are written to the output SD file.

High throughput docking, option 2
"""""""""""""""""""""""""""""""""

Alternatively, if the argument following ``-t`` is non-numeric it is interpreted
as a filter definition file. The filter definition file can be used to define
multiple termination filters and multiple output pose filters in a generic way.
Any docking score component can be used in the filter definitions.
``run_rbscreen.pl`` generates a filter definition file for multi-stage,
high-throughput docking, with progressive score thresholds for early termination
of poorly performing ligands. The use of filter definition files is preferred
over the more limited ``SCORE.INTER`` filtering described above, whose use is
now deprecated.

Automated ligand protonation/deprotonation
""""""""""""""""""""""""""""""""""""""""""

The ``-ap`` option activates the automated protonation of ligand positive
ionisable centres, notably amines, guanidines, imidazoles, and amidines. The
``-an`` option activates the automated deprotonation of ligand negative
ionisable centres, notably carboxylic acids, phosphates, phosphonates,
sulphates, and sulphonates. The precise rules used by |Dock| for protonation
and deprotonation are quite crude, and are not user-customisable. Therefore
these flags are not recommended for detailed validation experiments, in which
care should be taken that the ligand protonation states are set correctly in the
input SD file. Note that |Dock| is not capable of converting ionised centres back
to the neutral form; these are unidirectional transformations.

Control of ligand non-polar hydrogens
"""""""""""""""""""""""""""""""""""""

By default, |Dock| uses an implicit non-polar hydrogen model for receptor and
ligand, and all of the scoring function validation has been performed on this
basis. If the ``-allH`` option is not defined (recommended), all explicit
non-polar hydrogens encountered in the ligand input SD file are removed, and
only the polar hydrogens (bonded to O, N, or S) are retained. If the ``-allH``
option is defined (not recommended), no hydrogens are removed from the ligand.
Note that |Dock| is not capable of adding explicit non-polar hydrogens, if none
exist. In other words, the ``-allH`` option disables hydrogen removal, it does
not activate hydrogen addition. You should always make sure that polar hydrogens
are defined explicitly. If the ligand input SD file contains no explicit
non-polar hydrogens, the ``-allH`` option has no effect. Receptor protonation is
controlled by the system definition prm file.

rbcavity
^^^^^^^^

``rbcavity`` -- Cavity mapping and preparation of docking site (``.as``) file.

.. code-block:: bash

   $RBT_ROOT/bin/rbcavity
   {-r system definition .prm file}
   [-ras -was -d -v -s]
   [-l distance from cavity]
   [-b border]

Exploration of cavity mapping parameters
""""""""""""""""""""""""""""""""""""""""

``rbcavity -r .prm file``

You can run rbcavity with just the ``-r`` argument when first preparing a new
receptor for docking. This allows you to explore rapidly the impact of the
cavity mapping parameters on the generated cavities, whilst avoiding the
overhead of actually writing the docking site (``.as``) file to disk. The number
of cavities and volume of each cavity are written to standard output.

Visualisation of cavities
"""""""""""""""""""""""""

``rbcavity -r .prm file -d``

If you have access to InsightII you can use the ``-d`` option to dump the cavity
volumes in InsightII grid file format. There is no need to write the docking
site (``.as``) file first. The InsightII grid files should be loaded into the
reference coordinate space of the receptor and contoured at a contour level of
0.99.

Writing the docking site (.as) file
"""""""""""""""""""""""""""""""""""

``rbcavity -r .prm file -was``

When you are happy the mapping parameters, use the ``-was`` option to write the
docking site (``.as``) file to disk. The docking site file is a binary file that
contains the cavity volumes in a compact format, and a pre-calculated cuboid
grid extending over the cavities. The grid represents the distance from each
point in space to the nearest cavity grid point, and is used by the cavity
penalty scoring function. Calculating the distance grid can take a long time
(whereas the cavity mapping itself is usually very fast), hence the ``-was``
option should be used sparingly.

Analysis of cavity atoms
""""""""""""""""""""""""

``rbcavity -r .prm file -ras -l distance``

Use the ``-l`` options to list the receptor atoms within a given distance of any
of the cavity volumes, for example to determine which receptor OH/NH3+ groups
should be flexible. This option requires access to the pre-calculated distance
grid embedded within the docking site (.as) file, and is best used in
combination with the ``-ras`` option, which loads a previously generated docking
site file. This avoids the time consuming step of generating the cavity distance
grid again. If ``-l`` is used without ``-ras``, the cavity distance grid will be
calculated on-the-fly each time.

Miscellaneous options
"""""""""""""""""""""

The ``-s`` option writes out various statistics on the cavity and on the
receptor atoms in the vicinity of the cavity. These values have been used in
genetic programming model building for docking pose false positive removal. The
``-v`` option writes out the receptor coordinates in PSF/CRD format for use by
the rDock Viewer (not documented here). Note that the PSF/CRD files are not
suitable for simulation purposes, only for visualisation, as the atom types are
not set correctly. The ``-b`` option controls the size of the cavity distance
grid, and represents the border beyond the actual cavity volumes. It should not
be necessary to vary this parameter (default = 8 Å) unless longer-range scoring
functions are implemented.

rbcalcgrid
^^^^^^^^^^

``rbcalcgrid`` -- Calculation of vdW grid files (usually called by
``make_grid.csh`` wrapper script).

.. code-block:: bash

   $RBT_ROOT/bin/rbcalcgrid
   {-r system definition file}
   {-o output suffix for generated grids}
   {-p vdW scoring function prm file}
   [-g grid step]
   [-b border]

Note that, unlike ``rbdock`` and ``rbcavity``, spaces are not tolerated between
the command-line options and their corresponding arguments. See
``$RBT_ROOT/bin/make_grid.csh`` for common usage.

make_grid.csh
^^^^^^^^^^^^^

Creates vdW grids for all receptor ``.prm`` files listed on command line.
Front-end to ``rbcalcgrid``.

rbconvgrid
^^^^^^^^^^

rbmoegrid
^^^^^^^^^

``rbmoegrid`` -- Calculates grids for a given atom type.

.. code-block:: bash

   rbmoegrid -o <OutputRoot> -r <ReceptorPrmFile> -p <SFPrmFile> [-g <GridStep> -b <border> -t <tripos_type>]

   -o <OutFileName> (.grd is suffiexed)
   -r <ReceptorPrmFile> - receptor param file (contains active site params)
   -p <SFPrmFile> - scoring function param file (default calcgrid_vdw.prm)
   -g <GridStep> - grid step (default = 0.5A)
   -b <Border> - grid border around docking site (default = 1.0A)
   -t <AtomType> - Tripos atom type (default is C.3)

sdrmsd
^^^^^^

``sdrmsd`` -- calculation of ligand root mean squared displacement (RMSD)
between reference and docked poses. It takes into account molecule topological
symmetry. Requires Open Babel Python bindings.

.. code-block:: bash

   $RBT_ROOT/bin/sdrmsd [options] {reference SD file} {input SD file}

With two arguments
""""""""""""""""""

``sdrmsd`` calculates the RMSD between each record in the input SD file and the
first record of the reference SD file. If there is a mismatch in the number of
atoms, the record is skipped and the RMSD is not calculated. The RMSD is
calculated over the heavy (non-hydrogen) atoms only. Results are output to
standard output. If some record was skipped, a warning message will be printed
to standard error.

With fitting
""""""""""""

A molecular superposition will be done before calculation of the RMSD. The
output will specify an RMSD FIT calculation was done.

.. code-block:: bash

   sdrmsd -o output.sdf reference.sdf input.sdf
   sdrmsd --out=output.sdf reference.sdf input.sdf

Output a SD file
""""""""""""""""

This option will write an output SD file with the input molecules adding an
extra RMSD field to the file. If fitting was done, the molecule coordinates
will also be fitted to the reference.

.. code-block:: bash

   sdrmsd -o output.sdf reference.sdf input.sdf
   sdrmsd --out=output.sdf reference.sdf input.sdf

sdtether
^^^^^^^^

``sdtether`` -- Prepares a ligand SD file for tethered scaffold docking.
Requires Open Babel Python bindings. Annotates ligand SD file with tethered
substructure atom indices.

.. code-block:: bash

   $RBT_ROOT/bin/sdtether {ref. SDfile} {in SDfile} {out SDfile} "{SMARTS query}"

``sdtether`` performs the following actions:

* Runs the SMARTS query against the reference SD file to determine the tethered
  substructure atom indices and coordinates.
* If more than one substructure match is retrieved (e.g. due to topological
  symmetry, or if the query is too simple) all substructure matchs are retained
  as the reference and all ligands will be tethered according to all possible
  matches.
* Runs the SMARTS query against each record of the input ligand SD file in turn.
* For each substructure match, the ligand coordinates are transformed such that
  the principal axes of the matching substructure coordinates are aligned with
  the reference substructure coordinates.
* In addition, an SD data field is added to the ligand record which lists the
  atom indices of the substructure match, for later retrieval by |Dock|.
* Each transformed ligand is written to the output SD file.
* Note that if the SMARTS query returns more than one substructure match for a
  ligand, that ligand is written multiple times to the output file, once for
  each match, each of which will be docked independently with different
  tethering information.

sdfilter
^^^^^^^^

``sdfilter`` -- Post-process an SD file by filtering the records according to
data fields or attributes.
::

   sdfilter -f '$<DataField> <Operator> <Value>' [-s <DataField>] [sdFiles]

or ::

   sdfilter -f <filename> [-s <DataField>] [sdFiles]

.. note::

   Multiple filters are allowed and are OR'd together. Filters can be provided
   in a file, one per line. Standard Perl operators should be used. e.g.

   .. code-block:: perl

      eq ne lt gt le ge # for strings
      == != < > <= >=   # for numeric

   ``_REC`` (record #) is provided as a pseudo-data field. If ``-s`` option is used, ``_COUNT`` (#occurrences of DataField) is provided as a pseudo-data field. If SD file list not given, reads from standard input. Output is to standard output.

For example, if ``results.sd`` constains multiple ligands each having multiple poses (ordered by score), then running
::

   sdfilter -f'$_COUNT == 1' results.sd

will get you the first entry for each ligand.

sdreport
^^^^^^^^

``sdreport`` -- Produces text summaries of SD records.
::

   sdreport [-l] [-t [<FieldName, FieldName...>]] [-c <FieldName, FieldName...>] [-id <IDField>] [-nh] [-o] [-s] [-sup] [sdFiles]

   -l (list format) output all data fields for each record as processed
   -t (tab format) tabulate selected fields for each record as processed
   -c (csv format) comma delimited output of selected fields for each record as processed
   -s (summary format) output summary statistics for each unique value of ligand ID
   -sup (supplier format) tabulate supplier details (from Catalyst)
   -id <IDField> data field to use as ligand ID
   -nh don't output column headings in -t and -c formats
   -o use old (v3.00) score field names as default columns in -t and -c formats, else use v4.00 field names
   -norm use normalised score filed names as default columns in -t and -c formats (normalised = score / #ligand heavy atoms)

.. note::

   If ``-l``, ``-t`` or ``-c`` are combined with ``-s``, the listing/table is
   output withing each ligand summary. ``-sup`` should not be combined with
   other options. Default field names for ``-t`` and ``-c`` are RiboDock score
   field names. Default ID field name is ``Name``. If ``sdFiles`` not given,
   reads from standard input. Output is to standard output.

sdsplit
^^^^^^^

``sdsplit`` -- Splits SD records into multiple files of equal size.
::

   sdsplit [-<RecSize>] [-o <OutputRoot>] [sdFiles]

   -<RecSize> record size to split into (default = 1000 records)
   -o <OutputRoot> Root name for output files (default = tmp)

.. note::

   If SD file list not given, reads from standard input.

sdsort
^^^^^^

``sdsort`` -- Sorts SD records by given data field.
::

   sdsort [-n] [-r] [-f <DataField>] [sdFiles]

   -n numeric sort (default is text sort)
   -r descending sort (default is ascending sort)
   -f <DataField> specifies sort field
   -s fast mode. Sorts the records for each named compound independently (must be consecutive)
   -id <NameField> specifies compound name field (default = 1st title line)

.. note::

   ``_REC`` (record #) is provided as a pseudo-data field. If SD file list not
   given, reads from standard input. Output is to standard output. Fast mode can
   be safely used for partial sorting of huge SD files of raw docking hits
   without running into memory problems.

sdmodify
^^^^^^^^

``sdmodify`` -- Script to set the first title line equal to a given data field.
::

   sdmodify -f <DataField> [sdFiles]

.. note::

   If ``sdFiles`` not given, reads from standard input. Output is to standard
   output.

rbhtfinder
^^^^^^^^^^

``rbhtfinder`` -- Script that simulates the result of a high throughput
protocol.
::

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

An example file would look like as follows:
::

   # 3 steps as the running filters (set by the "3" in next line)
   3
   if - -10 SCORE.INTER 1.0 if - SCORE.NRUNS 9 0.0 -1.0,
   if - -20 SCORE.INTER 1.0 if - SCORE.NRUNS 14 0.0 -1.0,
   if - SCORE.NRUNS 49 0.0 -1.0,
   # 1 writing filter (defined by the "1" in next line)
   1
   - SCORE.INTER -10,

In other (more understandable) words.

First, |Dock| runs 3 consecutive steps:

1. Run 10 runs and check if the ``SCORE.INTER`` is lower than -10, if it is the
   case:
2. Then run 5 more runs (until 15 runs) to see if the ``SCORE.INTER``
   reaches -20. If it is the case:
3. Run up to 50 runs to freely sample the different conformations the molecule
   displays.

And, second:

For the printing information, only print out all those poses where
``SCORE.INTER`` is better than -10 (for avoiding excessive printing).

rblist
^^^^^^

``rblist`` -- Output interaction center info for ligands in SD file (with optional autoionisation).

.. code-block:: bash

   rblist -i <InputSDFile> [-o <OutputSDFile>] [-ap ] [-an] [-allH]

   -i <InputSDFile> - input ligand SD file
   -o <OutputSDFilde> - output SD file with descriptors (default = no output)
   -ap - protonate all neutral amines, guanidines, imidazoles (default = disabled)
   -an - deprotonate all carboxylic, sulphur and phosphorous acid groups (default = disabled)
   -allH - read all hydrogens present (default = polar hydrogens only)
   -tr - rotate all secondary amides to trans (default = leave alone)
   -l - verbose listing of ligand atoms and rotable bonds (default = compact table format)
