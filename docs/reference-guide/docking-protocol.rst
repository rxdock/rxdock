.. _docking-protocol:

Docking protocol
================

Protocol summary
----------------

Pose generation
^^^^^^^^^^^^^^^

|Dock| uses a combination of stochastic and deterministic search techniques to
generate low energy ligand poses. The standard docking protocol to generate a
single ligand pose uses 3 stages of genetic algorithm search (GA1, GA2, GA3),
followed by low temperature Monte Carlo (MC) and simplex minimisation (MIN)
stages.

Several scoring function parameters are varied between the stages to promote
efficient sampling. The ECUT parameter of the :math:`S_{\text{inter}}` vdW
potential (defining the hardness of the intermolecular close range potential)
is increased from 1 in the first GA stage (GA1) to a maximum of 120 in the MC
and MIN stages, with intermediate values of 5 in GA2 and 25 in GA3. The
functional form of the :math:`S_{\text{inter}}` vdW potential is switched from
a 4-8 potential in GA1 and GA2 to a 6-12 potential in GA3, MC and MIN.

In a similar fashion, the overall weight of the :math:`S_{\text{intra}}`
dihedral potential is ramped up from an initial value of 0.1 in GA1 to a final
value of 0.5 in the MC and MIN stages, with intermediate values of 0.2 in GA2
and 0.3 in GA3. In contrast, the :math:`S_{\text{intra}}` vdW parameters (as
used for the ligand intramolecular potential) remain fixed at the final, hard
values throughout the calculation (ECUT = 120, 6-12 potential).

Overall, we found this combination of parameter changes allows for efficient
sampling of the very poor starting poses, whilst minimising the likelihood that
poor ligand internal conformations are artificially favoured and trapped early
in the search, and ensures that physically realistic potentials are used for
final optimisation and analysis.

.. _docking-protocol-genetic-algorithm:

Genetic algorithm
^^^^^^^^^^^^^^^^^

The GA chromosome consists of the ligand centre of mass (com), the ligand
orientation, as represented by the quaternion (q) required to rotate the ligand
principal axes from the Cartesian reference axes, the ligand rotatable dihedral
angles, and the receptor rotatable dihedral angles. The ligand centre of mass
and orientation descriptors, although represented by multiple floating point
values (com.x, com.y, com.z, and q.s, q.x, q.y, q.z respectively), are operated
on as intact entities by the GA mutation and crossover operators.

For so-called free docking, in which no external restraints other than the
cavity penalty are imposed, the initial population is generated such that the
ligand centre of mass is constrained to lie on a randomly selected grid point
within the defined docking volume, and the ligand orientation and all dihedral
angles are randomised completely. Mutations to the ligand centre of mass are by
a random distance along a randomly oriented unit vector. Mutations to the ligand
orientation are performed by rotating the ligand principal axes by a random
angle around a randomly oriented unit vector. Mutations to the ligand and
receptor dihedral angles are by a random angle. All mutation distances and
angles are randomly selected from rectangular distributions of defined width.

A generation is considered to have passed when the number of new individuals
created is equal to the population size. Instead of having a fixed number of
generations, the GA is allowed to continue until the population converges. The
population is considered converged when the score of the best scoring pose fails
to improve by more than 0.1 over the last three generations. This allows early
termination of poorly performing runs for which the initial population is not
able to generate a good solution.

During initial testing the impact of a wide variety of GA parameters
(:numref:`Table %s <table-ga-parameter-space>`) were explored on a small,
representative set of protein-ligand complexes (3ptb, 1rbp, 1stp, 3dfr). We
measured the frequency that the algorithm was able to find the experimental
conformation, and the average run time. Optimum results were obtained with a
steady state GA, roulette wheel selection, a single population of size
(:math:`100 \times \text{(number of rotatable bonds)}`), a crossover:mutation
ratio of 40:60, and mutation distribution widths of ligand translation 2 Å,
ligand rotation of 30 degrees and dihedral angle of 30 degrees. These parameters
have been found to be generally robust across a wide variety of systems.

.. _table-ga-parameter-space:

.. table:: Summary of GA parameter space explored, and final values

   +-----------------------------+--------------------------------------------------+----------------------------------+
   | Parameter                   | Values Explored                                  | Final Values                     |
   +=============================+==================================================+==================================+
   | Number of populations       | 1, 2, 3, 4, 5                                    | 1                                |
   +-----------------------------+--------------------------------------------------+----------------------------------+
   | Selection operator          | Roulette wheel, Rank                             | Roulette wheel                   |
   +-----------------------------+--------------------------------------------------+----------------------------------+
   | Mutation                    | Rectangular Cauchy                               | Rectangular                      |
   +-----------------------------+--------------------------------------------------+----------------------------------+
   | GA                          | Generational, Steady state                       | Steady state                     |
   +-----------------------------+--------------------------------------------------+----------------------------------+
   | Elitism                     | Yes, No                                          | No                               |
   +-----------------------------+--------------------------------------------------+----------------------------------+
   | No of individuals modified  | All values from 1 to population size             | 0.5 \* population size           |
   | in each generation          |                                                  |                                  |
   +-----------------------------+--------------------------------------------------+----------------------------------+
   | Population size             | 50, 75, 100, 125, 150, 200, 400, 800 \* number   | 100 \* number of rotatable bonds |
   |                             | of rotatable bonds                               |                                  |
   +-----------------------------+--------------------------------------------------+----------------------------------+
   | Probability of choosing     | 0.0, 0.05, 0.1 ... 0.9, 0.95, 1.0                | 0.4                              |
   | Crossover vs. Mutation      |                                                  |                                  |
   +-----------------------------+--------------------------------------------------+----------------------------------+
   | Torsion step                | 3, 12, 21, 30 degrees                            | 30 degrees                       |
   +-----------------------------+--------------------------------------------------+----------------------------------+
   | Rotational step             | 3, 12, 21, 30 degrees                            | 30 degrees                       |
   +-----------------------------+--------------------------------------------------+----------------------------------+
   | Translation step            | 0.1, 0.8, 1.4, 2.0 Å                             | 2.0 Å                            |
   +-----------------------------+--------------------------------------------------+----------------------------------+

Monte Carlo
^^^^^^^^^^^

The method and parameters for low temperature Monte Carlo are similar to those
described for phase 4 of the RiboDock simulated annealing search protocol. The
overall number of trials is scaled according to the number of rotatable bonds in
the ligand, from a minimum of 500 (:math:`N_{\text{rot}} = 0`) to a maximum of
2000 (:math:`N_{\text{rot}} = 15`). Maximum step sizes are: translation 0.1 Å,
ligand rotation of 10 degrees and dihedral angle of 10 degrees. Step sizes are
halved if the Metropolis acceptance rate falls below 0.25.

Simplex
^^^^^^^

The Nelder-Mead's simplex minimisation routine operates on the same chromosome
representation as the GA, with the exception that the composite descriptors
(centre of mass and orientation) are decomposed into their constituent floating
point values.

Code implementation
-------------------

Docking protocols are constructed at run-time (by ``RbtTransformFactory`` class)
from docking protocol definition files (|Dock| .prm format). The default location
for docking protocol files is ``$RBT_ROOT/data/scripts/``. The docking protocol
definition file defines the sequence of search algorithms that constitute a
single docking run for a single ligand record. Each search algorithm component
operates either on a single chromosome representing the system degrees of
freedom, or on a population of such chromosomes. The chromosome is constructed
(by ``RbtChromFactory`` class) as an aggregate of individual chromosome elements
for the receptor, ligand and explicit solvent degrees of freedom, as defined by
the flexibility parameters in the system definition file.

.. table:: Chromosome elements

   +-------------+--------------------+-------------------------------+--------+
   | Element     | Defined by         | Class                         | Length |
   +=============+====================+===============================+========+
   | Position    | Centre of mass     | ``RbtChromPositionElement``   | 3      |
   +-------------+--------------------+-------------------------------+--------+
   | Orientation | Euler angles for   | ``RbtChromPositionElement``   | 3      |
   |             | principal axes     |                               |        |
   +-------------+--------------------+-------------------------------+--------+
   | Dihedral    | Dihedral angle for | ``RbtChromDihedralElement``   | 1 per  |
   |             | rotatable bond     |                               | bond   |
   +-------------+--------------------+-------------------------------+--------+
   | Occupancy   | Explicit water     | ``RbtChromOccupancylElement`` | 1 per  |
   |             | occupancy state    |                               | water  |
   +-------------+--------------------+-------------------------------+--------+

Standard |Dock| docking protocol (``dock.prm``)
-----------------------------------------------

As stated above in this section, the standard |Dock| docking protocol consists of
three phases of a genetic algorithm search, followed by low-temperature Monte
Carlo and simplex minimisation.

.. table:: Search algorithm components and C++ implementation classes

   +-------------------------+-------------------------+-----------------------+
   | Component               | Class                   | Operates on           |
   +=========================+=========================+=======================+
   | Randomise population    | ``RbtRandPopTransform`` | Chromosome population |
   +-------------------------+-------------------------+-----------------------+
   | Genetic algorithm       | ``RbtGATransform``      | Chromosome population |
   | search                  |                         |                       |
   +-------------------------+-------------------------+-----------------------+
   | Monte Carlo             | ``RbtSimAnnTransform``  | Single chromosome     |
   | simulated annealing     |                         |                       |
   +-------------------------+-------------------------+-----------------------+
   | Simplex minimisation    | ``RbtSimplexTransform`` | Single chromosome     |
   +-------------------------+-------------------------+-----------------------+
   | Null operation          | ``RbtNullTransform``    | N/A                   |
   +-------------------------+-------------------------+-----------------------+

.. table:: Docking protocol data files

   +------------------------+--------------------------------------------------+
   | File                   | Description                                      |
   +========================+==================================================+
   | ``score.prm``          | Calculates score only for initial conformation   |
   |                        | (standard scoring function)                      |
   +------------------------+--------------------------------------------------+
   | ``scole_solv.prm``     | As above, but uses desolvation scoring function  |
   +------------------------+--------------------------------------------------+
   | ``minimise.prm``       | Simplex minimisation of initial conformation     |
   |                        | (standard scoring function)                      |
   +------------------------+--------------------------------------------------+
   | ``minimise_solv.prm``  | As above, but uses desolvation scoring function  |
   +------------------------+--------------------------------------------------+
   | ``dock.prm``           | Full docking search (standard scoring function)  |
   +------------------------+--------------------------------------------------+
   | ``dock_solv.prm``      | As above, but uses desolvation scoring function  |
   +------------------------+--------------------------------------------------+
   | ``dock_grid.prm``      | Full docking search (standard scoring function,  |
   |                        | grid-based vdW term)                             |
   +------------------------+--------------------------------------------------+
   | ``dock_solv_grid.prm`` | Full docking search (desolvation scoring         |
   |                        | function, grid-based vdW term)                   |
   +------------------------+--------------------------------------------------+

By way of example, the ``dock.prm`` script is documented in detail. The other
scripts are very similar.

.. code-block:: python

   SECTION SCORE
      INTER RbtInterIdxSF.prm
      INTRA RbtIntraSF.prm
      SYSTEM RbtTargetSF.prm
   END_SECTION

**Scoring Function** The scoring function definition is referenced within the
docking protocol definition file itself, in the ``SCORE`` section. This section
contains entries for the ``INTER``, ``INTRA`` and ``SYSTEM`` scoring function
definition files.

.. code-block:: python

   SECTION SETSLOPE_1
      TRANSFORM RbtNullTransform
      # Dock with a high penalty for leaving the cavity
      WEIGHT@SCORE.RESTR.CAVITY 5.0
      # Gradually ramp up dihedral weight from 0.1-->0.5
      WEIGHT@SCORE.INTRA.DIHEDRAL 0.1
      # Gradually ramp up energy cut off for switching to quadratic
      ECUT@SCORE.INTER.VDW 1.0
      # Start docking with a 4-8 vdW potential
      USE 4_8@SCORE.INTER.VDW TRUE
      # Broader angular dependence
      DA1MAX@SCORE.INTER.POLAR 180.0
      # Broader angular dependence
      DA2MAX@SCORE.INTER.POLAR 180.0
      # Broader distance range
      DR12MAX@SCORE.INTER.POLAR 1.5
   END_SECTION

**Genetic algorithm** All sections that contain the ``TRANSFORM`` parameter are
interpreted as a search algorithm component. The value of the ``TRANSFORM``
parameter is the C++ implementation class name for that transform. An
``RbtNullTransform`` can be used to send messages to the scoring function to
modify key scoring function parameters in order to increase search efficiency.
All parameter names that contain the ``@`` symbol are interpreted as scoring
function messages, where the string before the ``@`` is the scoring function
parameter name, the string after the ``@`` is the scoring function term, and the
parameter value is the new value for the scoring function parameter. Messages
are sent blind, with no success feedback, as the docking protocol has no
knowledge of the composition of the scoring function terms.

Here, we start the docking with a soft 4-8 vdW potential, a reduced dihedral
potential, and extended polar ranges (distances and angles) for the
intermolecular polar potential. These changes are all designed to aid sampling
efficiency by not overpenalising bad contacts in the initial, randomised
population, and by encouraging the formation of intermolecular hydrogen bonds.

.. code-block:: python

   SECTION RANDOM_POP
      TRANSFORM RbtRandPopTransform
      POP_SIZE 50
      SCALE_CHROM_LENGTH TRUE
   END_SECTION

Creates an initial, randomised chromosome population. If ``SCALE_CHROM_LENGTH``
is false, the population is of fixed size, defined by ``POP_SIZE``. If
``SCALE_CHROM_LENGTH`` is true, the population is proportional to the overall
chromosome length, defined by ``POP_SIZE`` multiplied by the chromosome length.

.. code-block:: python

   SECTION GA_SLOPE1
      TRANSFORM RbtGATransform
      PCROSSOVER 0.4 # Prob. of crossover
      XOVERMUT TRUE # Cauchy mutation after each crossover
      CMUTATE FALSE # True = Cauchy; False = Rectang. for regular mutations
      STEP_SIZE 1.0 # Max relative mutation
   END_SECTION

First round of GA.

.. code-block:: python

   SECTION SETSLOPE_3
      TRANSFORM RbtNullTransform
      WEIGHT@SCORE.INTRA.DIHEDRAL 0.2
      ECUT@SCORE.INTER.VDW 5.0
      DA1MAX@SCORE.INTER.POLAR 140.0
      DA2MAX@SCORE.INTER.POLAR 140.0
      DR12MAX@SCORE.INTER.POLAR 1.2
   END_SECTION

Increases the ligand dihedral weight, increases the short-range intermolecular
vdW hardness (ECUT), and decreases the range of the intermolecular polar
distances and angles.

.. code-block:: python

   SECTION GA_SLOPE3
      TRANSFORM RbtGATransform
      PCROSSOVER 0.4 # Prob. of crossover
      XOVERMUT TRUE # Cauchy mutation after each crossover
      CMUTATE FALSE # True = Cauchy ; False = Rectang. for regular mutations
      STEP_SIZE 1.0 # Max relative mutation
   END_SECTION

Second round of GA with revised scoring function parameters.

.. code-block:: python

   SECTION SETSLOPE_5
      TRANSFORM RbtNullTransform
      WEIGHT@SCORE.INTRA.DIHEDRAL 0.3
      ECUT@SCORE.INTER.VDW 25.0
      # Now switch to a convential 6-12 for final GA, MC, minimisation
      USE 4_8@SCORE.INTER.VDW FALSE
      DA1MAX@SCORE.INTER.POLAR 120.0
      DA2MAX@SCORE.INTER.POLAR 120.0
      DR12MAX@SCORE.INTER.POLAR 0.9
   END_SECTION

Further increases the ligand dihedral weight, further increases the short-range
intermolecular vdW hardness (ECUT), and further decreases the range of the
intermolecular polar distances and angles. Also switches from softer 4-8 vdW
potential to a harder 6-12 potential for final round of GA, MC and minimisation.

.. code-block:: python

   SECTION GA_SLOPE5
      TRANSFORM RbtGATransform
      PCROSSOVER 0.4 # Prob. of crossover
      XOVERMUT TRUE # Cauchy mutation after each crossover
      CMUTATE FALSE # True = Cauchy ; False = Rectang. for regular mutations
      STEP_SIZE 1.0 # Max relative mutation
   END_SECTION

Final round of GA with revised scoring function parameters.

.. code-block:: python

   SECTION SETSLOPE_10
      TRANSFORM RbtNullTransform
      WEIGHT@SCORE.INTRA.DIHEDRAL 0.5 # Final dihedral weight matches SF file
      ECUT@SCORE.INTER.VDW 120.0 # Final ECUT matches SF file
      DA1MAX@SCORE.INTER.POLAR 80.0
      DA2MAX@SCORE.INTER.POLAR 100.0
      DR12MAX@SCORE.INTER.POLAR 0.6
   END_SECTION

Resets all the modified scoring function parameters to their final values,
corresponding to the values in the scoring function definition files. It is
important that the final scoring function optimised by the docking search can be
compared directly with the score-only and minimisation-only protocols, in which
the scoring function parameters are not modified.

.. code-block:: python

   SECTION MC_10K
      TRANSFORM RbtSimAnnTransform
      START_T 10.0
      FINAL_T 10.0
      NUM_BLOCKS 5
      STEP_SIZE 0.1
      MIN_ACC_RATE 0.25
      PARTITION_DIST 8.0
      PARTITION_FREQ 50
      HISTORY_FREQ 0
   END_SECTION

**Monte Carlo** Low temperature Monte Carlo sampling, starting from fittest
chromosome from final round of GA.

.. code-block:: python

   SECTION SIMPLEX
      TRANSFORM RbtSimplexTransform
      MAX_CALLS 200
      NCYCLES 20
      STOPPING_STEP_LENGTH 10e-4
      PARTITION_DIST 8.0
      STEP_SIZE 1.0
      CONVERGENCE 0.001
   END_SECTION

**Minimisation** Simplex minimisation, starting from fittest chromosome from low
temperature Monte Carlo sampling.

.. code-block:: python

   SECTION FINAL
      TRANSFORM RbtNullTransform
      WEIGHT@SCORE.RESTR.CAVITY 1.0 # revert to standard cavity penalty
   END_SECTION

Finally, we reset the cavity restraint penalty to 1. The penalty has been held
at a value of 5 throughout the search, to strongly discourage the ligand from
leaving the docking site.
