Scoring functions
=================

Component scoring functions
---------------------------

The |Dock| master scoring function (:math:`S_{\text{total}}`) is a weighted sum
of intermolecular (:math:`S_{\text{inter}}`), ligand intramolecular
(:math:`S_{\text{intra}}`), site intramolecular (:math:`S_{\text{site}}`), and
external restraint terms (:math:`S_{\text{restraint}}`) :eq:`(%s) <eqnStotal>`.
:math:`S_{\text{inter}}` is the main term of interest as it represents the
protein-ligand (or RNA-ligand) interaction score :eq:`(%s) <eqnSinter>`.
:math:`S_{\text{intra}}` represents the relative energy of the ligand
conformation :eq:`(%s) <eqnSintra>`. Similarly, :math:`S_{\text{site}}`
represents the relative energy of the flexible regions of the active site
:eq:`(%s) <eqnSsite>`. In the current implementation, the only flexible bonds in
the active site are terminal OH and NH3+ bonds. :math:`S_{\text{restraint}}` is
a collection of non-physical restraint functions that can be used to bias the
docking calculation in several useful ways :eq:`(%s) <eqnSrestraint>`.

.. math::
   :label: eqnStotal

   S^{\text{total}} = S^{\text{inter}} + S^{\text{intra}} + S^{\text{site}} + S^{\text{restraint}}

.. math::
   :label: eqnSinter

   S^{\text{inter}} = W_{\text{vdW}}^{\text{inter}} \cdot S_{\text{vdW}}^{\text{inter}}
                    + W_{\text{polar}}^{\text{inter}} \cdot S_{\text{polar}}^{\text{inter}}
		    + W_{\text{repul}}^{\text{inter}} \cdot S_{\text{repul}}^{\text{inter}}
		    + W_{\text{arom}}^{\text{inter}} \cdot S_{\text{arom}}^{\text{inter}} +\\
   + W_{\text{solv}} \cdot S_{\text{solv}} + W_{\text{rot}} \cdot N_{\text{rot}} + W_{\text{const}}

.. math::
   :label: eqnSintra

   S^{\text{intra}} = W_{\text{vdW}}^{\text{intra}} \cdot S_{\text{vdW}}^{\text{intra}}
                    + W_{\text{polar}}^{\text{intra}} \cdot S_{\text{polar}}^{\text{intra}}
		    + W_{\text{repul}}^{\text{intra}} \cdot S_{\text{repul}}^{\text{intra}}
		    + W_{\text{dihedral}}^{\text{intra}} \cdot S_{\text{dihedral}}^{\text{intra}}

.. math::
   :label: eqnSsite

   S^{\text{site}} = W_{\text{vdW}}^{\text{site}} \cdot S_{\text{vdW}}^{\text{site}}
                   + W_{\text{polar}}^{\text{site}} \cdot S_{\text{polar}}^{\text{site}}
		   + W_{\text{repul}}^{\text{site}} \cdot S_{\text{repul}}^{\text{site}}
		   + W_{\text{dihedral}}^{\text{site}} \cdot S_{\text{dihedral}}^{\text{site}}

.. math::
   :label: eqnSrestraint

   S^{\text{restraint}} = W_{\text{cavity}} \cdot S_{\text{cavity}} + W_{\text{tether}} \cdot S_{\text{tether}}
                        + W_{\text{nmr}} \cdot S_{\text{nmr}} + W_{\text{ph4}} \cdot S_{\text{ph4}}

:math:`S_{\text{inter}}`, :math:`S_{\text{intra}}`, and :math:`S_{\text{site}}`
are built from a common set of constituent potentials, which are described
below. The main changes to the original RiboDock scoring function
[RiboDock2004]_ are:

i.   the replacement of the crude steric potentials (:math:`S_{\text{lipo}}` and
     :math:`S_{\text{rep}}`) with a true van der Waals potential,
     :math:`S_{\text{vdW}}`
ii.  the introduction of two generalised terms for all short range attractive
     (:math:`S_{\text{polar}}`) and repulsive (:math:`S_{\text{repul}}`) polar
     interactions
iii. the implementation of a fast weighed solvent accessible surface (WSAS) area
     solvation term
iv.  the addition of explicit dihedral potentials

van der Waals potential
^^^^^^^^^^^^^^^^^^^^^^^

We have replaced the :math:`S_{\text{lipo}}` and :math:`S_{\text{rep}}`
empirical potentials used in RiboDock with a true vdW potential similar to that
used by GOLD [GOLD2005]_. Atom types and vdW radii were taken from the
Tripos 5.2 force field and are listed in the :ref:`Appendix section
<user-guide-appendix>` (:numref:`Table %s <table-vdw-parameters-tripos-52>`).
Energy well depths are switchable between the original Tripos 5.2 values and
those used by GOLD, which are calculated from the atomic polarisability and
ionisation potentials of the atom types involved. Additional atom types were
created for carbons with implicit hydrogens, as the Tripos force field uses an
all-atom representation. vdW radii for implicit hydrogen types are increased by
0.1 Å for each implicit hydrogen, with well depths unchanged. The functional
form is switchable between a softer 4-8 and a harder 6-12 potential. A quadratic
potential is used at close range to prevent excessive energy penalties for
atomic clashes. The potential is truncated at longer range
(:math:`1.5 \cdot r_{\min}`, the sum of the vdW radii).

The quadratic potential is used at repulsive energies between
:math:`e_{\text{cutoff}}` and :math:`e_0`, where :math:`e_{\text{cutoff}}` is
defined as a multiple of the energy well depth
(:math:`e_{\text{cutoff}} = \text{ECUT} \cdot e_{\min}`), and :math:`e_0` is the
energy at zero separation, defined as a multiple of :math:`e_{\text{cutoff}}`
(:math:`e_0 = \text{E0} \cdot e_{\text{cutoff}}`). ECUT can vary between 1 and
120 during the docking search (see :ref:`Genetic algorithm subsection
<docking-protocol-genetic-algorithm>` in :ref:`Docking protocol section
<docking-protocol>`), whereas E0 is fixed at 1.5.

Empirical attractive and repulsive polar potentials
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We continue to use an empirical Bohm-like potential to score hydrogen-bonding
and other short-range polar interactions. The original RiboDock polar terms
(:math:`S_{\text{H-bond}}`, :math:`S_{\text{posC-acc}}`,
:math:`S_{\text{acc-acc}}`, :math:`S_{\text{don-don}}`) are generalised and
condensed into two scoring functions, :math:`S_{\text{polar}}` and
:math:`S_{\text{repul}}` (:eq:`(%s) <eqnSpolar>` and :eq:`(%s) <eqnSrepul>`,
also taking into account :eq:`(%s) <eqnf1DeltaX>`, :eq:`(%s) <eqnf2i>`,
:eq:`(%s) <eqnsgni>`, :eq:`(%s) <eqnci>`, :eq:`(%s) <eqnf3DeltaX>`,
:eq:`(%s) <eqnNi>`, :eq:`(%s) <eqnSsolv>`, :eq:`(%s) <eqnSsolvrs>`,
:eq:`(%s) <eqnSsolvpij>`, :eq:`(%s) <eqnSsolvSi>`, :eq:`(%s) <eqnSsolvbij>`,
:eq:`(%s) <eqnSsolvAi>`, and :eq:`(%s) <eqnSsolvDeltaGWSAS>`), which deal with
attractive and repulsive interactions respectively. Six types of polar
interaction centres are considered: hydrogen bond donors (DON), metal ions (M+),
positively charged carbons (C+, as found at the centre of guanidinium,
amidinium and imidazole groups), hydrogen bond acceptors with pronounced lone
pair directionality (ACC_LP), acceptors with in-plane preference but limited
lone-pair directionality (ACC_PLANE), and all remaining acceptors (ACC). The
ACC_LP type is used for carboxylate oxygens and O sp2 atoms in RNA bases, with
ACC_PLANE used for other O sp2 acceptors. This distinction between acceptor
types was not made in RiboDock, in which all acceptors were implicitly of type
ACC.

.. math::
   :label: eqnSpolar

   S_{\text{polar}} = \sum_{\text{IC1-IC2}} f_1(|\Delta R_{12}|) \cdot \text{ANG}_{\text{IC1}}
                    \cdot \text{ANG}_{\text{IC2}} \cdot f_2(\text{IC1}) \cdot f_2(\text{IC2})
		    \cdot f_3(\text{IC1}) \cdot f_3(\text{IC2})

.. math::
   :label: eqnSrepul

   S_{\text{repul}} = \sum_{\text{IC1-IC2}} f_1(\Delta R_{12}) \cdot \text{ANG}_{\text{IC1}}
                    \cdot \text{ANG}_{\text{IC2}} \cdot f_2(\text{IC1}) \cdot f_2(\text{IC2})
		    \cdot f_3(\text{IC1}) \cdot f_3(\text{IC2})

.. math::
   :label: eqnf1DeltaX

   f_1(\Delta X) = \begin{cases}
                      1                                                                        & \Delta X \leq \Delta X_{\min}\\
                      1 - \frac{\Delta X - \Delta X_{\min}}{\Delta X_{\max} - \Delta X_{\min}} & \Delta X_{\min} < \Delta X \leq \Delta X_{\max}\\
                      0                                                                        & \Delta X > \Delta X_{\max}
                   \end{cases}

.. math::
   :label: eqnf2i

   f_2(i) = sgn(i)(1 + 0.5 |c_i|)

.. math::
   :label: eqnsgni

   sgn(i) = \begin{cases}
               -1   & \text{ACC}, \text{ACC}\_\text{LP}, \text{ACC}\_\text{PLANE}\\
               +0.5 & \text{C+}\\
               +1.0 & \text{DON}, \text{M+}
            \end{cases}

.. math::
   :label: eqnci

   c_i = \text{formal charge on primary atom of interaction centre \(i\)}

.. math::
   :label: eqnf3DeltaX

   f_3(\Delta X) = \begin{cases}
                      \sqrt{\frac{N_i}{25}} & \text{macromolecular interaction centres}\\
                      1                     & \text{ligand interaction centres}
                   \end{cases}

.. math::
   :label: eqnNi

   N_i = \text{number of non-hydrogen macromolecule atoms within}\\
         \text{5 Å radius of primary atom of interaction centre \(i\)}

Individual interaction scores are the product of simple scaling functions for
geometric variables, formal charges and local neighbour density. The scaling
functions themselves, and the formal charge assignment method, are retained from
RiboDock. Metals are assigned a uniform formal charge of +1. C+ is considered to
be a weak donor in this context and scores are scaled by 50 % relative to
conventional donors by the assignment of :math:`sgn(i) = 0.5` in
:eq:`(%s) <eqnf2i>`. Pseudo-formal charges are no longer assigned to selected
RNA base atoms. The geometric functions minimally include an interaction
distance term, with the majority also including angular terms dependent on the
type of the interaction centres. Geometric parameters and the angular functions
are summarised in :ref:`Appendix section <user-guide-appendix>` (:numref:`Table
%s <table-geometric-parameters-empirical-terms>` and :numref:`Table %s
<table-angular-functions-polar-interactions>`, respectively).

The most notable improvements to RiboDock are that attractive (hydrogen bond and
metal) interactions with ACC_LP and ACC_PLANE acceptors include terms for
:math:`\phi` and :math:`\theta` (as defined in [ref 3]) to enforce the
relevant lone pair directionality. These replace the :math:`\alpha_{\text{ACC}}`
dependence, which is retained for the ACC acceptor type. No distinction between
acceptor types is made for attractive interactions with C+ carbons, or for
repulsive interactions between acceptors. In these circumstances all acceptors
are treated as type ACC. Such C+-ACC interactions, which in RiboDock were
described by only a distance function, (:math:`S_{\text{posC-acc}}`) now include
angular functions around the carbon and acceptor groups. Repulsive interactions
between donors, and between acceptors, also have an angular dependence. This
allows a stronger weight, and a longer distance range, to be used to penalise
disallowed head-to-head interactions without forbidding allowable contacts. One
of the issues in RiboDock was that it was not possible to include neutral
acceptors in the acceptor-acceptor repulsion term with a simple distance
function.

Solvation potential
^^^^^^^^^^^^^^^^^^^

The desolvation potential in |Dock| combines a weighted solvent-accessible
surface area approach [WSAS2001]_ with a rapid probabilistic approximation to
the calculation of solvent-accessible surface areas [RASASA1988]_ based on
pairwise interatomic distances and radii (:eq:`(%s) <eqnSsolv>`, taking into
account :eq:`(%s) <eqnSsolvrs>`, :eq:`(%s) <eqnSsolvpij>`,
:eq:`(%s) <eqnSsolvSi>`, :eq:`(%s) <eqnSsolvbij>`, :eq:`(%s) <eqnSsolvAi>`,
and :eq:`(%s) <eqnSsolvDeltaGWSAS>`).

.. math::
   :label: eqnSsolv

   S_{\text{solv}} = (\Delta G_{\text{WSAS}}^{\text{site,bound}}
                   - \Delta G_{\text{WSAS}}^{\text{site\(_0\),unbound}})
                   + (\Delta G_{\text{WSAS}}^{\text{ligand,bound}}
		   - \Delta G_{\text{WSAS}}^{\text{ligand\(_0\),unbound}})

.. math::
   :label: eqnSsolvrs

   r_{\text{s}} = 0.6 \text{Å}

.. math::
   :label: eqnSsolvpij

   p_{ij} = \begin{cases}
               0.8875 & \quad \text{1-2 intramolecular connections}\\
               0.3516 & \quad \text{1-3 intramolecular connections}\\
               0.3156 & \quad \text{1-4 intramolecular connections and above}\\
               0.3156 & \quad \text{intermolecular interactions}
            \end{cases}

.. math::
   :label: eqnSsolvSi

   S_i = 4 \pi (r_i + r_s)^2

.. math::
   :label: eqnSsolvbij

   b_{ij} = \pi (r_i + r_s) (r_j + r_i + 2 r_s - d) \Big(1 - \frac{r_j - r_i}{d}\Big)

.. math::
   :label: eqnSsolvAi

   A_i = S_i \prod_j 1 - \frac{p_i p_{ij} b_{ij}}{S_i}

.. math::
   :label: eqnSsolvDeltaGWSAS

   \Delta G_{\text{WSAS}} = \sum_{i=1}^{n_i} w_i A_i

The calculation is fast enough therefore to be used in docking. We have
redefined the solvation atom types compared to the original work [WSAS2001]_ and
recalibrated the weights against the same training set of experimental solvation
free energies in water (395 molecules). The total number of atom types (50,
including 6 specifically for ionic groups and metals) is slightly lower than in
original work (54). Our atom types reflect the fact that |Dock| uses implicit
non-polar hydrogens. The majority of types are a combination of hybridisation
state and the number of implicit or explicit hydrogens. All solvation parameters
are listed in :ref:`Appendix section <user-guide-appendix>` (:numref:`Table %s
<table-solvation-parameters>`).

:math:`S_{\text{solv}}` is calculated as the change in solvation energy of the
ligand and the docking site upon binding of the ligand. The reference energies
are taken from the initial conformations of the ligand and site (as read from
file) and not from the current pose under evaluation. This is done to take into
account any changes to intramolecular solvation energy. Strictly speaking the
intramolecular components should be reported separately under
:math:`S_{\text{intra}}` and :math:`S_{\text{site}}` but this is not done for
reasons of computational efficiency.

Dihedral potential
^^^^^^^^^^^^^^^^^^

Dihedral energies are calculated using Tripos 5.2 dihedral parameters for all
ligand and site rotatable bonds. Corrections are made to account for the missing
contributions from the implicit non-polar hydrogens.

Intermolecular scoring functions under evaluation
-------------------------------------------------

Training sets
^^^^^^^^^^^^^

We constructed a combined set of protein-ligand and RNA-ligand complexes for
training of |Dock|. Molecular data files for the protein-ligand complexes were
extracted from the downloaded CCDC/Astex cleanlist [ASTEX2007]_ and used without
substantive modification. The only change was to convert ligand MOL2 files to
MDL SD format using Corina [CORINA1990]_, leaving the coordinates and
protonation states intact.

Protein MOL2 files were read directly. The ten RNA-ligand NMR structures from
the RiboDock validation set were supplemented with five RNA-ligand crystal
structures (1f1t, 1f27, 1j7t, 1lc4, 1mwl) prepared in a similar way. All 15
RNA-ligand structures have measured binding affinities.

58 complexes (43 protein-ligand and 15 RNA-ligand) were selected for the initial
fitting of component scoring function weights. Protein-ligand structures were
chosen (of any X-ray resolution) that had readily available experimental binding
affinities [PDBbind2004]_. 102 complexes were used for the main validation of
native docking accuracy for different scoring function designs, consisting of 87
of the 92 entries in the high resolution (R < 2 Å) clean-list (covalently bound
ligands removed -- 1aec, 1b59, 1tpp, 1vgc, 4est), and the 15 RNA-ligand
complexes.

Scoring functions design
^^^^^^^^^^^^^^^^^^^^^^^^

Component weights (:math:`W`) for each term in the intermolecular scoring
function (:math:`S_{\text{inter}}`) were obtained by least squares regression
of the component scores to :math:`\Delta G_{\text{bind}}` values for the binding
affinity training set described above (58 entries). Each ligand was subjected
first to simplex minimisation in the docking site, starting from the
crystallographic pose, to relieve any minor non-bonded clashes with the site.
The scoring function used for minimisation was initialised with reasonable
manually assigned weights. If the fitted weights deviated significantly from the
initial weights the procedure was repeated until convergence. Certain weights
(:math:`W_{\text{repul}}`, :math:`W_{\text{rot}}`, :math:`W_{\text{const}}`)
were constrained to have positive values to avoid non-physical, artefactual
models. Note that the presence of :math:`W_{\text{rot}}` and
:math:`W_{\text{const}}` in :math:`S_{\text{inter}}` improves the quality of the
fit to the binding affinities but does not impact on native ligand docking
accuracy.

Ten intermolecular scoring functions were derived with various combinations of
terms (:numref:`Table %s <table-intermolecular-scoring-function-weights>`). SF0
is a baseline scoring function that has the van der Waals potential only. SF1
adds a simplified polar potential, without f2 (formal charge) and f3 (neighbour
density) scaling functions, and with a single acceptor type (ACC) that lacks
lone-pair directionality. SF2 has the full polar potential (f2 and f3 scaling
functions, ACC, ACC_LP and ACC_PLANE acceptor types) and adds the repulsive
polar potential. SF3 has the same functional form as SF2 but with empirical
weights in regular use at RiboTargets. SF4 replaces the repulsive polar
potential with the WSAS desolvation potential described above. SF5 has the same
functional form as SF4 but with empirical weights in regular use at RiboTargets.
SF6 combines the repulsive polar and desolvation potentials. SF7 has the same
functional form as SF2 and SF3 but with weights for :math:`W_{\text{vdW}}` and
:math:`W_{\text{polar}}` taken from SF5. SF8 and SF9 add the crude aromatic term
from RiboDock [RiboDock2004]_ to SF3 and SF5 respectively. The
:math:`S_{\text{intra}}` functional form and weights were held constant, and
equivalent to SF3, to avoid any differences in ligand conformational energies
affecting the docking results. As the :math:`S_{\text{site}}` scores are
calculated simultaneously with :math:`S_{\text{inter}}` (for computational
reasons) the :math:`S_{\text{site}}` functional form and weights vary in line
with :math:`S_{\text{inter}}`. There is surprisingly little variation in
correlation coefficient (R) and root mean square error (RMSE) in predicted
binding energy over the ten scoring functions (:numref:`Table %s
<table-intermolecular-scoring-function-weights>`). The best results are obtained
with SF4 (R = 0.67, RMSE = 9.6 kJ/mol).

.. _table-intermolecular-scoring-function-weights:

.. table:: Intermolecular scoring function weights under evaluation (a =
	   constrained to be > zero; b = fixed values; c = correlation
	   coefficient (R), and root mean squared error (RMSE) between
	   :math:`S_{\text{inter}}` and :math:`\Delta G_{\text{bind}}`, for
	   minimised experimental ligand poses, over binding affinity validation
	   set (58 entries)).

   +----+------------------------+--------------------------+-------------------------+---------------------------+-------------------------+-------------------------+---------------------------+----------+----------+
   | SF | :math:`W_{\text{vdW}}` | :math:`W_{\text{polar}}` | :math:`W_{\text{solv}}` | :math:`W_{\text{repul}}`\ | :math:`W_{\text{arom}}` | :math:`W_{\text{rot}}`\ | :math:`W_{\text{const}}`\ | R\       | RMSE\    |
   |    |                        |                          |                         | :sup:`a`                  |                         | :sup:`a`                | :sup:`a`                  | :sup:`c` | :sup:`c` |
   +====+========================+==========================+=========================+===========================+=========================+=========================+===========================+==========+==========+
   | 0  | 1.4                    | -\                       | -\                      | -\                        | -\                      | 0                       | 0                         | 0.62     | 10.9     |
   +----+------------------------+--------------------------+-------------------------+---------------------------+-------------------------+-------------------------+---------------------------+----------+----------+
   | 1  | 1.126                  | 2.36                     | -\                      | -\                        | -\                      | 0.217                   | 0                         | 0.64     | 10.2     |
   +----+------------------------+--------------------------+-------------------------+---------------------------+-------------------------+-------------------------+---------------------------+----------+----------+
   | 2  | 1.192                  | 2.087                    | -\                      | 2.984                     | -\                      | 0                       | 0                         | 0.63     | 10.4     |
   +----+------------------------+--------------------------+-------------------------+---------------------------+-------------------------+-------------------------+---------------------------+----------+----------+
   | 3  | 1.000\ :sup:`b`        | 3.400\ :sup:`b`          | -\                      | 5.000\ :sup:`b`           | -\                      | 0                       | 0                         | 0.59     | 10.9     |
   +----+------------------------+--------------------------+-------------------------+---------------------------+-------------------------+-------------------------+---------------------------+----------+----------+
   | 4  | 1.317                  | 3.56                     | 0.449                   | -\                        | -\                      | 0                       | 4.\                       | 0.67     | 9.6      |
   +----+------------------------+--------------------------+-------------------------+---------------------------+-------------------------+-------------------------+---------------------------+----------+----------+
   | 5  | 1.500\ :sup:`b`        | 5.000\ :sup:`b`          | 0.500\ :sup:`b`         | -\                        | -\                      | 0.568                   | 4.782                     | 0.62     | 10.7     |
   +----+------------------------+--------------------------+-------------------------+---------------------------+-------------------------+-------------------------+---------------------------+----------+----------+
   | 6  | 1.314                  | 4.447                    | 0.500\ :sup:`b`         | 5.000\ :sup:`b`           | -\                      | 0                       | 0                         | 0.62     | 10.4     |
   +----+------------------------+--------------------------+-------------------------+---------------------------+-------------------------+-------------------------+---------------------------+----------+----------+
   | 7  | 1.500\ :sup:`b`        | 5.000\ :sup:`b`          | -\                      | 5.000\ :sup:`b`           | -\                      | 0.986                   | 12.046                    | 0.55     | 12.9     |
   +----+------------------------+--------------------------+-------------------------+---------------------------+-------------------------+-------------------------+---------------------------+----------+----------+
   | 8  | 1.000\ :sup:`b`        | 3.400\ :sup:`b`          | -\                      | 5.000\ :sup:`b`           | -1.6\ :sup:`b`          | 0                       | 0                         | 0.53     | 11.8     |
   +----+------------------------+--------------------------+-------------------------+---------------------------+-------------------------+-------------------------+---------------------------+----------+----------+
   | 9  | 1.500\ :sup:`b`        | 5.000\ :sup:`b`          | 0.500\ :sup:`b`         | -\                        | -1.6\ :sup:`b`          | 0.647                   | 5.056                     | 0.58     | 11.5     |
   +----+------------------------+--------------------------+-------------------------+---------------------------+-------------------------+-------------------------+---------------------------+----------+----------+

Scoring functions validation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ability of the ten intermolecular scoring functions (SF0 to SF9) to
reproduce known ligand binding modes was determined on the combined test set of
102 protein-ligand and RNA-ligand complexes. The intra-ligand scoring function
(:math:`S_{\text{intra}}`) was kept constant, with component weights equivalent
to SF3, and a dihedral weight of 0.5. Terminal OH and NH3 groups on the receptor
in the vicinity of the docking site were fully flexible during docking. Ligand
pose populations of size :math:`N_{\text{pop}} = 300` were collected for each
complex and intermolecular scoring function combination. The population size was
increased to :math:`N_{\text{pop}} = 1000` for two of the most promising scoring
functions (SF3 and SF5).

Protein-ligand docking accuracy is remarkably insensitive to scoring function
changes. Almost half of the ligand binding modes can be reproduced with a vdW
potential only (SF0). The addition of a simplified polar potential (SF1)
increases the accuracy to over 70 % of protein-ligand test cases predicted to
within 2 Å RMSD. The success rate increases further to 78 % with SF3, which has
the full attractive and repulsive polar potentials, and empirically adjusted
weights relative to SF2. Subsequent changes to the component terms and weights,
including the addition of the desolvation potential, have little or no impact on
the protein-ligand RMSD metric.

The nucleic acid set shows a much greater sensitivity to scoring function
changes. This can in part be explained by the smaller sample size that amplifies
the percentage changes in the RMSD metric, but nevertheless the trends are
clear. There is a gradual increase in docking accuracy from SF0 (37 %) to SF3
(52 %), but absolute performance is much lower than for the protein-ligand test
set. This level of docking accuracy for nucleic acid-ligand complexes is broadly
consistent with the original RiboDock scoring function, despite the fact that
the original steric term (LIPO) has been replaced by a true vdW potential. The
introduction of the desolvation potential in place of the empirical repulsive
polar potential (in SF4 and SF5) results in a substantial improvement in
accuracy, to around 70 % of test cases within 2 Å RMSD. Subsequent changes (SF6
to SF9) degrade the accuracy. The lower performance of SF7, which has the higher
weights for the VDW and POLAR terms taken from SF5 but lacks the desolvation
potential, demonstrates that it is the desolvation term itself that is having
the beneficial effect, and not merely the reweighting of the other terms. The
inclusion of the geometric aromatic term in SF8 and SF9 has a detrimental
impact on the performance of SF3 and SF5 respectively.

Overall, SF5 achieves optimum performance across proteins and nucleic acids
(76.7 % within 2 Å RMSD). SF3 (no desolvation potential) and SF5 (with
desolvation potential) were selected as the best intermolecular scoring
functions. Finally, these two scoring functions, SF3 and SF5, were the ones
implemented in |Dock| with the names of ``dock.prm`` and ``dock_solv.prm``,
respectively.

.. note::

   In virtual screening campaigns, or in experiments where score of different
   ligands is compared, the best scoring poses for each molecule (as defined by
   the lowest :math:`S_{\text{total}}` within the sample) are sorted and ranked
   by :math:`S_{\text{inter}}`. In other words, the contributions to
   :math:`S_{\text{total}}` from :math:`S_{\text{intra}}`,
   :math:`S_{\text{site}}` and :math:`S_{\text{restraint}}` are ignored when
   comparing poses between different ligands against the same target. The
   rationale for this is that, in particular, the ligand intramolecular scores
   are not on an absolute scale and can differ markedly between different
   ligands.

Code implementation
-------------------

Scoring functions for docking are constructed at run-time (by ``RbtSFFactory``
class) from scoring function definition files (|Dock| ``.prm`` format). The
default location for scoring function definition files is
``$RBT_ROOT/data/sf/``.

The total score is an aggregate of intermolecular ligand-receptor and
ligand-solvent interactions (branch ``SCORE.INTER``), intra-ligand interactions
(branch ``SCORE.INTRA``), intra-receptor, intra-solvent and receptor-solvent
interactions (branch ``SCORE.SYSTEM``), and external restraint penalties (branch
``SCORE.RESTR``).

The ``SCORE.INTER``, ``SCORE.INTRA`` and ``SCORE.SYSTEM`` branches consist of
weighted sums of interaction terms as shown below. Note that not all terms
appear in all branches. See the rDock draft paper [rDock2014]_ for more details
on the implementation of these terms.

.. table:: Scoring function terms and C++ implementation classes.

   +-------+---------------------------------------------+-------------------+---------------------+-------------------+
   | Term  | Description                                 | ``INTER``         | ``INTRA``           | ``SYSTEM``        |
   +=======+=============================================+===================+=====================+===================+
   | VDW   | van der Waals                               | ``RbtVdWIdxSF``   | ``RbtVdwIntraSF``   | ``RbtVdwIdxSF``   |
   +-------+---------------------------------------------+-------------------+---------------------+-------------------+
   | VDW   | van der Waals (grid based)                  | ``RbtVdwGridSF``  | N/A                 | N/A               |
   +-------+---------------------------------------------+-------------------+---------------------+-------------------+
   | POLAR | Attractive polar                            | ``RbtPolarIdxSF`` | ``RbtPolarIntraSF`` | ``RbtPolarIdxSF`` |
   +-------+---------------------------------------------+-------------------+---------------------+-------------------+
   | REPUL | Repulsive polar                             | ``RbtPolarIdxSF`` | ``RbtPolarIntraSF`` | ``RbtPolarIdxSF`` |
   +-------+---------------------------------------------+-------------------+---------------------+-------------------+
   | SOLV  | Desolvation                                 | ``RbtSAIdxSF``    | ``RbtSAIdxSF``      | ``RbtSAIdxSF``    |
   +-------+---------------------------------------------+-------------------+---------------------+-------------------+
   | CONST | Translation/rotational binding entropy      | ``RbtConstSF``    | N/A                 | ``RbtConstSF``    |
   |       | penalty                                     |                   |                     |                   |
   +-------+---------------------------------------------+-------------------+---------------------+-------------------+
   | ROT   | Torsional binding entropy penalty           | ``RbtRotSF``      | N/A                 | N/A               |
   +-------+---------------------------------------------+-------------------+---------------------+-------------------+

Two intermolecular scoring functions (``SCORE.INTER`` branch) have been
validated. These are known informally as the standard scoring function and the
desolvation scoring function (referred to as SF3 and SF5 respectively in the
rDock draft paper [rDock2014]_). The standard intermolecular scoring function
consists of VDW, POLAR and REPUL terms. In the desolvation scoring function, the
REPUL term is replaced by a more finely parameterised desolvation potential
(SOLV term) based on a weighted solvent-accessible surface (WSAS) area model.
The ligand intramolecular scoring function (``SCORE.INTRA`` branch) remains
constant in both cases, and has the same terms and weights as the standard
intermolecular scoring function.

.. table:: Scoring function data files.

   +----------------------------+--------------------------------------------------------------------------------------+
   | File                       | Description                                                                          |
   +============================+======================================================================================+
   | ``RbtInterIdxSF.prm``      | Intermolecular scoring function definition (standard scoring function, SF3)          |
   +----------------------------+--------------------------------------------------------------------------------------+
   | ``RbtInterGridSF.prm``     | As above, but vdW term uses a precalculated grid                                     |
   +----------------------------+--------------------------------------------------------------------------------------+
   | ``RbtSolvIdxSF.prm``       | Intermolecular scoring function definition (desolvation scoring function, SF5)       |
   +----------------------------+--------------------------------------------------------------------------------------+
   | ``calcgrid_vdw1.prm``      | vdW term only (``ECUT = 1``), for calculating vdW grid (used by ``rbcalcgrid``)      |
   +----------------------------+--------------------------------------------------------------------------------------+
   | ``calcgrid_vdw5.prm``      | vdW term only (``ECUT = 5``), for calculating vdW grid (used by ``rbcalcgrid``)      |
   +----------------------------+--------------------------------------------------------------------------------------+
   | ``Tripos52_vdw.prm``       | vdW term parameter file                                                              |
   +----------------------------+--------------------------------------------------------------------------------------+
   | ``Tripos52_dihedrals.prm`` | Dihedral term parameter file                                                         |
   +----------------------------+--------------------------------------------------------------------------------------+
   | ``solvation_asp.prm``      | Desolvation term parameter file                                                      |
   +----------------------------+--------------------------------------------------------------------------------------+

.. note::

   External restraint penalty terms are defined by the user in the system
   definition ``.prm`` file. Originally, rDock did not support flexible receptor
   dihedrals or explicit structural waters, and the overall scoring function
   consisted of just the ``SCORE.INTER`` and ``SCORE.INTRA`` branches. At that
   time, the intermolecular scoring function definition file (e.g.
   ``RbtInterIdxSF.prm``) represented precisely the ``SCORE.INTER`` terms, and
   the intramolecular definition file (``RbtIntraSF.prm``) represented precisely
   the ``SCORE.INTRA`` terms. With the introduction of receptor flexibility and
   explicit structural waters (and hence the need for the ``SCORE.SYSTEM``
   branch), the situation is slightly more complex. For implementation reasons,
   many of the terms reported under ``SCORE.SYSTEM`` (with the exception of the
   dihedral term) are calculated simultaneously with the terms reported under
   ``SCORE.INTER``, and hence their parameterisation is defined implicitly in
   the intermolecular scoring function definition file. In contrast, the ligand
   intramolecular scoring function terms can be controlled independently.

References
----------

.. [ASTEX2007] Hartshorn, M.J., Verdonk, M.L., Chessari, G., Brewerton, S.C.,
	       Mooij, W.T.M., Mortenson, P.N., and Murray, C.W. (2007). Diverse,
	       High-Quality Test Set for the Validation of Protein-Ligand
	       Docking Performance. J. Med. Chem. 50, 726--741.
	       `doi:10.1021/jm061277y <https://doi.org/10.1021/jm061277y>`__
.. [GOLD2005] Verdonk, M.L., Chessari, G., Cole, J.C., Hartshorn, M.J.,
	      Murray, C.W., Nissink, J.W.M., Taylor, R.D., and
	      Taylor, R. (2005). Modeling Water Molecules in Protein-Ligand
	      Docking Using GOLD. J. Med. Chem. 48, 6504--6515.
	      `doi:10.1021/jm050543p <https://doi.org/10.1021/jm050543p>`__
.. [PDBbind2004] Wang, R., Fang, X., Lu, Y., and Wang, S. (2004). The PDBbind
		 Database: Collection of Binding Affinities for Protein-Ligand
		 Complexes with Known Three-Dimensional Structures. J. Med.
		 Chem. 47, 2977--2980. `doi:10.1021/jm030580l
		 <https://doi.org/10.1021/jm030580l>`__
.. [WSAS2001] Wang, J., Wang, W., Huo, S., Lee, M., and Kollman, P.A. (2001).
	      Solvation Model Based on Weighted Solvent Accessible Surface Area.
	      J. Phys. Chem. B 105, 5055--5067. `doi:10.1021/jp0102318
	      <https://doi.org/10.1021/jp0102318>`__
.. [CORINA1990] Gasteiger, J., Rudolph, C., and Sadowski, J. (1990). Automatic
		generation of 3D-atomic coordinates for organic molecules.
		Tetrahedron Computer Methodology 3, 537--547.
		`10.1016/0898-5529(90)90156-3
		<https://doi.org/10.1016/0898-5529(90)90156-3>`__
.. [RASASA1988] Hasel, W., Hendrickson, T.F., and Still, W.C. (1988). A rapid
		approximation to the solvent accessible surface areas of atoms.
		Tetrahedron Computer Methodology 1, 103--116.
		`doi:10.1016/0898-5529(88)90015-2
		<https://doi.org/10.1016/0898-5529(88)90015-2>`__
