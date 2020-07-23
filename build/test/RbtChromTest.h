// Unit tests for RbtChromElement and subclasses
//
// Required input files:
// 1koc.prm  rDock receptor file
// 1koc.psf  Receptor topology file
// 1koc.crd  Receptor coordinate file
// 1koc_c.sd Ligand coordinate file
// 1koc.as   Docking site
//
// Required environment:
// Make sure the above files are colocated in a single directory
// and define RBT_HOME env. variable to point at this directory
#ifndef RBTCHROMTEST_H_
#define RBTCHROMTEST_H_

#include <gtest/gtest.h>

#include "RbtBiMolWorkSpace.h"
#include "RbtChrom.h"
#include "RbtDockingSite.h"
#include "RbtModel.h"
#include "RbtSFAgg.h"

namespace rxdock {

namespace unittest {

class RbtChromTest : public ::testing::Test {
protected:
  static double TINY;
  // TextFixture methods
  void SetUp() override;
  void TearDown() override;

  // extra setup for building a simple scoring function and workspace
  void setupWorkSpace();

  // rdock helper methods
  // RMSD calculation between two coordinate lists
  double rmsd(const RbtCoordList &rc, const RbtCoordList &c);
  // Measures the variation in chromosome vector values for repeated
  // mutations or randomisations
  void measureRandOrMutateDiff(RbtChromElement *chrom, int nTrials,
                               bool bMutate, // true=mutate; false=randomise
                               double &meanDiff, double &minDiff,
                               double &maxDiff);
  // Measures the variation in chromosome vector values for repeated
  // crossovers
  void measureCrossoverDiff(RbtChromElement *chrom, int nTrials,
                            double &meanDiff, double &minDiff, double &maxDiff);

  RbtModelPtr m_recep_1koc;
  RbtModelPtr m_lig_1koc;
  RbtDockingSitePtr m_site_1koc;
  RbtChromElementPtr m_chrom_1koc;
  RbtAtomList m_atomList;           // All atoms in receptor and ligand
  RbtBiMolWorkSpacePtr m_workSpace; // simple workspace
  RbtSFAggPtr m_SF;                 // simple scoring function
};

} // namespace unittest

} // namespace rxdock

#endif /*RBTCHROMTEST_H_*/
