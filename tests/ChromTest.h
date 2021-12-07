// Unit tests for ChromElement and subclasses
//
// Required input files:
// 1koc.json RxDock receptor file
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

#include "rxdock/BiMolWorkSpace.h"
#include "rxdock/Chrom.h"
#include "rxdock/DockingSite.h"
#include "rxdock/Model.h"
#include "rxdock/SFAgg.h"

namespace rxdock {

namespace unittest {

class ChromTest : public ::testing::Test {
protected:
  static double TINY;
  // TextFixture methods
  void SetUp() override;
  void TearDown() override;

  // extra setup for building a simple scoring function and workspace
  void setupWorkSpace();

  // rdock helper methods
  // RMSD calculation between two coordinate lists
  double rmsd(const CoordList &rc, const CoordList &c);
  // Measures the variation in chromosome vector values for repeated
  // mutations or randomisations
  void measureRandOrMutateDiff(ChromElement *chrom, int nTrials,
                               bool bMutate, // true=mutate; false=randomise
                               double &meanDiff, double &minDiff,
                               double &maxDiff);
  // Measures the variation in chromosome vector values for repeated
  // crossovers
  void measureCrossoverDiff(ChromElement *chrom, int nTrials, double &meanDiff,
                            double &minDiff, double &maxDiff);

  ModelPtr m_recep_1koc;
  ModelPtr m_lig_1koc;
  DockingSitePtr m_site_1koc;
  ChromElementPtr m_chrom_1koc;
  AtomList m_atomList;           // All atoms in receptor and ligand
  BiMolWorkSpacePtr m_workSpace; // simple workspace
  SFAggPtr m_SF;                 // simple scoring function
};

} // namespace unittest

} // namespace rxdock

#endif /*RBTCHROMTEST_H_*/
