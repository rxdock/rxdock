// Unit tests for rDock search algorithms
//
// Required input files:
// 1YET.prm  rDock receptor file
// 1YET.psf  Receptor topology file
// 1YET.crd  Receptor coordinate file
// 1YET_c.sd Ligand coordinate file
// 1YET.as   Docking site
//
// Required environment:
// Make sure the above files are colocated in a single directory
// and define RBT_HOME env. variable to point at this directory
#ifndef SEARCHTEST_H_
#define SEARCHTEST_H_

#include <gtest/gtest.h>

#include "RbtBiMolWorkSpace.h"
#include "RbtChrom.h"
#include "RbtDockingSite.h"
#include "RbtModel.h"
#include "RbtSFAgg.h"

namespace rxdock {

namespace unittest {

class SearchTest : public ::testing::Test {
protected:
  // TextFixture methods
  void SetUp() override;
  void TearDown() override;

  // rdock helper methods
  // RMSD calculation between two coordinate lists
  double rmsd(const RbtCoordList &rc, const RbtCoordList &c);

  RbtAtomList m_atomList;           // All atoms in receptor, ligand and solvent
  RbtBiMolWorkSpacePtr m_workSpace; // simple workspace
  RbtSFAggPtr m_SF;                 // simple scoring function
};

} // namespace unittest

} // namespace rxdock

#endif /*SEARCHTEST_H_*/
