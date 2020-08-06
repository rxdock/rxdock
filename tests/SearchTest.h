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

#include "rxdock/BiMolWorkSpace.h"
#include "rxdock/Chrom.h"
#include "rxdock/DockingSite.h"
#include "rxdock/Model.h"
#include "rxdock/SFAgg.h"

namespace rxdock {

namespace unittest {

class SearchTest : public ::testing::Test {
protected:
  // TextFixture methods
  void SetUp() override;
  void TearDown() override;

  // rdock helper methods
  // RMSD calculation between two coordinate lists
  double rmsd(const CoordList &rc, const CoordList &c);

  AtomList m_atomList;           // All atoms in receptor, ligand and solvent
  BiMolWorkSpacePtr m_workSpace; // simple workspace
  SFAggPtr m_SF;                 // simple scoring function
};

} // namespace unittest

} // namespace rxdock

#endif /*SEARCHTEST_H_*/
