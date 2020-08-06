// Unit tests for solvent occupancy
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
#ifndef OCCUPANCY_H_
#define OCCUPANCY_H_

#include <gtest/gtest.h>

#include "rxdock/BiMolWorkSpace.h"
#include "rxdock/Model.h"

namespace rxdock {

namespace unittest {

class OccupancyTest : public ::testing::Test {
protected:
  static double TINY;
  // TextFixture methods
  void SetUp() override;
  void TearDown() override;

  // Helper functions
  double CompareScoresForDisabledAndNoSolvent();
  double CompareScoresForSolventModes();
  bool CheckFlexAtomFactory(Model *pModel, int expected[9][3]);
  BiMolWorkSpacePtr m_workSpace;
  ModelList m_solventList;
};

} // namespace unittest

} // namespace rxdock

#endif /*OCCUPANCY_TEST_H_*/
