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

#include "RbtBiMolWorkSpace.h"
#include "RbtModel.h"

class OccupancyTest : public ::testing::Test {
protected:
  static double TINY;
  // TextFixture methods
  void SetUp() override;
  void TearDown() override;

  // Helper functions
  double CompareScoresForDisabledAndNoSolvent();
  double CompareScoresForSolventModes();
  bool CheckFlexAtomFactory(RbtModel *pModel, int expected[9][3]);
  RbtBiMolWorkSpacePtr m_workSpace;
  RbtModelList m_solventList;
};

#endif /*OCCUPANCY_TEST_H_*/
