#include "OccupancyTest.h"
#include "RbtFlexAtomFactory.h"
#include "RbtLigandFlexData.h"
#include "RbtMdlFileSink.h"
#include "RbtMdlFileSource.h"
#include "RbtPRMFactory.h"
#include "RbtPolarIdxSF.h"
#include "RbtSAIdxSF.h"
#include "RbtSFAgg.h"
#include "RbtSetupPolarSF.h"
#include "RbtVariant.h"
#include "RbtVdwIdxSF.h"

double OccupancyTest::TINY = 1E-4;

void OccupancyTest::SetUp() {
  try {
    // Create the docking site, receptor, and ligand objects
    const std::string &wsName = "1YET";
    std::string prmFileName = Rbt::GetRbtFileName("", wsName + ".prm");
    std::string ligFileName = Rbt::GetRbtFileName("", wsName + "_c.sd");
    std::string asFileName = Rbt::GetRbtFileName("", wsName + ".as");
    RbtParameterFileSourcePtr spPrmSource(
        new RbtParameterFileSource(prmFileName));
    RbtMolecularFileSourcePtr spMdlFileSource(
        new RbtMdlFileSource(ligFileName, true, true, true));
    m_workSpace = new RbtBiMolWorkSpace();
    std::ifstream istr(asFileName.c_str(),
                       std::ios_base::in | std::ios_base::binary);
    m_workSpace->SetDockingSite(new RbtDockingSite(istr));
    istr.close();
    RbtPRMFactory prmFactory(spPrmSource, m_workSpace->GetDockingSite());
    m_workSpace->SetReceptor(prmFactory.CreateReceptor());
    m_workSpace->SetLigand(prmFactory.CreateLigand(spMdlFileSource));
    m_workSpace->SetSolvent(prmFactory.CreateSolvent());
    m_solventList = m_workSpace->GetSolvent();
  } catch (RbtError &e) {
    std::cout << e << std::endl;
  }
}

void OccupancyTest::TearDown() {
  m_workSpace.SetNull();
  m_solventList.clear();
}

// 1) Check that the total desolvation score for receptor / ligand / (disabled)
// solvent matches the total score with no solvent present
TEST_F(OccupancyTest, SolvationSF) {
  RbtBaseSF *sfAgg = new RbtSFAgg("SCORE");
  RbtBaseSF *sfInter = new RbtSFAgg("INTER");
  RbtBaseSF *sfSolv = new RbtSAIdxSF("SOLV");
  sfAgg->Add(sfInter);
  sfInter->Add(sfSolv);
  m_workSpace->SetSF(sfAgg);

  ASSERT_LT(CompareScoresForDisabledAndNoSolvent(), TINY);
}

// 2) Check that the total polar score for receptor / ligand / (disabled)
// solvent matches the total score with no solvent present
TEST_F(OccupancyTest, PolarSF) {
  RbtBaseSF *sfAgg = new RbtSFAgg("SCORE");
  RbtBaseSF *sfInter = new RbtSFAgg("INTER");
  RbtBaseSF *sfSetupPolar = new RbtSetupPolarSF("SETUP");
  RbtBaseSF *sfPolar = new RbtPolarIdxSF("POLAR");
  sfAgg->Add(sfInter);
  // For this term we have to set the interaction range and increment manually
  // These values match those in the standard scoring function
  sfPolar->SetRange(5.31);
  sfPolar->SetParameter(RbtPolarIdxSF::GetIncr(), 3.36);
  sfInter->Add(sfSetupPolar);
  sfInter->Add(sfPolar);
  m_workSpace->SetSF(sfAgg);

  ASSERT_LT(CompareScoresForDisabledAndNoSolvent(), TINY);
}

// 3) Check that the total vdW score for receptor / ligand / (disabled) solvent
// matches the total score with no solvent present
TEST_F(OccupancyTest, VdwSF) {
  RbtBaseSF *sfAgg = new RbtSFAgg("SCORE");
  RbtBaseSF *sfInter = new RbtSFAgg("INTER");
  RbtBaseSF *sfVdw = new RbtVdwIdxSF("VDW");
  sfAgg->Add(sfInter);
  sfInter->Add(sfVdw);
  m_workSpace->SetSF(sfAgg);

  ASSERT_LT(CompareScoresForDisabledAndNoSolvent(), TINY);
}

TEST_F(OccupancyTest, VdwSFSolventModes) {
  RbtBaseSF *sfAgg = new RbtSFAgg("SCORE");
  RbtBaseSF *sfInter = new RbtSFAgg("INTER");
  RbtBaseSF *sfSystem = new RbtSFAgg("SYSTEM");
  RbtBaseSF *sfVdw = new RbtVdwIdxSF("VDW");
  sfVdw->SetTrace(0);
  sfAgg->Add(sfInter);
  sfAgg->Add(sfSystem);
  sfInter->Add(sfVdw);
  m_workSpace->SetSF(sfAgg);

  ASSERT_LT(CompareScoresForSolventModes(), TINY);
}

// 5) Checks the results of RbtFlexAtomFactory for receptor flexibility modes
TEST_F(OccupancyTest, FlexAtomFactoryReceptor) {
  int expected[2][3] = {{2034, 10, 0}, {2044, 0, 0}};
  RbtModelPtr spReceptor = m_workSpace->GetReceptor();
  bool isTestOK = true;
  // Test two modes:
  // mode 0 is as read from PRM file (flex OH/NH3+)
  // mode 1 is fixed receptor
  for (int mode = 0; mode < 2; ++mode) {
    if (mode == 1) {
      spReceptor->SetFlexData(nullptr);
    }
    RbtFlexAtomFactory flexAtomFactory(spReceptor);
    RbtAtomRList fixedAtomList = flexAtomFactory.GetFixedAtomList();
    RbtAtomRList tetheredAtomList = flexAtomFactory.GetTetheredAtomList();
    RbtAtomRList freeAtomList = flexAtomFactory.GetFreeAtomList();
    int nFixed = fixedAtomList.size();
    int nTethered = tetheredAtomList.size();
    int nFree = freeAtomList.size();
    // std::cout << "Receptor: #fixed = " << nFixed << "; #tethered = " <<
    // nTethered
    // << "; #free = " << nFree << std::endl; for (RbtAtomRListConstIter iter =
    // tetheredAtomList.begin(); iter != tetheredAtomList.end(); ++iter) {
    //  std::cout << "Tethered: " << (*iter)->GetAtomName() << ": " <<
    //  (*iter)->GetUser1Value() << std::endl;
    //}
    if ((nFixed != expected[mode][0]) || (nTethered != expected[mode][1]) ||
        (nFree != expected[mode][2])) {
      isTestOK = false;
    }
  }
  ASSERT_TRUE(isTestOK);
}

// 6) Checks the results of RbtFlexAtomFactory for solvent flexibility modes
TEST_F(OccupancyTest, FlexAtomFactorySolvent) {
  // Expected sizes for fixed, tethered and free atom lists
  // for each value of model flexibility mode
  int expected[9][3] = {{3, 0, 0}, {0, 3, 0}, {0, 3, 0}, {0, 3, 0}, {0, 3, 0},
                        {0, 3, 0}, {0, 0, 3}, {0, 0, 3}, {0, 0, 3}};
  bool isTestOK = false;
  RbtModelList solventList = m_workSpace->GetSolvent();
  if (!solventList.empty()) {
    RbtModelPtr solvent = solventList.front();
    isTestOK = CheckFlexAtomFactory(solventList.front(), expected);
  }
  ASSERT_TRUE(isTestOK);
}

// 7) Checks the results of RbtFlexAtomFactory for ligand flexibility modes
TEST_F(OccupancyTest, FlexAtomFactoryLigand) {
  // Expected sizes for fixed, tethered and free atom lists
  // for each value of model flexibility mode
  int expected[9][3] = {{0, 0, 44}, {0, 0, 44}, {0, 0, 44},
                        {0, 0, 44}, {0, 0, 44}, {0, 0, 44},
                        {0, 0, 44}, {0, 0, 44}, {0, 0, 44}};
  bool isTestOK = false;
  RbtModelPtr spLigand = m_workSpace->GetLigand();
  if (!spLigand.Null()) {
    isTestOK = CheckFlexAtomFactory(spLigand, expected);
  }
  ASSERT_TRUE(isTestOK);
}

double OccupancyTest::CompareScoresForDisabledAndNoSolvent() {
  double retVal = 0.0;
  RbtBaseSF *pSF = m_workSpace->GetSF();
  if (pSF) {
    // A) No solvent present
    m_workSpace->RemoveSolvent();
    double scoreNoSolvent = pSF->Score();
    // B) With disabled solvent present
    for (RbtModelListIter iter = m_solventList.begin();
         iter != m_solventList.end(); ++iter) {
      (*iter)->SetOccupancy(0.0, 0.5);
    }
    m_workSpace->SetSolvent(m_solventList);
    double scoreDisabledSolvent = pSF->Score();
    // C) With enabled solvent present (not used at present)
    for (RbtModelListIter iter = m_solventList.begin();
         iter != m_solventList.end(); ++iter) {
      (*iter)->SetOccupancy(1.0, 0.5);
    }
    m_workSpace->SetSolvent(m_solventList);
    double scoreEnabledSolvent = pSF->Score();
    std::cout << "Score no solvent = " << scoreNoSolvent << std::endl;
    std::cout << "Score enabled solvent = " << scoreEnabledSolvent << std::endl;
    std::cout << "Score disabled solvent = " << scoreDisabledSolvent
              << std::endl;
    retVal = std::fabs(scoreNoSolvent - scoreDisabledSolvent);
  }
  return retVal;
}

bool OccupancyTest::CheckFlexAtomFactory(RbtModel *pModel, int expected[9][3]) {
  bool retVal = true;
  if (pModel) {
    RbtFlexData *pFlexData = pModel->GetFlexData();
    RbtFlexAtomFactory flexAtomFactory;
    for (int mode = 0; mode < 9; ++mode) {
      RbtChromElement::eMode eTransMode = (RbtChromElement::eMode)(mode / 3);
      RbtChromElement::eMode eRotMode = (RbtChromElement::eMode)(mode % 3);
      std::string strTransMode = RbtChromElement::ModeToStr(eTransMode);
      std::string strRotMode = RbtChromElement::ModeToStr(eRotMode);
      pFlexData->SetParameter(RbtLigandFlexData::GetTransMode(), strTransMode);
      pFlexData->SetParameter(RbtLigandFlexData::GetRotMode(), strRotMode);
      pModel->SetFlexData(pFlexData);
      RbtFlexAtomFactory flexAtomFactory(pModel);
      RbtAtomRList fixedAtomList = flexAtomFactory.GetFixedAtomList();
      RbtAtomRList tetheredAtomList = flexAtomFactory.GetTetheredAtomList();
      RbtAtomRList freeAtomList = flexAtomFactory.GetFreeAtomList();
      int nFixed = fixedAtomList.size();
      int nTethered = tetheredAtomList.size();
      int nFree = freeAtomList.size();
      // std::cout << std::endl << "Trans = " << strTransMode << ", Rot = " <<
      // strRotMode
      // << std::endl; for (RbtAtomRListConstIter iter = fixedAtomList.begin();
      // iter
      // != fixedAtomList.end(); ++iter) {
      //  std::cout << "Fixed: " << (*iter)->GetAtomName() << ": " <<
      //  (*iter)->GetUser1Value() << std::endl;
      //}
      // for (RbtAtomRListConstIter iter = tetheredAtomList.begin(); iter !=
      // tetheredAtomList.end(); ++iter) {
      //  std::cout << "Tethered: " << (*iter)->GetAtomName() << ": " <<
      //  (*iter)->GetUser1Value() << std::endl;
      //}
      if ((nFixed != expected[mode][0]) || (nTethered != expected[mode][1]) ||
          (nFree != expected[mode][2])) {
        retVal = false;
      }
    }
  }
  return retVal;
}

double OccupancyTest::CompareScoresForSolventModes() {
  double retVal = 0.0;
  double lastScore = 0.0;
  RbtBaseSF *pSF = m_workSpace->GetSF();
  if (pSF) {
    RbtModelList solventList = m_workSpace->GetSolvent();
    for (int mode = 0; mode < 9; ++mode) {
      RbtChromElement::eMode eTransMode = (RbtChromElement::eMode)(mode / 3);
      RbtChromElement::eMode eRotMode = (RbtChromElement::eMode)(mode % 3);
      std::string strTransMode = RbtChromElement::ModeToStr(eTransMode);
      std::string strRotMode = RbtChromElement::ModeToStr(eRotMode);
      // Change the modes of all except the first solvent model, to ensure a mix
      // of tethered / free solvent
      for (RbtModelListIter iter = solventList.begin() + 1;
           iter != solventList.end(); ++iter) {
        RbtFlexData *pFlexData = (*iter)->GetFlexData();
        pFlexData->SetParameter(RbtLigandFlexData::GetTransMode(),
                                strTransMode);
        pFlexData->SetParameter(RbtLigandFlexData::GetRotMode(), strRotMode);
        (*iter)->SetFlexData(pFlexData);
      }
      m_workSpace->RemoveSolvent();
      m_workSpace->SetSolvent(solventList);
      double score = pSF->Score();
      // std::cout << "mode = " << mode << "; score = " << score << std::endl;
      if (mode > 0) {
        double diff = std::fabs(score - lastScore);
        retVal = std::max(diff, retVal);
      }
      lastScore = score;
      // RbtStringVariantMap scoreMap;
      // pSF->ScoreMap(scoreMap);
      // for (RbtStringVariantMapConstIter vIter = scoreMap.begin(); vIter !=
      // scoreMap.end(); vIter++) { std::cout << (*vIter).first << " = " <<
      // (*vIter).second << std::endl;
      //}
    }
  }
  return retVal;
}
