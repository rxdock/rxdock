#include "SearchTest.h"
#include "RbtBiMolWorkSpace.h"
#include "RbtCavityGridSF.h"
#include "RbtGATransform.h"
#include "RbtMdlFileSink.h"
#include "RbtMdlFileSource.h"
#include "RbtPRMFactory.h"
#include "RbtRandPopTransform.h"
#include "RbtSimAnnTransform.h"
#include "RbtSimplexTransform.h"
#include "RbtTransformAgg.h"
#include "RbtVdwIdxSF.h"
#include "RbtVdwIntraSF.h"

void SearchTest::SetUp() {
  try {
    // Create the docking site, receptor, ligand and solvent objects
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
    // Combine the atom lists of receptor, ligand and solvent
    int nModels = m_workSpace->GetNumModels();
    for (int i = 0; i < nModels; i++) {
      RbtAtomList atomList = m_workSpace->GetModel(i)->GetAtomList();
      std::copy(atomList.begin(), atomList.end(),
                std::back_inserter(m_atomList));
    }
    // Set up a minimal workspace and scoring function for docking
    m_SF = new RbtSFAgg("SCORE");
    RbtBaseSF *sfInter = new RbtVdwIdxSF("INTER_VDW");
    sfInter->SetParameter(RbtVdwSF::GetEcut(), 1.0);
    m_SF->Add(sfInter);
    RbtBaseSF *sfIntra = new RbtVdwIntraSF("INTRA_VDW");
    sfIntra->SetParameter(RbtVdwSF::GetEcut(), 1.0);
    m_SF->Add(sfIntra);
    m_workSpace->SetSF(m_SF);
  } catch (RbtError &e) {
    std::cout << e << std::endl;
  }
}

void SearchTest::TearDown() {
  m_atomList.clear();
  m_SF.SetNull();
  m_workSpace.SetNull();
}

// RMSD calculation between two coordinate lists
double SearchTest::rmsd(const RbtCoordList &rc, const RbtCoordList &c) {
  double retVal(0.0);
  unsigned int nCoords = rc.size();
  if (c.size() != nCoords) {
    retVal = 999.9;
  } else {
    for (unsigned int i = 0; i < nCoords; i++) {
      retVal += Rbt::Length2(rc[i], c[i]);
    }
    retVal = std::sqrt(retVal / float(nCoords));
  }
  return retVal;
}

// 1 Check that receptor, ligand and solvent models are loaded into workspace
// Should be 6 models in total (4 solvent)
TEST_F(SearchTest, PRMFactory) { ASSERT_EQ(m_workSpace->GetNumModels(), 6); }

// 2 Check RbtFlexDataVisitor subclass correctly identifies movable heavy atoms
// in cavity
TEST_F(SearchTest, HeavyAtomFactory) {
  RbtAtomRList heavyAtomList;
  // find all the movable heavy atoms in the receptor, ligand and solvent
  if (m_workSpace) {
    RbtCavityGridSF::HeavyAtomFactory atomFactory(m_workSpace->GetModels());
    heavyAtomList = atomFactory.GetAtomList();
  }
  ASSERT_EQ(heavyAtomList.size(), 42);
}

// 3 Run a sample GA
TEST_F(SearchTest, GA) {
  RbtTransformAggPtr spTransformAgg(new RbtTransformAgg());
  RbtBaseTransform *pRandPop = new RbtRandPopTransform();
  RbtBaseTransform *pGA = new RbtGATransform();
  pRandPop->SetParameter(RbtBaseObject::GetTraceStr(), 4);
  pGA->SetParameter(RbtBaseObject::GetTraceStr(), 4);
  spTransformAgg->Add(pRandPop);
  spTransformAgg->Add(pGA);
  m_workSpace->SetTransform(spTransformAgg);
  bool isOK(true);
  try {
    m_workSpace->Run();
  } catch (RbtError &e) {
    std::cout << e.Message() << std::endl;
    isOK = false;
  }
  ASSERT_TRUE(isOK);
}

// 4 Run a sample Simplex
TEST_F(SearchTest, Simplex) {
  RbtTransformAggPtr spTransformAgg(new RbtTransformAgg());
  // RbtBaseTransform* pRandPop = new RbtRandPopTransform();
  RbtBaseTransform *pSimplex = new RbtSimplexTransform();
  // pRandPop->SetParameter(RbtBaseObject::GetTraceStr(), 4);
  pSimplex->SetParameter(RbtBaseObject::GetTraceStr(), 1);
  pSimplex->SetParameter(RbtSimplexTransform::GetMaxCalls(), 500);
  pSimplex->SetParameter(RbtSimplexTransform::GetNCycles(), 100);
  pSimplex->SetParameter(RbtSimplexTransform::GetStepSize(), 1.0);
  // spTransformAgg->Add(pRandPop);
  spTransformAgg->Add(pSimplex);
  m_workSpace->SetTransform(spTransformAgg);
  bool isOK(true);
  try {
    m_workSpace->Run();
  } catch (RbtError &e) {
    std::cout << e.Message() << std::endl;
    isOK = false;
  }
  ASSERT_TRUE(isOK);
}

// 5 Run a sample simulated annealing
TEST_F(SearchTest, SimAnn) {
  RbtBaseTransform *pSimAnn = new RbtSimAnnTransform();
  pSimAnn->SetParameter(RbtBaseObject::GetTraceStr(), 2);
  pSimAnn->SetParameter(RbtSimAnnTransform::GetBlockLength(), 100);
  pSimAnn->SetParameter(RbtSimAnnTransform::GetNumBlocks(), 50);
  pSimAnn->SetParameter(RbtSimAnnTransform::GetStartT(), 300.0);
  pSimAnn->SetParameter(RbtSimAnnTransform::GetFinalT(), 50.0);
  pSimAnn->SetParameter(RbtSimAnnTransform::GetStepSize(), 0.5);
  pSimAnn->SetParameter(RbtSimAnnTransform::GetPartitionDist(), 0.0);
  pSimAnn->SetParameter(RbtSimAnnTransform::GetPartitionFreq(), 0);
  m_workSpace->SetTransform(pSimAnn);
  bool isOK(true);
  try {
    m_workSpace->Run();
  } catch (RbtError &e) {
    std::cout << e.Message() << std::endl;
    isOK = false;
  }
  delete pSimAnn;
  ASSERT_TRUE(isOK);
}

// 6 Check we can reload solvent coords from ligand SD file
TEST_F(SearchTest, Restart) {
  RbtTransformAggPtr spTransformAgg(new RbtTransformAgg());
  RbtBaseTransform *pSimplex = new RbtSimplexTransform();
  pSimplex->SetParameter(RbtBaseObject::GetTraceStr(), 1);
  pSimplex->SetParameter(RbtSimplexTransform::GetMaxCalls(), 500);
  pSimplex->SetParameter(RbtSimplexTransform::GetNCycles(), 100);
  pSimplex->SetParameter(RbtSimplexTransform::GetStepSize(), 1.0);
  spTransformAgg->Add(pSimplex);
  m_workSpace->SetTransform(spTransformAgg);
  bool isOK(true);
  double finalScore(0.0);
  double restartScore(0.0);
  try {
    // create an output sink for the minimised ligand/solvent
    RbtMolecularFileSinkPtr spMdlFileSink(
        new RbtMdlFileSink("restart.sd", RbtModelPtr()));
    m_workSpace->SetSink(spMdlFileSink);
    m_workSpace->Run();
    m_workSpace->Save();
    finalScore = m_workSpace->GetSF()->Score();
    // Reload the receptor, minimised ligand and solvent
    std::string prmFileName = Rbt::GetRbtFileName("", "1YET.prm");
    RbtParameterFileSourcePtr spPrmSource(
        new RbtParameterFileSource(prmFileName));
    RbtMolecularFileSourcePtr spMdlFileSource(
        new RbtMdlFileSource("restart.sd", true, true, true));
    // Ligand segment is always called H, solvent will be H2, H3 etc.
    spMdlFileSource->SetSegmentFilterMap(Rbt::ConvertStringToSegmentMap("H"));
    RbtPRMFactory prmFactory(spPrmSource, m_workSpace->GetDockingSite());
    prmFactory.SetTrace(1);
    m_workSpace->SetReceptor(prmFactory.CreateReceptor());
    m_workSpace->SetLigand(prmFactory.CreateLigand(spMdlFileSource));
    m_workSpace->SetSolvent(prmFactory.CreateSolvent());
    m_workSpace->UpdateModelCoordsFromChromRecords(spMdlFileSource, 1);
    restartScore = m_workSpace->GetSF()->Score();
    // RbtStringVariantMap scoreMap;
    // m_workSpace->GetSF()->ScoreMap(scoreMap);
    // for (RbtStringVariantMapConstIter iter = scoreMap.begin(); iter !=
    // scoreMap.end(); ++iter) {
    //  std::cout << iter->first << " = " << iter->second << std::endl;
    //}
    std::cout << "Final score = " << finalScore << std::endl;
    std::cout << "Restart score = " << restartScore << std::endl;
  } catch (RbtError &e) {
    std::cout << e.Message() << std::endl;
    isOK = false;
  }
  ASSERT_TRUE(isOK);
  ASSERT_LT(std::fabs(restartScore - finalScore), 0.01);
}
