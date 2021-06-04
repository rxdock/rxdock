#include "SearchTest.h"
#include "rxdock/BiMolWorkSpace.h"
#include "rxdock/CavityGridSF.h"
#include "rxdock/GATransform.h"
#include "rxdock/MdlFileSink.h"
#include "rxdock/MdlFileSource.h"
#include "rxdock/PRMFactory.h"
#include "rxdock/RandPopTransform.h"
#include "rxdock/SimAnnTransform.h"
#include "rxdock/SimplexTransform.h"
#include "rxdock/TransformAgg.h"
#include "rxdock/VdwIdxSF.h"
#include "rxdock/VdwIntraSF.h"

using namespace rxdock;
using namespace rxdock::unittest;

void SearchTest::SetUp() {
  try {
    // Create the docking site, receptor, ligand and solvent objects
    const std::string &wsName = "1YET";
    std::string prmFileName = GetDataFileName("", wsName + ".prm");
    std::string ligFileName = GetDataFileName("", wsName + "_c.sd");
    std::string dockingSiteFileName =
        GetDataFileName("", wsName + "-docking-site.json");
    ParameterFileSourcePtr spPrmSource(new ParameterFileSource(prmFileName));
    MolecularFileSourcePtr spMdlFileSource(
        new MdlFileSource(ligFileName, true, true, true));
    m_workSpace = new BiMolWorkSpace();
    std::ifstream dockingSiteFile(dockingSiteFileName.c_str());
    json siteData;
    dockingSiteFile >> siteData;
    dockingSiteFile.close();
    m_workSpace->SetDockingSite(new DockingSite(siteData.at("docking-site")));
    PRMFactory prmFactory(spPrmSource, m_workSpace->GetDockingSite());
    m_workSpace->SetReceptor(prmFactory.CreateReceptor());
    m_workSpace->SetLigand(prmFactory.CreateLigand(spMdlFileSource));
    m_workSpace->SetSolvent(prmFactory.CreateSolvent());
    // Combine the atom lists of receptor, ligand and solvent
    int nModels = m_workSpace->GetNumModels();
    for (int i = 0; i < nModels; i++) {
      AtomList atomList = m_workSpace->GetModel(i)->GetAtomList();
      std::copy(atomList.begin(), atomList.end(),
                std::back_inserter(m_atomList));
    }
    // Set up a minimal workspace and scoring function for docking
    m_SF = new SFAgg(GetMetaDataPrefix() + "score");
    BaseSF *sfInter = new VdwIdxSF("inter.vdw");
    sfInter->SetParameter(VdwSF::GetEcut(), 1.0);
    m_SF->Add(sfInter);
    BaseSF *sfIntra = new VdwIntraSF("intra.vdw");
    sfIntra->SetParameter(VdwSF::GetEcut(), 1.0);
    m_SF->Add(sfIntra);
    m_workSpace->SetSF(m_SF);
  } catch (Error &e) {
    std::cout << e.what() << std::endl;
  }
}

void SearchTest::TearDown() {
  m_atomList.clear();
  m_SF.SetNull();
  m_workSpace.SetNull();
}

// RMSD calculation between two coordinate lists
double SearchTest::rmsd(const CoordList &rc, const CoordList &c) {
  double retVal(0.0);
  unsigned int nCoords = rc.size();
  if (c.size() != nCoords) {
    retVal = 999.9;
  } else {
    for (unsigned int i = 0; i < nCoords; i++) {
      retVal += Length2(rc[i], c[i]);
    }
    retVal = std::sqrt(retVal / float(nCoords));
  }
  return retVal;
}

// 1 Check that receptor, ligand and solvent models are loaded into workspace
// Should be 6 models in total (4 solvent)
TEST_F(SearchTest, PRMFactory) { ASSERT_EQ(m_workSpace->GetNumModels(), 6); }

// 2 Check FlexDataVisitor subclass correctly identifies movable heavy atoms
// in cavity
TEST_F(SearchTest, HeavyAtomFactory) {
  AtomRList heavyAtomList;
  // find all the movable heavy atoms in the receptor, ligand and solvent
  if (m_workSpace) {
    CavityGridSF::HeavyAtomFactory atomFactory(m_workSpace->GetModels());
    heavyAtomList = atomFactory.GetAtomList();
  }
  ASSERT_EQ(heavyAtomList.size(), 42);
}

// 3 Run a sample GA
TEST_F(SearchTest, GA) {
  TransformAggPtr spTransformAgg(new TransformAgg());
  BaseTransform *pRandPop = new RandPopTransform();
  BaseTransform *pGA = new GATransform();
  spTransformAgg->Add(pRandPop);
  spTransformAgg->Add(pGA);
  m_workSpace->SetTransform(spTransformAgg);
  try {
    ASSERT_NO_THROW(m_workSpace->Run());
  } catch (Error &e) {
    std::cout << e.Message() << std::endl;
  }
}

// 4 Run a sample Simplex
TEST_F(SearchTest, Simplex) {
  TransformAggPtr spTransformAgg(new TransformAgg());
  // BaseTransform* pRandPop = new RandPopTransform();
  BaseTransform *pSimplex = new SimplexTransform();
  pSimplex->SetParameter(SimplexTransform::GetMaxCalls(), 500);
  pSimplex->SetParameter(SimplexTransform::GetNCycles(), 100);
  pSimplex->SetParameter(SimplexTransform::GetStepSize(), 1.0);
  // spTransformAgg->Add(pRandPop);
  spTransformAgg->Add(pSimplex);
  m_workSpace->SetTransform(spTransformAgg);
  try {
    ASSERT_NO_THROW(m_workSpace->Run());
  } catch (Error &e) {
    std::cout << e.Message() << std::endl;
  }
}

// 5 Run a sample simulated annealing
TEST_F(SearchTest, SimAnn) {
  BaseTransform *pSimAnn = new SimAnnTransform();
  pSimAnn->SetParameter(SimAnnTransform::GetBlockLength(), 100);
  pSimAnn->SetParameter(SimAnnTransform::GetNumBlocks(), 50);
  pSimAnn->SetParameter(SimAnnTransform::GetStartT(), 300.0);
  pSimAnn->SetParameter(SimAnnTransform::GetFinalT(), 50.0);
  pSimAnn->SetParameter(SimAnnTransform::GetStepSize(), 0.5);
  pSimAnn->SetParameter(SimAnnTransform::GetPartitionDist(), 0.0);
  pSimAnn->SetParameter(SimAnnTransform::GetPartitionFreq(), 0);
  m_workSpace->SetTransform(pSimAnn);
  try {
    ASSERT_NO_THROW(m_workSpace->Run());
  } catch (Error &e) {
    std::cout << e.Message() << std::endl;
  }
  delete pSimAnn;
}

// 6 Check we can reload solvent coords from ligand SD file
TEST_F(SearchTest, Restart) {
  TransformAggPtr spTransformAgg(new TransformAgg());
  BaseTransform *pSimplex = new SimplexTransform();
  pSimplex->SetParameter(SimplexTransform::GetMaxCalls(), 500);
  pSimplex->SetParameter(SimplexTransform::GetNCycles(), 100);
  pSimplex->SetParameter(SimplexTransform::GetStepSize(), 1.0);
  spTransformAgg->Add(pSimplex);
  m_workSpace->SetTransform(spTransformAgg);
  double finalScore(0.0);
  double restartScore(0.0);
  try {
    // create an output sink for the minimised ligand/solvent
    MolecularFileSinkPtr spMdlFileSink(
        new MdlFileSink("restart.sd", ModelPtr()));
    m_workSpace->SetSink(spMdlFileSink);
    ASSERT_NO_THROW(m_workSpace->Run());
    ASSERT_NO_THROW(m_workSpace->Save());
    finalScore = m_workSpace->GetSF()->Score();
    // Reload the receptor, minimised ligand and solvent
    std::string prmFileName = GetDataFileName("", "1YET.prm");
    ParameterFileSourcePtr spPrmSource(new ParameterFileSource(prmFileName));
    MolecularFileSourcePtr spMdlFileSource(
        new MdlFileSource("restart.sd", true, true, true));
    // Ligand segment is always called H, solvent will be H2, H3 etc.
    ASSERT_NO_THROW(
        spMdlFileSource->SetSegmentFilterMap(ConvertStringToSegmentMap("H")));
    PRMFactory prmFactory(spPrmSource, m_workSpace->GetDockingSite());
    ASSERT_NO_THROW(m_workSpace->SetReceptor(prmFactory.CreateReceptor()));
    ASSERT_NO_THROW(
        m_workSpace->SetLigand(prmFactory.CreateLigand(spMdlFileSource)));
    ASSERT_NO_THROW(m_workSpace->SetSolvent(prmFactory.CreateSolvent()));
    ASSERT_NO_THROW(
        m_workSpace->UpdateModelCoordsFromChromRecords(spMdlFileSource));
    restartScore = m_workSpace->GetSF()->Score();
    // StringVariantMap scoreMap;
    // m_workSpace->GetSF()->ScoreMap(scoreMap);
    // for (StringVariantMapConstIter iter = scoreMap.begin(); iter !=
    // scoreMap.end(); ++iter) {
    //  std::cout << iter->first << " = " << iter->second << std::endl;
    //}
    std::cout << "Final score = " << finalScore << std::endl;
    std::cout << "Restart score = " << restartScore << std::endl;
  } catch (Error &e) {
    std::cout << e.Message() << std::endl;
  }
  ASSERT_LT(std::fabs(restartScore - finalScore), 0.01);
}
