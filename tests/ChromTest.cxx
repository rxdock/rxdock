#include "ChromTest.h"
#include "rxdock/BiMolWorkSpace.h"
#include "rxdock/ChromDihedralElement.h"
#include "rxdock/ChromFactory.h"
#include "rxdock/ChromOccupancyElement.h"
#include "rxdock/ChromPositionElement.h"
#include "rxdock/Euler.h"
#include "rxdock/LigandFlexData.h"
#include "rxdock/MdlFileSource.h"
#include "rxdock/PRMFactory.h"
#include "rxdock/Rand.h"
#include "rxdock/ReceptorFlexData.h"
#include "rxdock/SFAgg.h"
#include "rxdock/VdwIdxSF.h"
#include "rxdock/VdwIntraSF.h"

using namespace rxdock;
using namespace rxdock::unittest;

double ChromTest::TINY = 1E-4;

void ChromTest::SetUp() {
  try {
    // Create a receptor model, ligand model, and docking site
    const std::string &wsName = "1koc";
    std::string prmFileName =
        GetDataFileName("data/receptors", wsName + ".json");
    std::string ligFileName = GetDataFileName("", wsName + "_c.sd");
    std::string dockingSiteFileName =
        GetDataFileName("data/grids", wsName + "-docking-site.json");
    ParameterFileSourcePtr spPrmSource(new ParameterFileSource(prmFileName));
    MolecularFileSourcePtr spMdlFileSource(
        new MdlFileSource(ligFileName, true, true, true));
    std::ifstream dockingSiteFile(dockingSiteFileName);
    json siteData;
    dockingSiteFile >> siteData;
    dockingSiteFile.close();
    m_site_1koc = new DockingSite(siteData.at("docking-site"));
    PRMFactory prmFactory(spPrmSource);
    m_recep_1koc = prmFactory.CreateReceptor();
    m_lig_1koc = prmFactory.CreateLigand(spMdlFileSource);
    // Create a ligand chromosome for free docking
    m_chrom_1koc = new Chrom();
    FlexDataPtr ligFlexData = new LigandFlexData(m_site_1koc);
    ligFlexData->SetModel(m_lig_1koc);
    ChromFactory chromFactory1;
    ligFlexData->Accept(chromFactory1);
    m_chrom_1koc->Add(chromFactory1.GetChrom());
    // Create a receptor chromosome for flexible terminal bonds
    FlexDataPtr recFlexData = new ReceptorFlexData(m_site_1koc);
    recFlexData->SetModel(m_recep_1koc);
    ChromFactory chromFactory2;
    recFlexData->Accept(chromFactory2);
    m_chrom_1koc->Add(chromFactory2.GetChrom());
    // Combine the atom lists of ligand and receptor
    m_atomList = m_lig_1koc->GetAtomList();
    AtomList recepAtomList = m_recep_1koc->GetAtomList();
    std::copy(recepAtomList.begin(), recepAtomList.end(),
              std::back_inserter(m_atomList));
  } catch (Error &e) {
    std::cout << e.what() << std::endl;
  }
}

void ChromTest::TearDown() {
  m_recep_1koc.SetNull();
  m_lig_1koc.SetNull();
  m_site_1koc.SetNull();
  m_chrom_1koc.SetNull();
  m_atomList.clear();
  m_SF.SetNull();
  m_workSpace.SetNull();
}

// extra setup/teardown for building a simple scoring function and workspace
void ChromTest::setupWorkSpace() {
  // Set up a minimal workspace and scoring function for docking
  m_SF = new SFAgg(GetMetaDataPrefix() + "score");
  BaseSF *sfInter = new VdwIdxSF("inter.vdw");
  sfInter->SetParameter(VdwSF::GetEcut(), 1.0);
  m_SF->Add(sfInter);
  BaseSF *sfIntra = new VdwIntraSF("intra.vdw");
  sfIntra->SetParameter(VdwSF::GetEcut(), 1.0);
  m_SF->Add(sfIntra);
  m_workSpace = new BiMolWorkSpace();
  m_workSpace->SetSF(m_SF);
  m_workSpace->SetDockingSite(m_site_1koc);
  m_workSpace->SetReceptor(m_recep_1koc);
  m_workSpace->SetLigand(m_lig_1koc);
}

// RMSD calculation between two coordinate lists
double ChromTest::rmsd(const CoordList &rc, const CoordList &c) {
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

// The unit tests
// 1) Does ChromFactory create a valid chromosome?
TEST_F(ChromTest, ChromFactory) { ASSERT_NE(m_chrom_1koc.Ptr(), nullptr); }

// 2) Can we delete a cloned chromosome without affecting the original?
TEST_F(ChromTest, CloneDestructor) {
  ChromElementPtr clone = m_chrom_1koc->clone();
  ASSERT_NE(m_chrom_1koc.Ptr(), nullptr);
}

// 3) Is the length of a cloned chromosome equal to the original?
TEST_F(ChromTest, CloneLengthEquals) {
  ChromElementPtr clone = m_chrom_1koc->clone();
  ASSERT_EQ(m_chrom_1koc->GetLength(), clone->GetLength());
}

// 4) Does operator== return true for a clone and an original (before mutation)?
TEST_F(ChromTest, ChromOperatorEquals) {
  ChromElementPtr clone = m_chrom_1koc->clone();
  ASSERT_EQ(*m_chrom_1koc, *clone);
}

// 5) Does operator!= return true for a mutated clone and an original?
TEST_F(ChromTest, ChromOperatorEqualsAfterMutate) {
  ChromElementPtr clone = m_chrom_1koc->clone();
  clone->Mutate(1.0);
  ASSERT_NE(*m_chrom_1koc, *clone);
}

// 6) Does SyncToModel leave the model coords unchanged if chromosome has
// not been mutated?
TEST_F(ChromTest, SyncToModel) {
  CoordList coordsBefore, coordsAfter;
  GetCoordList(m_atomList, coordsBefore);
  m_chrom_1koc->SyncToModel();
  GetCoordList(m_atomList, coordsAfter);
  double rms = rmsd(coordsBefore, coordsAfter);
  ASSERT_LT(rms, TINY);
}

// 7) Does SyncToModel change the model coords after a mutation?
TEST_F(ChromTest, SyncToModelAfterMutate) {
  CoordList coordsBefore, coordsAfter;
  GetCoordList(m_atomList, coordsBefore);
  m_chrom_1koc->Mutate(1.0);
  m_chrom_1koc->SyncToModel();
  GetCoordList(m_atomList, coordsAfter);
  double rms = rmsd(coordsBefore, coordsAfter);
  ASSERT_GT(rms, TINY);
}

// 8) Does SyncToModel leave the model coords unchanged if chromosome has
// been mutated then reset?
TEST_F(ChromTest, SyncToModelAfterReset) {
  CoordList coordsBefore, coordsAfter;
  GetCoordList(m_atomList, coordsBefore);
  m_chrom_1koc->Mutate(1.0);
  m_chrom_1koc->SyncToModel();
  m_chrom_1koc->Reset();
  m_chrom_1koc->SyncToModel();
  GetCoordList(m_atomList, coordsAfter);
  double rms = rmsd(coordsBefore, coordsAfter);
  ASSERT_LT(rms, TINY);
}

// 9) If we make a clone BEFORE the original has been mutated, check that
// Reset() on the clone will reset to the original (unmutated) model coords
TEST_F(ChromTest, CloneReset) {
  CoordList coordsBefore, coordsAfter;
  GetCoordList(m_atomList, coordsBefore);
  ChromElementPtr clone = m_chrom_1koc->clone();
  clone->Mutate(1.0);
  clone->SyncToModel();
  clone->Reset();
  clone->SyncToModel();
  GetCoordList(m_atomList, coordsAfter);
  double rms = rmsd(coordsBefore, coordsAfter);
  ASSERT_LT(rms, TINY);
}

// 10) If we make a clone AFTER the original has been mutated, check that
// Reset() on the clone still resets to the original (unmutated) model coords
TEST_F(ChromTest, MutatedCloneReset) {
  CoordList coordsBefore, coordsAfter;
  GetCoordList(m_atomList, coordsBefore);
  // Mutate the original before creating the clone
  m_chrom_1koc->Mutate(1.0);
  m_chrom_1koc->SyncToModel();
  ChromElement *clone = m_chrom_1koc->clone();
  clone->Mutate(1.0);
  clone->SyncToModel();
  clone->Reset();
  clone->SyncToModel();
  GetCoordList(m_atomList, coordsAfter);
  double rms = rmsd(coordsBefore, coordsAfter);
  ASSERT_LT(rms, TINY);
  delete clone;
}

// 11) Does operator!= return true for a randomised clone and an original?
TEST_F(ChromTest, ChromOperatorEqualsAfterRandomise) {
  ChromElementPtr clone = m_chrom_1koc->clone();
  clone->Randomise();
  ASSERT_NE(*m_chrom_1koc, *clone);
}

// 12) Does SyncToModel change the model coords after a randomisation?
TEST_F(ChromTest, SyncToModelAfterRandomise) {
  CoordList coordsBefore, coordsAfter;
  GetCoordList(m_atomList, coordsBefore);
  m_chrom_1koc->Randomise();
  m_chrom_1koc->SyncToModel();
  GetCoordList(m_atomList, coordsAfter);
  double rms = rmsd(coordsBefore, coordsAfter);
  ASSERT_GT(rms, TINY);
}

// 13) Does SyncToModel leave the model coords unchanged if chromosome has
// been randomised then reset?
TEST_F(ChromTest, SyncToModelAfterRandomiseReset) {
  CoordList coordsBefore, coordsAfter;
  GetCoordList(m_atomList, coordsBefore);
  m_chrom_1koc->Randomise();
  m_chrom_1koc->SyncToModel();
  m_chrom_1koc->Reset();
  m_chrom_1koc->SyncToModel();
  GetCoordList(m_atomList, coordsAfter);
  double rms = rmsd(coordsBefore, coordsAfter);
  ASSERT_LT(rms, TINY);
}

// 14) Does operator== return true for a mutated, reset clone and an original?
TEST_F(ChromTest, ChromOperatorEqualsAfterMutateReset) {
  ChromElementPtr clone = m_chrom_1koc->clone();
  clone->Mutate(1.0);
  clone->Reset();
  ASSERT_EQ(*m_chrom_1koc, *clone);
}

// 15) Does operator== return true for a randomised, reset clone and an
// original?
TEST_F(ChromTest, ChromOperatorEqualsAfterRandomiseReset) {
  ChromElementPtr clone = m_chrom_1koc->clone();
  clone->Randomise();
  clone->Reset();
  ASSERT_EQ(*m_chrom_1koc, *clone);
}

// 16) Check that chromosomes are still equal after crossover of unmutated
// clones
TEST_F(ChromTest, NullCrossover) {
  ChromElementPtr clone1 = m_chrom_1koc->clone();
  ChromElementPtr clone2 = m_chrom_1koc->clone();
  ChromElementPtr clone3 = m_chrom_1koc->clone();
  ChromElementPtr clone4 = m_chrom_1koc->clone();
  // is equal before crossover
  ASSERT_EQ(*clone1, *clone2);
  ASSERT_EQ(*clone1, *clone3);
  ASSERT_EQ(*clone1, *clone4);
  Crossover(clone1, clone2, clone3, clone4);
  // is equal after crossover
  ASSERT_EQ(*clone1, *clone2);
  ASSERT_EQ(*clone1, *clone3);
  ASSERT_EQ(*clone1, *clone4);
}

// 17) Check that conversion from quaternion to Euler angles and back
// is a null operation
TEST_F(ChromTest, Euler) {
  ChromElementPtr clone = m_chrom_1koc->clone();
  AtomList atomList = m_lig_1koc->GetAtomList();
  PrincipalAxes prAxes = GetPrincipalAxesOfAtoms(atomList);
  PrincipalAxes cartesianAxes;
  Quat q = GetQuatFromAlignAxes(cartesianAxes, prAxes);
  Euler euler(q);
  Quat q1 = euler.ToQuat();
  ASSERT_LT(Length(q1 - q), TINY);
}

// 18) Check that repeated sync to/from model leaves the genotype and phenotype
// unchanged
TEST_F(ChromTest, RepeatedSync) {
  // Keep clone unchanged as a reference
  ChromElementPtr clone = m_chrom_1koc->clone();
  CoordList coordsBefore, coordsAfter;
  GetCoordList(m_atomList, coordsBefore);
  for (int i = 0; i < 100; i++) {
    m_chrom_1koc->SyncToModel();
    m_chrom_1koc->SyncFromModel();
  }
  GetCoordList(m_atomList, coordsAfter);
  double rms = rmsd(coordsBefore, coordsAfter);
  double cmp = m_chrom_1koc->Compare(*clone);
  ASSERT_LT(rms, TINY);
  ASSERT_LT(cmp, TINY);
}

// 19) Check that an invalid Compare() returns -1 (i.e. for unequal chromosome
// lengths)
TEST_F(ChromTest, CompareWithNullChrom) {
  ChromElementPtr nullChrom = new Chrom();
  double cmp = m_chrom_1koc->Compare(*nullChrom);
  ASSERT_LT(cmp, 0.0);
  ASSERT_NE(*m_chrom_1koc, *nullChrom);
}

// 20) Check that Compare() after a single mutation never returns greater than
// the relative mutation distance (loops over 10000 repeats)
TEST_F(ChromTest, CompareAfterMutate) {
  Rand &rand = GetRandInstance();
  for (int i = 0; i < 10000; i++) {
    ChromElementPtr clone = m_chrom_1koc->clone();
    double mutationDistance = rand.GetRandom01();
    clone->Mutate(mutationDistance);
    double cmp = m_chrom_1koc->Compare(*clone);
    ASSERT_GE(cmp, 0.0);
    ASSERT_LE(cmp, mutationDistance);
  }
}

// 21) Test operator== with modified _THRESHOLD
TEST_F(ChromTest, OperatorEqualsWithModifiedThreshold) {
  ChromElementPtr clone = m_chrom_1koc->clone();
  clone->SyncToModel();
  clone->SyncFromModel();
  double &_THRESHOLD = ChromElement::GetThreshold();
  // Very fine threshold - should fail
  _THRESHOLD = 1E-20;
  ASSERT_NE(*m_chrom_1koc, *clone);

  // Fine threshold, but should pass
  // Note: adjusted for MSVC, should be 1E-10
  _THRESHOLD = 1E-7;
  ASSERT_EQ(*m_chrom_1koc, *clone);

  // Crude threshold - should pass even after a single mutation
  _THRESHOLD = 1.0;
  clone->Mutate(0.5);
  ASSERT_EQ(*m_chrom_1koc, *clone);
}

// 22) Test population constructor
TEST_F(ChromTest, PopulationConstructor) {
  setupWorkSpace();
  int popSize = 100;
  PopulationPtr pop = new Population(m_chrom_1koc, popSize, m_SF);
  ASSERT_EQ(pop->GetMaxSize(), popSize);
  ASSERT_EQ(pop->GetActualSize(), popSize);
}

// 23) Check that population constructor with zero size throws an error
TEST_F(ChromTest, PopulationConstructorZeroSize) {
  setupWorkSpace();
  int popSize = 0;
  ASSERT_THROW(PopulationPtr pop = new Population(m_chrom_1koc, popSize, m_SF),
               Error);
}

// 24) Check that population constructor with null chromosome throws an error
TEST_F(ChromTest, PopulationConstructorNullChromosome) {
  setupWorkSpace();
  int popSize = 100;
  ASSERT_THROW(PopulationPtr pop = new Population(nullptr, popSize, m_SF),
               Error);
}

// 25) Check that population constructor with null scoring function throws an
// error
TEST_F(ChromTest, PopulationConstructorNullSF) {
  setupWorkSpace();
  int popSize = 100;
  ASSERT_THROW(PopulationPtr pop =
                   new Population(m_chrom_1koc, popSize, nullptr),
               Error);
}

// 26) Check that roulette wheel fitness values are calculated OK
TEST_F(ChromTest, PopulationRWFitness) {
  setupWorkSpace();
  int popSize = 100;
  PopulationPtr pop = new Population(m_chrom_1koc, popSize, m_SF);
  const GenomeList &genomeList = pop->GetGenomeList();
  double lastValue = 0.0;
  for (GenomeListConstIter iter = genomeList.begin(); iter != genomeList.end();
       ++iter) {
    double value = (*iter)->GetRWFitness();
    // std::cout << (*iter)->GetScore() << "\t" << (*iter)->GetRWFitness() <<
    // std::endl;
    ASSERT_GE(value, lastValue);
    lastValue = value;
  }
  ASSERT_LT(std::fabs(lastValue - 1.0), TINY);
}

// 27) Check that roulette wheel selection will select higher scoring genomes on
// average
TEST_F(ChromTest, PopulationRWSelect) {
  setupWorkSpace();
  int popSize = 100;
  PopulationPtr pop = new Population(m_chrom_1koc, popSize, m_SF);
  double popMean = pop->GetScoreMean();
  double selectionMean(0.0);
  int nToSelect(100);
  for (int i = 0; i < nToSelect; i++) {
    GenomePtr genome = pop->RouletteWheelSelect();
    double score = genome->GetScore();
    selectionMean += score;
    // std::cout << "RWSelect score=" << score << std::endl;
  }
  selectionMean /= nToSelect;
  // std::cout << "popMean = " << popMean << "; selectionMean = " <<
  // selectionMean << std::endl;
  ASSERT_GT(selectionMean, popMean);
}

// 28) Check that population GAstep does not decrease the score of the best
// genome and that population size remains constant
TEST_F(ChromTest, PopulationGAstep) {
  setupWorkSpace();
  int popSize = 100;
  int nReplicates = 50;
  int nIter = 100;
  double equalityThreshold = 1.0E-2;
  double relStepSize = 1.0;
  double pcross = 0.4;
  bool xovermut = true;
  bool cmutate = false;
  PopulationPtr pop = new Population(m_chrom_1koc, popSize, m_SF);
  double lastScore = pop->Best()->GetScore();
  bool isOK = true;
  try {
    for (int i = 0; (i < nIter) && isOK; ++i) {
      ASSERT_NO_THROW(pop->GAstep(nReplicates, relStepSize, equalityThreshold,
                                  pcross, xovermut, cmutate));
      double score = pop->Best()->GetScore();
      ASSERT_GE(score, lastScore);
      lastScore = score;
    }
  } catch (Error &e) {
    std::cout << e.Message() << std::endl;
  }
}

// 29) Checks the behaviour of Model::GetChrom if flex data has not been
// defined Should return a zero length chromosome for rigid model
TEST_F(ChromTest, ModelGetChromUndefinedFlexData) {
  ChromElementPtr chrom = new Chrom();
  chrom->Add(m_lig_1koc->GetChrom());
  chrom->Add(m_recep_1koc->GetChrom());
  ASSERT_EQ(chrom->GetLength(), 0);
}

// 30) Checks the behaviour of Model::SetFlexData(nullptr)
// Should return a zero length chromosome for rigid model
TEST_F(ChromTest, ModelGetChromNullFlexData) {
  m_lig_1koc->SetFlexData(nullptr);
  m_recep_1koc->SetFlexData(nullptr);
  ChromElementPtr chrom = new Chrom();
  chrom->Add(m_lig_1koc->GetChrom());
  chrom->Add(m_recep_1koc->GetChrom());
  ASSERT_EQ(chrom->GetLength(), 0);
}

// 31) Checks that Model::SetFlexData creates the same chromosome
// as the manual method
TEST_F(ChromTest, ModelGetChrom) {
  m_lig_1koc->SetFlexData(new LigandFlexData(m_site_1koc));
  m_recep_1koc->SetFlexData(new ReceptorFlexData(m_site_1koc));
  ChromElementPtr chrom = new Chrom();
  chrom->Add(m_lig_1koc->GetChrom());
  chrom->Add(m_recep_1koc->GetChrom());
  ASSERT_EQ(*chrom, *m_chrom_1koc);
}

// 32) Checks that the variation in randomised values remains in bounds
// for tethered orientations
TEST_F(ChromTest, RandomiseTetheredOrientation) {
  double transStepSize = 0.1;
  double rotStepSize = 10.0;
  ChromElement::eMode transMode = ChromElement::FIXED;
  ChromElement::eMode rotMode = ChromElement::TETHERED;
  double maxTrans = 1.0;
  double maxRot = 45.0;
  rotStepSize *= M_PI / 180.0;
  maxRot *= M_PI / 180.0;
  ChromElementPtr chrom = new ChromPositionElement(
      m_lig_1koc, m_site_1koc, transStepSize, rotStepSize, transMode, rotMode,
      maxTrans, maxRot);
  double meanDiff, minDiff, maxDiff;
  // Measure 1M randomisations of the tethered orientation
  measureRandOrMutateDiff(chrom, 1000000, false, meanDiff, minDiff, maxDiff);
  // The expected mean rotation should be half the maximum step size
  // as CompareVector always returns a positive difference.
  double expectedMean = 0.5 * maxRot / rotStepSize;
  // Check that min, max rotation are in range [0,maxRot]
  // and that mean rotation is within 1% of the expected mean.
  ASSERT_GE(minDiff, 0.0);
  ASSERT_LE(maxDiff, (2.0 * expectedMean));
  ASSERT_GT((meanDiff / expectedMean), 0.99);
  ASSERT_LT((meanDiff / expectedMean), 1.01);
}

// 33) Checks that the variation in randomised values remains in bounds
// for tethered COM
TEST_F(ChromTest, RandomiseTetheredCOM) {
  double transStepSize = 0.1;
  double rotStepSize = 10.0;
  ChromElement::eMode transMode = ChromElement::TETHERED;
  ChromElement::eMode rotMode = ChromElement::FIXED;
  double maxTrans = 1.0;
  double maxRot = 45.0;
  rotStepSize *= M_PI / 180.0;
  maxRot *= M_PI / 180.0;
  ChromElementPtr chrom = new ChromPositionElement(
      m_lig_1koc, m_site_1koc, transStepSize, rotStepSize, transMode, rotMode,
      maxTrans, maxRot);
  double meanDiff, minDiff, maxDiff;
  // Measure 1M randomisations of the tethered COM
  measureRandOrMutateDiff(chrom, 1000000, false, meanDiff, minDiff, maxDiff);
  // The expected mean displacement should be half the maximum step size
  // as CompareVector always returns a positive difference.
  double expectedMean = 0.5 * maxTrans / transStepSize;
  // Check that min, max displacement are in range [0,maxTrans]
  // and that mean rotation is within 1% of the expected mean.
  ASSERT_GE(minDiff, 0.0);
  ASSERT_LE(maxDiff, (2.0 * expectedMean));
  ASSERT_GT((meanDiff / expectedMean), 0.99);
  ASSERT_LT((meanDiff / expectedMean), 1.01);
}

// 34) Checks that the variation in mutated values remains in bounds
// for tethered orientations
TEST_F(ChromTest, RandomiseTetheredDihedral) {
  double stepSize = 10.0;
  ChromElement::eMode mode = ChromElement::TETHERED;
  double maxDelta = 45.0;
  BondList rotBondList =
      GetBondListWithPredicate(m_lig_1koc->GetBondList(), isBondRotatable());
  AtomList noTetheredAtoms;
  BondListConstIter iter = rotBondList.begin();
  ChromElementPtr chrom = new ChromDihedralElement(*iter, noTetheredAtoms,
                                                   stepSize, mode, maxDelta);
  double meanDiff, minDiff, maxDiff;
  // Measure 1M randomisations of the tethered dihedrals
  measureRandOrMutateDiff(chrom, 1000000, false, meanDiff, minDiff, maxDiff);
  // The expected mean displacement should be half the maximum step size
  // as CompareVector always returns a positive difference.
  double expectedMean = 0.5 * maxDelta / stepSize;
  // Check that min, max displacement are in range [0,maxTrans]
  // and that mean rotation is within 1% of the expected mean.
  ASSERT_GE(minDiff, 0.0);
  ASSERT_LE(maxDiff, (2.0 * expectedMean));
  ASSERT_GT((meanDiff / expectedMean), 0.99);
  ASSERT_LT((meanDiff / expectedMean), 1.01);
}

// 35) Checks that the variation in mutated values remains in bounds
// for tethered COM
TEST_F(ChromTest, MutateTetheredOrientation) {
  double transStepSize = 0.1;
  double rotStepSize = 10.0;
  ChromElement::eMode transMode = ChromElement::FIXED;
  ChromElement::eMode rotMode = ChromElement::TETHERED;
  double maxTrans = 1.0;
  double maxRot = 45.0;
  rotStepSize *= M_PI / 180.0;
  maxRot *= M_PI / 180.0;
  ChromElementPtr chrom = new ChromPositionElement(
      m_lig_1koc, m_site_1koc, transStepSize, rotStepSize, transMode, rotMode,
      maxTrans, maxRot);
  double meanDiff, minDiff, maxDiff;
  measureRandOrMutateDiff(chrom, 1000000, true, meanDiff, minDiff, maxDiff);
  ASSERT_GE(minDiff, 0.0);
  ASSERT_LE(maxDiff, 1.01 * (maxRot / rotStepSize));
}

// 36) Checks that the variation in crossed over values remains in bounds
// for tethered orientations
TEST_F(ChromTest, MutateTetheredCOM) {
  double transStepSize = 0.1;
  double rotStepSize = 10.0;
  ChromElement::eMode transMode = ChromElement::TETHERED;
  ChromElement::eMode rotMode = ChromElement::FIXED;
  double maxTrans = 1.0;
  double maxRot = 45.0;
  rotStepSize *= M_PI / 180.0;
  maxRot *= M_PI / 180.0;
  ChromElementPtr chrom = new ChromPositionElement(
      m_lig_1koc, m_site_1koc, transStepSize, rotStepSize, transMode, rotMode,
      maxTrans, maxRot);
  double meanDiff, minDiff, maxDiff;
  measureRandOrMutateDiff(chrom, 1000000, true, meanDiff, minDiff, maxDiff);
  ASSERT_GE(minDiff, 0.0);
  ASSERT_LE(maxDiff, 1.01 * (maxTrans / transStepSize));
}

// 37) Checks that the variation in crossed over values remains in bounds
// for tethered COM
TEST_F(ChromTest, MutateTetheredDihedral) {
  double stepSize = 10.0;
  ChromElement::eMode mode = ChromElement::TETHERED;
  double maxDelta = 45.0;
  BondList rotBondList =
      GetBondListWithPredicate(m_lig_1koc->GetBondList(), isBondRotatable());
  AtomList noTetheredAtoms;
  BondListConstIter iter = rotBondList.begin();
  ChromElementPtr chrom = new ChromDihedralElement(*iter, noTetheredAtoms,
                                                   stepSize, mode, maxDelta);
  double meanDiff, minDiff, maxDiff;
  measureRandOrMutateDiff(chrom, 1000000, true, meanDiff, minDiff, maxDiff);
  ASSERT_GE(minDiff, 0.0);
  ASSERT_LE(maxDiff, 1.01 * (maxDelta / stepSize));
}

// 38) Checks that the variation in randomised values remains in bounds
// for tethered dihedral
TEST_F(ChromTest, CrossoverTetheredOrientation) {
  double transStepSize = 0.1;
  double rotStepSize = 10.0;
  ChromElement::eMode transMode = ChromElement::FIXED;
  ChromElement::eMode rotMode = ChromElement::TETHERED;
  double maxTrans = 1.0;
  double maxRot = 45.0;
  rotStepSize *= M_PI / 180.0;
  maxRot *= M_PI / 180.0;
  ChromElementPtr chrom = new ChromPositionElement(
      m_lig_1koc, m_site_1koc, transStepSize, rotStepSize, transMode, rotMode,
      maxTrans, maxRot);
  double meanDiff, minDiff, maxDiff;
  measureCrossoverDiff(chrom, 100, meanDiff, minDiff, maxDiff);
  ASSERT_GE(minDiff, 0.0);
  ASSERT_LE(maxDiff, 1.01 * (maxRot / rotStepSize));
}

// 39) Checks that the variation in mutated values remains in bounds
// for tethered dihedral
TEST_F(ChromTest, CrossoverTetheredCOM) {
  double transStepSize = 0.1;
  double rotStepSize = 10.0;
  ChromElement::eMode transMode = ChromElement::TETHERED;
  ChromElement::eMode rotMode = ChromElement::FIXED;
  double maxTrans = 1.0;
  double maxRot = 45.0;
  rotStepSize *= M_PI / 180.0;
  maxRot *= M_PI / 180.0;
  ChromElementPtr chrom = new ChromPositionElement(
      m_lig_1koc, m_site_1koc, transStepSize, rotStepSize, transMode, rotMode,
      maxTrans, maxRot);
  double meanDiff, minDiff, maxDiff;
  measureCrossoverDiff(chrom, 100, meanDiff, minDiff, maxDiff);
  ASSERT_GE(minDiff, 0.0);
  ASSERT_LE(maxDiff, 1.01 * (maxTrans / transStepSize));
}

// 40) Checks that the variation in crossed over values remains in bounds
// for tethered dihedral
TEST_F(ChromTest, CrossoverTetheredDihedral) {
  double stepSize = 10.0;
  ChromElement::eMode mode = ChromElement::TETHERED;
  double maxDelta = 45.0;
  BondList rotBondList =
      GetBondListWithPredicate(m_lig_1koc->GetBondList(), isBondRotatable());
  AtomList noTetheredAtoms;
  BondListConstIter iter = rotBondList.begin();
  ChromElementPtr chrom = new ChromDihedralElement(*iter, noTetheredAtoms,
                                                   stepSize, mode, maxDelta);
  double meanDiff, minDiff, maxDiff;
  measureCrossoverDiff(chrom, 100, meanDiff, minDiff, maxDiff);
  ASSERT_GE(minDiff, 0.0);
  ASSERT_LE(maxDiff, 1.01 * (maxDelta / stepSize));
}

// 41) Checks that randomised occupancy is from a rectangular distribution
TEST_F(ChromTest, RandomiseOccupancy) {
  double stepSize = 0.1;
  double threshold = 0.5;
  ChromElementPtr chrom =
      new ChromOccupancyElement(m_lig_1koc, stepSize, threshold);
  double meanDiff, minDiff, maxDiff;
  measureRandOrMutateDiff(chrom, 1000000, false, meanDiff, minDiff, maxDiff);
  // Randomised occupancy should be a rectangular distribution
  // between 0 and 1, mean = 0.5.
  minDiff *= stepSize;
  meanDiff *= stepSize;
  maxDiff *= stepSize;
  ASSERT_LT(std::fabs(minDiff), 0.01);
  ASSERT_LT(std::fabs(meanDiff - 0.5), 0.01);
  ASSERT_LT(std::fabs(maxDiff - 1.0), 0.01);
}

// 42) Checks that actual occupancy probability matches the desired probability
TEST_F(ChromTest, OccupancyThreshold) {
  double stepSize = 0.1;
  double occupancyProb = 0.7;
  double threshold = 1.0 - occupancyProb;
  ChromElementPtr chrom =
      new ChromOccupancyElement(m_lig_1koc, stepSize, threshold);
  int nTrials = 1000000;
  int nEnabled = 0;
  AtomList atomList = m_lig_1koc->GetAtomList();
  Atom *pAtom0 = atomList.front();
  // Measure whether the actual probability of enabling the model atoms
  // is equal to the desired probability
  for (int i = 0; i < nTrials; i++) {
    chrom->Randomise();
    chrom->SyncToModel();
    if (pAtom0->GetEnabled()) {
      nEnabled++;
    }
  }
  double enabledProb = (double)nEnabled / (double)nTrials;
  ASSERT_LT(std::fabs(enabledProb - occupancyProb), 0.01);
}

void ChromTest::measureRandOrMutateDiff(ChromElement *chrom, int nTrials,
                                        bool bMutate, double &meanDiff,
                                        double &minDiff, double &maxDiff) {
  std::vector<double> refVec;
  chrom->GetVector(refVec);
  meanDiff = 0.0;
  minDiff = 0.0;
  maxDiff = 0.0;
  for (int iTrial = 0; iTrial < nTrials; iTrial++) {
    if (bMutate) {
      chrom->Mutate(1.0);
    } else {
      chrom->Randomise();
    }
    int i(0);
    double diff = chrom->CompareVector(refVec, i);
    meanDiff += diff;
    if (iTrial == 0) {
      minDiff = diff;
      maxDiff = diff;
    } else {
      minDiff = std::min(minDiff, diff);
      maxDiff = std::max(maxDiff, diff);
    }
  }
  if (nTrials > 0) {
    meanDiff /= nTrials;
  }
}

void ChromTest::measureCrossoverDiff(ChromElement *chrom, int nTrials,
                                     double &meanDiff, double &minDiff,
                                     double &maxDiff) {
  std::vector<double> refVec;
  chrom->GetVector(refVec);
  meanDiff = 0.0;
  minDiff = 0.0;
  maxDiff = 0.0;
  ChromElementPtr chrom2 = chrom->clone();
  ChromElementPtr chrom3 = chrom->clone();
  ChromElementPtr chrom4 = chrom->clone();
  for (int iTrial = 0; iTrial < nTrials; iTrial++) {
    chrom->Randomise();
    chrom2->Randomise();
    Crossover(chrom, chrom2, chrom3, chrom4);
    int i(0);
    double diff3 = chrom3->CompareVector(refVec, i);
    i = 0;
    double diff4 = chrom4->CompareVector(refVec, i);
    meanDiff += diff3;
    meanDiff += diff4;
    if (iTrial == 0) {
      minDiff = std::min(diff3, diff4);
      maxDiff = std::max(diff3, diff4);
    } else {
      minDiff = std::min(minDiff, diff3);
      minDiff = std::min(minDiff, diff4);
      maxDiff = std::max(maxDiff, diff3);
      maxDiff = std::max(maxDiff, diff4);
    }
  }
  if (nTrials > 0) {
    meanDiff /= (nTrials * 2);
  }
}
