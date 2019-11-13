/***********************************************************************
 * The rDock program was developed from 1998 - 2006 by the software team
 * at RiboTargets (subsequently Vernalis (R&D) Ltd).
 * In 2006, the software was licensed to the University of York for
 * maintenance and distribution.
 * In 2012, Vernalis and the University of York agreed to release the
 * program as Open Source software.
 * This version is licensed under GNU-LGPL version 3.0 with support from
 * the University of Barcelona.
 * http://rdock.sourceforge.net/
 ***********************************************************************/

// Simple, non-adaptive simulated annealing protocol
#include <iomanip>

#include "RbtBaseSF.h"
#include "RbtChrom.h"
#include "RbtSFRequest.h"
#include "RbtSimAnnTransform.h"
#include "RbtWorkSpace.h"

// Simple class to keep track of Monte Carlo sampling statistics
RbtMCStats::RbtMCStats() { Init(0.0); }

void RbtMCStats::Init(double score) {
  _min = _max = _initial = _final = score;
  InitBlock(score);
}

void RbtMCStats::InitBlock(double score) {
  _blockMin = _blockMax = _blockInitial = _blockFinal = score;
  _total = 0.0;
  _total2 = 0.0;
  _steps = 0;
  _accepted = 0;
}

void RbtMCStats::Accumulate(double score, bool bAccepted) {
  _steps++;
  if (bAccepted)
    _accepted++;
  _total += score;
  _total2 += score * score;
  _blockMin = std::min(_blockMin, score);
  _blockMax = std::max(_blockMax, score);
  _blockFinal = _final = score;
  _min = std::min(_min, score);
  _max = std::max(_max, score);
}

double RbtMCStats::Mean() const { return _total / _steps; }
double RbtMCStats::Variance() const {
  return _total2 / _steps - std::pow(Mean(), 2);
}
double RbtMCStats::AccRate() const { return float(_accepted) / float(_steps); }

// Static data member for class type
std::string RbtSimAnnTransform::_CT("RbtSimAnnTransform");
// Parameter names
std::string RbtSimAnnTransform::_START_T("START_T");
std::string RbtSimAnnTransform::_FINAL_T("FINAL_T");
std::string RbtSimAnnTransform::_BLOCK_LENGTH("BLOCK_LENGTH");
std::string RbtSimAnnTransform::_SCALE_CHROM_LENGTH("SCALE_CHROM_LENGTH");
std::string RbtSimAnnTransform::_NUM_BLOCKS("NUM_BLOCKS");
std::string RbtSimAnnTransform::_STEP_SIZE("STEP_SIZE");
std::string RbtSimAnnTransform::_MIN_ACC_RATE("MIN_ACC_RATE");
std::string RbtSimAnnTransform::_PARTITION_DIST("PARTITION_DIST");
std::string RbtSimAnnTransform::_PARTITION_FREQ("PARTITION_FREQ");
std::string RbtSimAnnTransform::_HISTORY_FREQ("HISTORY_FREQ");

std::string &RbtSimAnnTransform::GetStartT() { return _START_T; }

std::string &RbtSimAnnTransform::GetFinalT() { return _FINAL_T; }

std::string &RbtSimAnnTransform::GetBlockLength() { return _BLOCK_LENGTH; }

std::string &RbtSimAnnTransform::GetNumBlocks() { return _NUM_BLOCKS; }

std::string &RbtSimAnnTransform::GetStepSize() { return _STEP_SIZE; }

std::string &RbtSimAnnTransform::GetPartitionDist() { return _PARTITION_DIST; }

std::string &RbtSimAnnTransform::GetPartitionFreq() { return _PARTITION_FREQ; }

////////////////////////////////////////
// Constructors/destructors
RbtSimAnnTransform::RbtSimAnnTransform(const std::string &strName)
    : RbtBaseBiMolTransform(_CT, strName), m_rand(Rbt::GetRbtRand()) {
  // Add parameters
  AddParameter(_START_T, 1000.0);
  AddParameter(_FINAL_T, 300.0);
  AddParameter(_BLOCK_LENGTH, 50);
  AddParameter(_SCALE_CHROM_LENGTH, true);
  AddParameter(_NUM_BLOCKS, 25);
  AddParameter(_STEP_SIZE, 1.0);
  AddParameter(_MIN_ACC_RATE, 0.25);
  AddParameter(_PARTITION_DIST, 0.0);
  AddParameter(_PARTITION_FREQ, 0);
  AddParameter(_HISTORY_FREQ, 0);
  m_spStats = RbtMCStatsPtr(new RbtMCStats());
#ifdef _DEBUG
  std::cout << _CT << " parameterised constructor" << std::endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtSimAnnTransform::~RbtSimAnnTransform() {
#ifdef _DEBUG
  std::cout << _CT << " destructor" << std::endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

////////////////////////////////////////
// Protected methods
///////////////////
void RbtSimAnnTransform::SetupReceptor() {}

void RbtSimAnnTransform::SetupLigand() {}

void RbtSimAnnTransform::SetupTransform() {
  // Construct the overall chromosome for the system
  m_chrom.SetNull();
  m_lastGoodVector.clear();
  m_minVector.clear();
  RbtWorkSpace *pWorkSpace = GetWorkSpace();
  if (pWorkSpace) {
    m_chrom = new RbtChrom(pWorkSpace->GetModels());
    int chromLength = m_chrom->GetLength();
    m_lastGoodVector.reserve(chromLength);
    m_minVector.reserve(chromLength);
  }
}

////////////////////////////////////////
// Private methods
///////////////////
// Pure virtual in RbtBaseTransform
// Actually apply the transform
void RbtSimAnnTransform::Execute() {
  // Get the current scoring function from the workspace
  RbtWorkSpace *pWorkSpace = GetWorkSpace();
  if (pWorkSpace == nullptr) // Return if this transform is not registered
    return;
  RbtBaseSF *pSF = pWorkSpace->GetSF();
  if (pSF == nullptr) // Return if workspace does not have a scoring function
    return;
  int iTrace = GetTrace();

  pWorkSpace->ClearPopulation();
  // Get the cooling schedule params
  double t = GetParameter(_START_T);
  double tFinal = GetParameter(_FINAL_T);
  int nBlocks = GetParameter(_NUM_BLOCKS);
  int blockLen = GetParameter(_BLOCK_LENGTH);
  bool bScale = GetParameter(_SCALE_CHROM_LENGTH);
  double stepSize = GetParameter(_STEP_SIZE);
  double minAccRate = GetParameter(_MIN_ACC_RATE);

  if (bScale) {
    int chromLength = m_chrom->GetLength();
    blockLen *= chromLength;
  }
  if (iTrace > 0) {
    std::cout << _CT << ": Block length = " << blockLen << std::endl;
  }

  // Cooling factor (check for nBlocks=1)
  double tFac = (nBlocks > 1) ? std::pow(tFinal / t, 1.0 / (nBlocks - 1)) : 1.0;

  // DM 15 Feb 1999 - don't initialise the Monte Carlo stats each block
  // if we are doing a constant temperature run
  bool bInitBlock = (t != tFinal);

  // Send the partitioning request separately based on the current partition
  // distance If partDist is zero, the partitioning automatically gets removed
  double partDist = GetParameter(_PARTITION_DIST);
  m_spPartReq = RbtRequestPtr(new RbtSFPartitionRequest(partDist));
  pSF->HandleRequest(m_spPartReq);

  // Update the chromosome to match the current model coords
  // and keep track of the vector corresponding to minimum score
  m_chrom->SyncFromModel();
  m_minVector.clear();
  m_chrom->GetVector(m_minVector);
  double score = pSF->Score();
  m_spStats->Init(score);

  std::cout.precision(3);
  std::cout.setf(std::ios_base::fixed, std::ios_base::floatfield);
  std::cout.setf(std::ios_base::right, std::ios_base::adjustfield);

  if (iTrace > 0) {
    std::cout << _CT << ": Initial score = " << score << std::endl;
    std::cout << std::endl
              << std::endl
              << std::setw(5) << "BLOCK" << std::setw(10) << "TEMP"
              << std::setw(10) << "ACC.RATE" << std::setw(10) << "STEP"
              << std::setw(10) << "INITIAL" << std::setw(10) << "FINAL"
              << std::setw(10) << "MEAN" << std::setw(10) << "S.DEV"
              << std::setw(10) << "MIN" << std::setw(10) << "MAX" << std::endl;
  }

  for (int iBlock = 1; iBlock <= nBlocks; iBlock++, t *= tFac) {
    if (bInitBlock) {
      m_spStats->InitBlock(pSF->Score());
    }
    MC(t, blockLen, stepSize);
    if (iTrace > 0) {
      std::cout << std::setw(5) << iBlock << std::setw(10) << t << std::setw(10)
                << m_spStats->AccRate() << std::setw(10) << stepSize
                << std::setw(10) << m_spStats->_blockInitial << std::setw(10)
                << m_spStats->_blockFinal << std::setw(10) << m_spStats->Mean()
                << std::setw(10) << std::sqrt(m_spStats->Variance())
                << std::setw(10) << m_spStats->_blockMin << std::setw(10)
                << m_spStats->_blockMax << std::endl;
    }

    // Halve the maximum step sizes for all enabled modes
    // if the acceptance rate is less than the threshold
    if (m_spStats->AccRate() < minAccRate) {
      stepSize *= 0.5;
      // Reinitialise the stats (only need to do it here if bInitBlock is false
      // otherwise it will be done at the beginning of the next block)
      if (!bInitBlock) {
        m_spStats->InitBlock(pSF->Score());
      }
    }
  }
  // Update the model coords with the minimum score chromosome
  m_chrom->SetVector(m_minVector);
  m_chrom->SyncToModel();
  RbtRequestPtr spClearPartReq(new RbtSFPartitionRequest(0.0));
  pSF->HandleRequest(spClearPartReq); // Clear any partitioning
  if (iTrace > 0) {
    std::cout << std::endl
              << _CT << ": Final score = " << pSF->Score() << std::endl;
  }
}

void RbtSimAnnTransform::MC(double t, int blockLen, double stepSize) {
  RbtBaseSF *pSF = GetWorkSpace()->GetSF();
  int iTrace = GetTrace();
  double score = pSF->Score();

  int nHisFreq = GetParameter(_HISTORY_FREQ);
  bool bHistory = nHisFreq > 0;
  int nPartFreq = GetParameter(_PARTITION_FREQ);
  bool bPartition = (nPartFreq > 0);

  // Keep a record of the last good chromosome vector, for fast revert following
  // a failed Metropolic test
  m_lastGoodVector.clear();
  m_chrom->GetVector(m_lastGoodVector);
  // Main loop over number of MC steps
  for (int iStep = 1; iStep <= blockLen; iStep++) {
    m_chrom->Mutate(stepSize);
    m_chrom->SyncToModel();
    double newScore = pSF->Score();
    double delta = newScore - score;
    bool bMetrop = ((delta < 0.0) || (std::exp(-1000.0 * delta / (8.314 * t)) >
                                      m_rand.GetRandom01()));
    // PASSED
    if (bMetrop) {
      score = newScore;
      // Update the last good vector
      m_lastGoodVector.clear();
      m_chrom->GetVector(m_lastGoodVector);
      // Update the minimum score vector
      if (score < m_spStats->_min) {
        m_minVector.clear();
        m_chrom->GetVector(m_minVector);
      }
    }
    // FAILED
    else {
      // revert to old chromosome
      // No need to SyncToModel as this will be done after the next mutation
      m_chrom->SetVector(m_lastGoodVector);
    }
    // Gather the statistics
    m_spStats->Accumulate(score, bMetrop);
    // Render to the history file if appropriate (true = with component scores)
    if (bHistory && (iStep % nHisFreq) == 0) {
      GetWorkSpace()->SaveHistory(true);
    }

    // Update the interaction lists if appropriate
    // We update the lists every nth accepted trial (as rejected trials don't
    // change the coords)
    if (bPartition && ((m_spStats->_accepted % nPartFreq) == 0)) {
      pSF->HandleRequest(m_spPartReq);
      double oldScore = score;
      score = pSF->Score();
      if (std::fabs(score - oldScore) > 0.001 && (iTrace > 1)) {
        std::cout
            << "** WARNING - Interaction lists updated, change in score = "
            << score - oldScore << std::endl;
      }
    }
  }
  // Ensure the model is synchronised with the chromosome on exit
  // This would not be the case if the last Metropolis test were failed
  m_chrom->SyncToModel();
}
