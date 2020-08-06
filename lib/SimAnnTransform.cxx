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

#include "rxdock/BaseSF.h"
#include "rxdock/Chrom.h"
#include "rxdock/SFRequest.h"
#include "rxdock/SimAnnTransform.h"
#include "rxdock/WorkSpace.h"

#include <loguru.hpp>

using namespace rxdock;

// Simple class to keep track of Monte Carlo sampling statistics
MCStats::MCStats() { Init(0.0); }

void MCStats::Init(double score) {
  _min = _max = _initial = _final = score;
  InitBlock(score);
}

void MCStats::InitBlock(double score) {
  _blockMin = _blockMax = _blockInitial = _blockFinal = score;
  _total = 0.0;
  _total2 = 0.0;
  _steps = 0;
  _accepted = 0;
}

void MCStats::Accumulate(double score, bool bAccepted) {
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

double MCStats::Mean() const { return _total / _steps; }
double MCStats::Variance() const {
  return _total2 / _steps - std::pow(Mean(), 2);
}
double MCStats::AccRate() const { return float(_accepted) / float(_steps); }

// Static data member for class type
std::string SimAnnTransform::_CT("SimAnnTransform");
// Parameter names
std::string SimAnnTransform::_START_T("START_T");
std::string SimAnnTransform::_FINAL_T("FINAL_T");
std::string SimAnnTransform::_BLOCK_LENGTH("BLOCK_LENGTH");
std::string SimAnnTransform::_SCALE_CHROM_LENGTH("SCALE_CHROM_LENGTH");
std::string SimAnnTransform::_NUM_BLOCKS("NUM_BLOCKS");
std::string SimAnnTransform::_STEP_SIZE("STEP_SIZE");
std::string SimAnnTransform::_MIN_ACC_RATE("MIN_ACC_RATE");
std::string SimAnnTransform::_PARTITION_DIST("PARTITION_DIST");
std::string SimAnnTransform::_PARTITION_FREQ("PARTITION_FREQ");
std::string SimAnnTransform::_HISTORY_FREQ("HISTORY_FREQ");

std::string &SimAnnTransform::GetStartT() { return _START_T; }

std::string &SimAnnTransform::GetFinalT() { return _FINAL_T; }

std::string &SimAnnTransform::GetBlockLength() { return _BLOCK_LENGTH; }

std::string &SimAnnTransform::GetNumBlocks() { return _NUM_BLOCKS; }

std::string &SimAnnTransform::GetStepSize() { return _STEP_SIZE; }

std::string &SimAnnTransform::GetPartitionDist() { return _PARTITION_DIST; }

std::string &SimAnnTransform::GetPartitionFreq() { return _PARTITION_FREQ; }

////////////////////////////////////////
// Constructors/destructors
SimAnnTransform::SimAnnTransform(const std::string &strName)
    : BaseBiMolTransform(_CT, strName), m_rand(GetRandInstance()) {
  LOG_F(2, "SimAnnTransform parameterised constructor");
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
  m_spStats = MCStatsPtr(new MCStats());
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

SimAnnTransform::~SimAnnTransform() {
  LOG_F(2, "SimAnnTransform destructor");
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

////////////////////////////////////////
// Protected methods
///////////////////
void SimAnnTransform::SetupReceptor() {}

void SimAnnTransform::SetupLigand() {}

void SimAnnTransform::SetupTransform() {
  // Construct the overall chromosome for the system
  m_chrom.SetNull();
  m_lastGoodVector.clear();
  m_minVector.clear();
  WorkSpace *pWorkSpace = GetWorkSpace();
  if (pWorkSpace) {
    m_chrom = new Chrom(pWorkSpace->GetModels());
    int chromLength = m_chrom->GetLength();
    m_lastGoodVector.reserve(chromLength);
    m_minVector.reserve(chromLength);
  }
}

////////////////////////////////////////
// Private methods
///////////////////
// Pure virtual in BaseTransform
// Actually apply the transform
void SimAnnTransform::Execute() {
  // Get the current scoring function from the workspace
  WorkSpace *pWorkSpace = GetWorkSpace();
  if (pWorkSpace == nullptr) // Return if this transform is not registered
    return;
  BaseSF *pSF = pWorkSpace->GetSF();
  if (pSF == nullptr) // Return if workspace does not have a scoring function
    return;

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
  LOG_F(INFO, "SimAnnTransform::Execute: Block length = {}", blockLen);

  // Cooling factor (check for nBlocks=1)
  double tFac = (nBlocks > 1) ? std::pow(tFinal / t, 1.0 / (nBlocks - 1)) : 1.0;

  // DM 15 Feb 1999 - don't initialise the Monte Carlo stats each block
  // if we are doing a constant temperature run
  bool bInitBlock = (t != tFinal);

  // Send the partitioning request separately based on the current partition
  // distance If partDist is zero, the partitioning automatically gets removed
  double partDist = GetParameter(_PARTITION_DIST);
  m_spPartReq = RequestPtr(new SFPartitionRequest(partDist));
  pSF->HandleRequest(m_spPartReq);

  // Update the chromosome to match the current model coords
  // and keep track of the vector corresponding to minimum score
  m_chrom->SyncFromModel();
  m_minVector.clear();
  m_chrom->GetVector(m_minVector);
  double score = pSF->Score();
  m_spStats->Init(score);

  LOG_F(INFO, "SimAnnTransform::Execute: Initial score = {}", score);
  LOG_F(INFO, "BLOCK      TEMP  ACC.RATE      STEP   INITIAL     FINAL      "
              "MEAN     S.DEV       MIN       MAX");
  for (int iBlock = 1; iBlock <= nBlocks; iBlock++, t *= tFac) {
    if (bInitBlock) {
      m_spStats->InitBlock(pSF->Score());
    }
    MC(t, blockLen, stepSize);
    LOG_F(INFO,
          "{:5d}{:10.3f}{:10.3f}{:10.3f}{:10.3f}{:10.3f}{:10.3f}{:10.3f}"
          "{:10.3f}{:10.3f}",
          iBlock, t, m_spStats->AccRate(), stepSize, m_spStats->_blockInitial,
          m_spStats->_blockFinal, m_spStats->Mean(),
          std::sqrt(m_spStats->Variance()), m_spStats->_blockMin,
          m_spStats->_blockMax);

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
  RequestPtr spClearPartReq(new SFPartitionRequest(0.0));
  pSF->HandleRequest(spClearPartReq); // Clear any partitioning
  LOG_F(INFO, "SimAnnTransform: Final score = {}", pSF->Score());
}

void SimAnnTransform::MC(double t, int blockLen, double stepSize) {
  BaseSF *pSF = GetWorkSpace()->GetSF();
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
      if (std::fabs(score - oldScore) > 0.001) {
        LOG_F(WARNING, "Interaction lists updated, change in score = {}",
              score - oldScore);
      }
    }
  }
  // Ensure the model is synchronised with the chromosome on exit
  // This would not be the case if the last Metropolis test were failed
  m_chrom->SyncToModel();
}
