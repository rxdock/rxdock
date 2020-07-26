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

#include "GATransform.h"
#include "Population.h"
#include "SFRequest.h"
#include "WorkSpace.h"
#include <iomanip>

using namespace rxdock;

std::string GATransform::_CT("GATransform");
std::string GATransform::_NEW_FRACTION("NEW_FRACTION");
std::string GATransform::_PCROSSOVER("PCROSSOVER");
std::string GATransform::_XOVERMUT("XOVERMUT");
std::string GATransform::_CMUTATE("CMUTATE");
std::string GATransform::_STEP_SIZE("STEP_SIZE");
std::string GATransform::_EQUALITY_THRESHOLD("EQUALITY_THRESHOLD");
std::string GATransform::_NCYCLES("NCYCLES");
std::string GATransform::_NCONVERGENCE("NCONVERGENCE");
std::string GATransform::_HISTORY_FREQ("HISTORY_FREQ");

GATransform::GATransform(const std::string &strName)
    : BaseBiMolTransform(_CT, strName), m_rand(GetRandInstance()) {
  AddParameter(_NEW_FRACTION, 0.5);
  AddParameter(_PCROSSOVER, 0.4);
  AddParameter(_XOVERMUT, true);
  AddParameter(_CMUTATE, false);
  AddParameter(_STEP_SIZE, 1.0);
  AddParameter(_EQUALITY_THRESHOLD, 0.1);
  AddParameter(_NCYCLES, 100);
  AddParameter(_NCONVERGENCE, 6);
  AddParameter(_HISTORY_FREQ, 0);
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

GATransform::~GATransform() { _RBTOBJECTCOUNTER_DESTR_(_CT); }

void GATransform::SetupReceptor() {}

void GATransform::SetupLigand() {}

void GATransform::SetupTransform() {}

void GATransform::Execute() {
  WorkSpace *pWorkSpace = GetWorkSpace();
  if (pWorkSpace == nullptr) {
    return;
  }
  BaseSF *pSF = pWorkSpace->GetSF();
  if (pSF == nullptr) {
    return;
  }
  PopulationPtr pop = pWorkSpace->GetPopulation();
  if (pop.Null() || (pop->GetMaxSize() < 1)) {
    return;
  }
  // Remove any partitioning from the scoring function
  // Not appropriate for a GA
  pSF->HandleRequest(new SFPartitionRequest(0.0));
  // This forces the population to rescore all the individuals in case
  // the scoring function has changed
  pop->SetSF(pSF);

  double newFraction = GetParameter(_NEW_FRACTION);
  double pcross = GetParameter(_PCROSSOVER);
  bool xovermut = GetParameter(_XOVERMUT);
  bool cmutate = GetParameter(_CMUTATE);
  double relStepSize = GetParameter(_STEP_SIZE);
  double equalityThreshold = GetParameter(_EQUALITY_THRESHOLD);
  int nCycles = GetParameter(_NCYCLES);
  int nConvergence = GetParameter(_NCONVERGENCE);
  int nHisFreq = GetParameter(_HISTORY_FREQ);

  double popsize = static_cast<double>(pop->GetMaxSize());
  int nrepl = static_cast<int>(newFraction * popsize);
  bool bHistory = nHisFreq > 0;
  int iTrace = GetTrace();

  double bestScore = pop->Best()->GetScore();
  // Number of consecutive cycles with no improvement in best score
  int iConvergence = 0;

  if (iTrace > 0) {
    std::cout.precision(3);
    std::cout.setf(std::ios_base::fixed, std::ios_base::floatfield);
    std::cout.setf(std::ios_base::right, std::ios_base::adjustfield);
    std::cout << std::endl
              << std::setw(5) << "CYCLE" << std::setw(5) << "CONV"
              << std::setw(10) << "BEST" << std::setw(10) << "MEAN"
              << std::setw(10) << "VAR" << std::endl;

    std::cout << std::endl
              << std::setw(5) << "Init" << std::setw(5) << "-" << std::setw(10)
              << bestScore << std::setw(10) << pop->GetScoreMean()
              << std::setw(10) << pop->GetScoreVariance() << std::endl;
  }

  for (int iCycle = 0; (iCycle < nCycles) && (iConvergence < nConvergence);
       ++iCycle) {
    if (bHistory && ((iCycle % nHisFreq) == 0)) {
      pop->Best()->GetChrom()->SyncToModel();
      pWorkSpace->SaveHistory(true);
    }
    pop->GAstep(nrepl, relStepSize, equalityThreshold, pcross, xovermut,
                cmutate);
    double score = pop->Best()->GetScore();
    if (score > bestScore) {
      bestScore = score;
      iConvergence = 0;
    } else {
      iConvergence++;
    }
    if (iTrace > 0) {
      std::cout << std::setw(5) << iCycle << std::setw(5) << iConvergence
                << std::setw(10) << score << std::setw(10)
                << pop->GetScoreMean() << std::setw(10)
                << pop->GetScoreVariance() << std::endl;
    }
  }
  pop->Best()->GetChrom()->SyncToModel();
  int ri = GetReceptor()->GetCurrentCoords();
  GetLigand()->SetDataValue(GetMetaDataPrefix() + "ri", ri);
}
