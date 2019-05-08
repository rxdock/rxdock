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

#include "RbtGATransform.h"
#include "RbtPopulation.h"
#include "RbtSFRequest.h"
#include "RbtWorkSpace.h"
#include <iomanip>
using std::setw;

std::string RbtGATransform::_CT("RbtGATransform");
std::string RbtGATransform::_NEW_FRACTION("NEW_FRACTION");
std::string RbtGATransform::_PCROSSOVER("PCROSSOVER");
std::string RbtGATransform::_XOVERMUT("XOVERMUT");
std::string RbtGATransform::_CMUTATE("CMUTATE");
std::string RbtGATransform::_STEP_SIZE("STEP_SIZE");
std::string RbtGATransform::_EQUALITY_THRESHOLD("EQUALITY_THRESHOLD");
std::string RbtGATransform::_NCYCLES("NCYCLES");
std::string RbtGATransform::_NCONVERGENCE("NCONVERGENCE");
std::string RbtGATransform::_HISTORY_FREQ("HISTORY_FREQ");

RbtGATransform::RbtGATransform(const std::string &strName)
    : RbtBaseBiMolTransform(_CT, strName), m_rand(Rbt::GetRbtRand()) {
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

RbtGATransform::~RbtGATransform() { _RBTOBJECTCOUNTER_DESTR_(_CT); }

void RbtGATransform::SetupReceptor() {}

void RbtGATransform::SetupLigand() {}

void RbtGATransform::SetupTransform() {}

void RbtGATransform::Execute() {
  RbtWorkSpace *pWorkSpace = GetWorkSpace();
  if (pWorkSpace == NULL) {
    return;
  }
  RbtBaseSF *pSF = pWorkSpace->GetSF();
  if (pSF == NULL) {
    return;
  }
  RbtPopulationPtr pop = pWorkSpace->GetPopulation();
  if (pop.Null() || (pop->GetMaxSize() < 1)) {
    return;
  }
  // Remove any partitioning from the scoring function
  // Not appropriate for a GA
  pSF->HandleRequest(new RbtSFPartitionRequest(0.0));
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

  int popsize = pop->GetMaxSize();
  int nrepl = newFraction * popsize;
  bool bHistory = nHisFreq > 0;
  int iTrace = GetTrace();

  double bestScore = pop->Best()->GetScore();
  // Number of consecutive cycles with no improvement in best score
  int iConvergence = 0;

  if (iTrace > 0) {
    cout.precision(3);
    cout.setf(ios_base::fixed, ios_base::floatfield);
    cout.setf(ios_base::right, ios_base::adjustfield);
    cout << endl
         << setw(5) << "CYCLE" << setw(5) << "CONV" << setw(10) << "BEST"
         << setw(10) << "MEAN" << setw(10) << "VAR" << endl;

    cout << endl
         << setw(5) << "Init" << setw(5) << "-" << setw(10) << bestScore
         << setw(10) << pop->GetScoreMean() << setw(10)
         << pop->GetScoreVariance() << endl;
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
      cout << setw(5) << iCycle << setw(5) << iConvergence << setw(10) << score
           << setw(10) << pop->GetScoreMean() << setw(10)
           << pop->GetScoreVariance() << endl;
    }
  }
  pop->Best()->GetChrom()->SyncToModel();
  int ri = GetReceptor()->GetCurrentCoords();
  GetLigand()->SetDataValue("RI", ri);
}
