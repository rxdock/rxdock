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

#include <iomanip>
using std::setw;

#include "NMSearch.h"
#include "RbtChrom.h"
#include "RbtSFRequest.h"
#include "RbtSimplexTransform.h"
#include "RbtWorkSpace.h"

// Static data member for class type
std::string RbtSimplexTransform::_CT("RbtSimplexTransform");
// Parameter names
std::string RbtSimplexTransform::_MAX_CALLS("MAX_CALLS");
std::string RbtSimplexTransform::_NCYCLES("NCYCLES");
std::string RbtSimplexTransform::_STOPPING_STEP_LENGTH("STOPPING_STEP_LENGTH");
std::string RbtSimplexTransform::_PARTITION_DIST("PARTITION_DIST");
std::string RbtSimplexTransform::_STEP_SIZE("STEP_SIZE");
std::string RbtSimplexTransform::_CONVERGENCE("CONVERGENCE");

////////////////////////////////////////
// Constructors/destructors
RbtSimplexTransform::RbtSimplexTransform(const std::string &strName)
    : RbtBaseBiMolTransform(_CT, strName) {
  AddParameter(_MAX_CALLS, 200);
  AddParameter(_NCYCLES, 5);
  AddParameter(_STOPPING_STEP_LENGTH, 10e-4);
  AddParameter(_PARTITION_DIST, 0.0);
  AddParameter(_STEP_SIZE, 0.1);
  AddParameter(_CONVERGENCE, 0.001);
#ifdef _DEBUG
  cout << _CT << " parameterised constructor" << endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtSimplexTransform::~RbtSimplexTransform() {
#ifdef _DEBUG
  cout << _CT << " destructor" << endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

////////////////////////////////////////
// Protected methods
///////////////////
void RbtSimplexTransform::SetupReceptor() {}

void RbtSimplexTransform::SetupLigand() {}

void RbtSimplexTransform::SetupSolvent() {}

void RbtSimplexTransform::SetupTransform() {
  // Construct the overall chromosome for the system
  m_chrom.SetNull();
  RbtWorkSpace *pWorkSpace = GetWorkSpace();
  if (pWorkSpace) {
    m_chrom = new RbtChrom(pWorkSpace->GetModels());
  }
}

////////////////////////////////////////
// Private methods
///////////////////
// Pure virtual in RbtBaseTransform
// Actually apply the transform
void RbtSimplexTransform::Execute() {
  // Get the current scoring function from the workspace
  RbtWorkSpace *pWorkSpace = GetWorkSpace();
  if (pWorkSpace == NULL) // Return if this transform is not registered
    return;
  RbtBaseSF *pSF = pWorkSpace->GetSF();
  if (pSF == NULL) // Return if workspace does not have a scoring function
    return;
  int iTrace = GetTrace();

  pWorkSpace->ClearPopulation();
  int maxcalls = GetParameter(_MAX_CALLS);
  int ncycles = GetParameter(_NCYCLES);
  double stopping = GetParameter(_STOPPING_STEP_LENGTH);
  double convergence = GetParameter(_CONVERGENCE);
  double stepSize = GetParameter(_STEP_SIZE);
  double partDist = GetParameter(_PARTITION_DIST);
  RbtRequestPtr spPartReq(new RbtSFPartitionRequest(partDist));
  RbtRequestPtr spClearPartReq(new RbtSFPartitionRequest(0.0));
  pSF->HandleRequest(spPartReq);

  m_chrom->SyncFromModel();
  // If we are minimising all degrees of freedom simultaneuously
  // we have to compile a vector of variable step sizes for the NMSearch
  RbtDoubleList sv;
  m_chrom->GetStepVector(sv);
  int nsv = sv.size();
  double *steps = new double[nsv];
  for (int i = 0; i < nsv; ++i) {
    steps[i] = sv[i] * stepSize;
  }

  // Set up the Simplex search object
  NMSearch *ssearch;
  NMSearch::SetMaxCalls(maxcalls);
  NMSearch::SetStoppingLength(stopping);
  int calls = 0;
  double initScore = pSF->Score(); // Current score
  double min = initScore;
  RbtDoubleList vc; // Vector representation of chromosome
  int N = ncycles;
  // Energy change between cycles - initialise so as not to terminate loop
  // immediately
  double delta = -convergence - 1.0;

  if (iTrace > 0) {
    cout.precision(3);
    cout.setf(ios_base::fixed, ios_base::floatfield);
    cout.setf(ios_base::right, ios_base::adjustfield);
    cout << endl
         << _CT << endl
         << setw(5) << "CYCLE" << setw(5) << "MODE" << setw(5) << "DOF"
         << setw(10) << "CALLS" << setw(10) << "SCORE" << setw(10) << "DELTA"
         << endl;
    cout << endl
         << setw(5) << "Init" << setw(5) << "-" << setw(5) << "-" << setw(10)
         << calls << setw(10) << initScore << setw(10) << "-" << endl
         << endl;
    if (iTrace > 1) {
      cout << *m_chrom << endl;
    }
  }

  for (int i = 0; (i < N) && (delta < -convergence); i++) {
    if (partDist > 0.0) {
      pSF->HandleRequest(spPartReq);
    }
    // Use a variable length simplex
    vc.clear();
    m_chrom->GetVector(vc);
    ssearch = new NMSearch(m_chrom, pSF);
    ssearch->InitVariableLengthRightSimplex(&vc, steps);
    if (iTrace > 0) {
      cout << setw(5) << i << setw(5) << "ALL" << setw(5) << vc.size();
    }
    // Do the simplex search and retrieve the minimum
    ssearch->ExploratoryMoves();
    double newmin = ssearch->GetMinVal();
    delta = newmin - min;
    calls += ssearch->GetFunctionCalls();
    m_chrom->SetVector(ssearch->GetMinPoint());
    if (iTrace > 0) {
      cout << setw(10) << calls << setw(10) << newmin << setw(10) << delta
           << endl;
      if (iTrace > 1) {
        cout << *m_chrom << endl;
      }
    }
    delete ssearch; // DM 27 Jun 2002 - garbage collection
    min = newmin;
  }
  m_chrom->SyncToModel();
  pSF->HandleRequest(spClearPartReq); // Clear any partitioning
  delete[] steps;
  if (iTrace > 0) {
    min = pSF->Score();
    delta = min - initScore;
    cout << endl
         << setw(5) << "Final" << setw(5) << "-" << setw(5) << "-" << setw(10)
         << calls << setw(10) << min << setw(10) << delta << endl;
  }
}
