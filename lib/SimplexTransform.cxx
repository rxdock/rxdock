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

#ifdef __PGI
#define EIGEN_DONT_VECTORIZE
#endif
#include <Eigen/Core>
#include <iomanip>

#include "rxdock/Chrom.h"
#include "rxdock/NMCriteria.h"
#include "rxdock/NMSimplex.h"
#include "rxdock/SFRequest.h"
#include "rxdock/SimplexCostFunction.h"
#include "rxdock/SimplexTransform.h"
#include "rxdock/WorkSpace.h"

#include <fmt/ostream.h>
#include <loguru.hpp>

using namespace rxdock;

// Static data member for class type
std::string SimplexTransform::_CT("SimplexTransform");
// Parameter names
std::string SimplexTransform::_MAX_CALLS("MAX_CALLS");
std::string SimplexTransform::_NCYCLES("NCYCLES");
std::string SimplexTransform::_STOPPING_STEP_LENGTH("STOPPING_STEP_LENGTH");
std::string SimplexTransform::_PARTITION_DIST("PARTITION_DIST");
std::string SimplexTransform::_STEP_SIZE("STEP_SIZE");
std::string SimplexTransform::_CONVERGENCE("CONVERGENCE");

std::string &SimplexTransform::GetMaxCalls() { return _MAX_CALLS; }

std::string &SimplexTransform::GetNCycles() { return _NCYCLES; }

std::string &SimplexTransform::GetStepSize() { return _STEP_SIZE; }

////////////////////////////////////////
// Constructors/destructors
SimplexTransform::SimplexTransform(const std::string &strName)
    : BaseBiMolTransform(_CT, strName) {
  LOG_F(2, "SimplexTransform parameterised constructor");
  AddParameter(_MAX_CALLS, 200);
  AddParameter(_NCYCLES, 5);
  AddParameter(_STOPPING_STEP_LENGTH, 10e-4);
  AddParameter(_PARTITION_DIST, 0.0);
  AddParameter(_STEP_SIZE, 0.1);
  AddParameter(_CONVERGENCE, 0.001);
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

SimplexTransform::~SimplexTransform() {
  LOG_F(2, "SimplexTransform destructor");
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

////////////////////////////////////////
// Protected methods
///////////////////
void SimplexTransform::SetupReceptor() {}

void SimplexTransform::SetupLigand() {}

void SimplexTransform::SetupSolvent() {}

void SimplexTransform::SetupTransform() {
  // Construct the overall chromosome for the system
  m_chrom.SetNull();
  WorkSpace *pWorkSpace = GetWorkSpace();
  if (pWorkSpace) {
    m_chrom = new Chrom(pWorkSpace->GetModels());
  }
}

////////////////////////////////////////
// Private methods
///////////////////
// Pure virtual in BaseTransform
// Actually apply the transform
void SimplexTransform::Execute() {
  LOG_F(2, "SimplexTransform::Execute");
  // Get the current scoring function from the workspace
  WorkSpace *pWorkSpace = GetWorkSpace();
  if (pWorkSpace == nullptr) // Return if this transform is not registered
    return;
  BaseSF *pSF = pWorkSpace->GetSF();
  if (pSF == nullptr) // Return if workspace does not have a scoring function
    return;

  pWorkSpace->ClearPopulation();
  int maxcalls = GetParameter(_MAX_CALLS);
  int ncycles = GetParameter(_NCYCLES);
  double stopping = GetParameter(_STOPPING_STEP_LENGTH);
  double convergence = GetParameter(_CONVERGENCE);
  double stepSize = GetParameter(_STEP_SIZE);
  double partDist = GetParameter(_PARTITION_DIST);
  RequestPtr spPartReq(new SFPartitionRequest(partDist));
  RequestPtr spClearPartReq(new SFPartitionRequest(0.0));
  pSF->HandleRequest(spPartReq);

  m_chrom->SyncFromModel();
  // If we are minimising all degrees of freedom simultaneuously
  // we have to compile a vector of variable step sizes for the NMSearch
  std::vector<double> sv;
  m_chrom->GetStepVector(sv);
  int nsv = sv.size();
  double *steps = new double[nsv];
  for (int i = 0; i < nsv; ++i) {
    steps[i] = sv[i] * stepSize;
  }

  SimplexCostFunction costFunction(pSF, m_chrom);

  // Builder to generate the optimizer with a composite stoping criterion
  auto optimizer = neldermead::CreateSimplex(
      costFunction, neldermead::CreateAndCriteria(
                        neldermead::IterationCriterion(maxcalls),
                        neldermead::RelativeValueCriterion<double>(stopping)));

  int calls = 0;
  double initScore = pSF->Score(); // Current score
  double min = initScore;
  std::vector<double> vc; // Vector representation of chromosome
  // Energy change between cycles - initialise so as not to terminate loop
  // immediately
  double delta = -convergence - 1.0;

  LOG_F(INFO, "CYCLE MODE  DOF     CALLS     SCORE     DELTA");
  LOG_F(INFO, " Init    -    -{:10d}{:10f}         -", calls, initScore);
  LOG_F(1, "{}", *m_chrom);

  for (int i = 0; (i < ncycles) && (delta < -convergence); i++) {
    if (partDist > 0.0) {
      pSF->HandleRequest(spPartReq);
    }
    // Use a variable length simplex
    vc.clear();
    m_chrom->GetVector(vc);
    SimplexCostFunction::ParameterType start_point =
        Eigen::Map<SimplexCostFunction::ParameterType, Eigen::Unaligned>(
            vc.data(), vc.size());

    optimizer.SetStartPoint(start_point); // Starting parameters
    optimizer.SetDelta(delta);            // Simplex size
    optimizer.Optimize(costFunction);     // Optimization start

    double newmin = optimizer.GetBestValue();
    delta = newmin - min;
    calls += costFunction.nCalls;
    SimplexCostFunction::ParameterType best = optimizer.GetBestParameters();
    std::vector<double> best_vec(best.data(), best.data() + best.size());
    m_chrom->SetVector(best_vec);

    LOG_F(INFO, "{:5d}  ALL{:5d}{:10d}{:10f}{:10f}", i, vc.size(), calls,
          newmin, delta);
    LOG_F(1, "{}", *m_chrom);

    min = newmin;
  }
  m_chrom->SyncToModel();
  pSF->HandleRequest(spClearPartReq); // Clear any partitioning
  delete[] steps;
  LOG_F(INFO, "Final    -    -{:10d}{:10f}{:10f}", calls, pSF->Score(),
        pSF->Score() - initScore);
}
