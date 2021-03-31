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

#include <sstream>

#include "rxdock/Filter.h"
#include "rxdock/Parser.h"
#include "rxdock/StringTokenIter.h"
#include "rxdock/WorkSpace.h"

#include <loguru.hpp>

using namespace rxdock;

const int STOP = 0;
const int CONT = -1;
const int NEXT = 1;

// Static data members
const std::string Filter::_CT = "Filter";
////////////////////////////////////////
// Constructors/destructors

// filter is used to determining if strfilter actually contains the
// filter (filter = true) or if strfilter is the name of a file
// that contains the filter (filter = false) This is the most
// common, so filter's value by default is false
Filter::Filter(std::string strfilter, bool filter) : BaseObject(_CT, "Filter") {
  LOG_F(2, "Filter default constructor");
  //  std::string filterfilen = GetParameter("_FILTER_FILE");
  SmartPtr<std::istream> filterfile;
  if (filter) // filterfilen is the filter
    filterfile = new std::istringstream(strfilter);
  else // assume the strfilter is the name of the file
       // where the filter is
    filterfile = new std::ifstream(strfilter.c_str(), std::ios_base::in);
  (*filterfile) >> nTermFilters;
  filteridx = 0;
  contextp = ContextPtr(new StringContext(filterfile));
  Parser p;
  for (int i = 0; i < nTermFilters; i++) {
    LOG_F(INFO, "------------- Terminate filter {} -------------", i);
    std::string s;
    getline(*filterfile, s, ',');
    SmartPtr<std::istream> istrp(new std::istringstream(s));
    TokenIterPtr ti(new StringTokenIter(istrp, contextp));
    FilterExpressionPtr filterExpr = p.Parse(ti, contextp);
    PrettyPrintVisitor visitor1(contextp);
    filterExpr->Accept(visitor1);
    terminationFilters.push_back(filterExpr);
  }
  (*filterfile) >> nWriteFilters;
  for (int i = 0; i < nWriteFilters; i++) {
    LOG_F(INFO, "------------- Write filter {} -----------------", i);
    std::string s;
    getline(*filterfile, s, ',');
    SmartPtr<std::istream> istrp(new std::istringstream(s));
    TokenIterPtr ti(new StringTokenIter(istrp, contextp));
    FilterExpressionPtr filterExpr = p.Parse(ti, contextp);
    PrettyPrintVisitor visitor1(contextp);
    filterExpr->Accept(visitor1);
    writtingFilter.push_back(filterExpr);
  }
  maxnruns = 1000;
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

Filter::~Filter() {
  LOG_F(2, "Filter destructor");
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

////////////////////////////////////////
// Public methods
////////////////

ModelPtr Filter::GetReceptor() const { return m_spReceptor; }
ModelPtr Filter::GetLigand() const { return m_spLigand; }

// Override Observer pure virtual
// Notify observer that subject has changed
void Filter::Update(Subject *theChangedSubject) {
  WorkSpace *pWorkSpace = GetWorkSpace();
  if (theChangedSubject == pWorkSpace) {
    // Check if receptor has been updated (model #0)
    if (pWorkSpace->GetNumModels() >= 1) {
      ModelPtr spReceptor = GetWorkSpace()->GetModel(0);
      if (spReceptor != m_spReceptor) {
        LOG_F(1, "BaseInterSF::Update: Receptor has been updated");
        m_spReceptor = spReceptor;
        SetupReceptor();
      }
    }
    // Check if ligand has been updated (model #1)
    if (pWorkSpace->GetNumModels() >= 2) {
      ModelPtr spLigand = GetWorkSpace()->GetModel(1);
      if (spLigand != m_spLigand) {
        LOG_F(1, "BaseInterSF::Update: Ligand has been updated");
        m_spLigand = spLigand;
        SetupLigand();
      }
    }
  }
}

// Called by Update when receptor is changed
void Filter::SetupReceptor() {
  // read filter file
  ((StringContextPtr)contextp)
      ->UpdateSite(m_spReceptor, GetWorkSpace()->GetDockingSite());
}

// Called by Update when ligand is changed
void Filter::SetupLigand() {
  ((StringContextPtr)contextp)->UpdateLigs(m_spLigand);
  filteridx = 0;
  nruns = 1;
  contextp->Assign(GetMetaDataPrefix() + "score.NRUNS", nruns);
}

// Called by Update when either model has changed
void Filter::SetupScore() {
  ((StringContextPtr)contextp)
      ->UpdateScores(GetWorkSpace()->GetSF(), GetWorkSpace()->GetModel(1));
  // write down rxdock.score.NRUNS to make sure is getting the
  // right value
  //  std::exit(1);
}

// Finished with ligand?
bool Filter::Terminate() {
  SetupScore();
  bool bTerm;
  if (nTermFilters > 0) {
    EvaluateVisitor visitor2(contextp);
    terminationFilters[filteridx]->Accept(visitor2);
    LOG_F(1, "{}\t{}\tnruns: {}", filteridx,
          terminationFilters[filteridx]->GetValue(), nruns);
    double val = terminationFilters[filteridx]->GetValue();
    if (val == STOP) {
      LOG_F(INFO, "Terminate with this ligand");
      bTerm = true;
    } else if (val == NEXT) {
      // keep checking filters until no more or one of them is not passed
      // by this conformation. I do not check if it return stop
      // or cont, just next filter.
      while ((val == NEXT) && (filteridx < nTermFilters)) {
        filteridx++;
        if (filteridx < nTermFilters) {
          nruns = 0; // it should not stop because of NRUNS
          contextp->Assign(GetMetaDataPrefix() + "score.NRUNS", nruns);
          terminationFilters[filteridx]->Accept(visitor2);
          val = terminationFilters[filteridx]->GetValue();
          LOG_F(INFO, "Go to next phase");
        }
      }
      if (filteridx == nTermFilters) {
        LOG_F(INFO, "Terminate with this ligand");
        bTerm = true; // output of filter is next phase, but there are no
                      // more phases
      } else {
        nruns = 1;
        contextp->Assign(GetMetaDataPrefix() + "score.NRUNS", nruns);
        bTerm = false;
      }
    } else if (val == CONT) {
      // finish because it's reached limit runs
      if (nruns >= maxnruns) {
        LOG_F(INFO, "Terminate with this ligand");
        bTerm = true;
      } else {
        nruns++;
        contextp->Assign(GetMetaDataPrefix() + "score.NRUNS", nruns);
        LOG_F(INFO, "Continue in this phase");
        bTerm = false;
      }
    } else {
      std::ostringstream error;
      error << "Wrong output: " << val;
      error << "; Termination filters should return 0,1 or -1" << std::ends;
      throw Error(_WHERE_, error.str());
    }
  } else
    bTerm = true; // if there are no Term filters, no need to write anything
                  // Simply terminate.
  return bTerm;
}

// Output conformation?
bool Filter::Write() {

  bool bWrite = true;
  for (int i = 0; i < nWriteFilters; i++) {
    EvaluateVisitor visitor2(contextp);
    writtingFilter[i]->Accept(visitor2);
    LOG_F(1, "Filter::Write: {}", writtingFilter[i]->GetValue());
    double val = writtingFilter[i]->GetValue();
    if (val >= 0.0) {
      bWrite = false;
      break;
    }
    // else it evaluates next filter. If there are no more, all
    // have returned a negative value and it returns true
  }
  if (bWrite) {
    LOG_F(INFO, "Write this ligand");
  } else {
    LOG_F(INFO, "Do not write this ligand");
  }
  return bWrite;
}
