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

#include "RbtFilter.h"
#include "RbtParser.h"
#include "RbtStringTokenIter.h"
#include "RbtWorkSpace.h"

using namespace rxdock;

const int STOP = 0;
const int CONT = -1;
const int NEXT = 1;

// Static data members
std::string RbtFilter::_CT("RbtFilter");
////////////////////////////////////////
// Constructors/destructors

// filter is used to determining if strfilter actually contains the
// filter (filter = true) or if strfilter is the name of a file
// that contains the filter (filter = false) This is the most
// common, so filter's value by default is false
RbtFilter::RbtFilter(std::string strfilter, bool filter)
    : RbtBaseObject(_CT, "Filter") {
#ifdef _DEBUG
  std::cout << _CT << " default constructor" << std::endl;
#endif //_DEBUG
       //  std::string filterfilen = GetParameter("_FILTER_FILE");
  SmartPtr<std::istream> filterfile;
  if (filter) // filterfilen is the filter
    filterfile = new std::istringstream(strfilter);
  else // assume the strfilter is the name of the file
       // where the filter is
    filterfile = new std::ifstream(strfilter.c_str(), std::ios_base::in);
  (*filterfile) >> nTermFilters;
  filteridx = 0;
  contextp = RbtContextPtr(new RbtStringContext(filterfile));
  RbtParser p;
  for (int i = 0; i < nTermFilters; i++) {
    std::cout << "\n------------- Terminate filter " << i << "------------"
              << std::endl;
    std::string s;
    getline(*filterfile, s, ',');
    SmartPtr<std::istream> istrp(new std::istringstream(s));
    RbtTokenIterPtr ti(new RbtStringTokenIter(istrp, contextp));
    RbtFilterExpressionPtr filterExpr = p.Parse(ti, contextp);
    PrettyPrintVisitor visitor1(contextp);
    filterExpr->Accept(visitor1);
    terminationFilters.push_back(filterExpr);
  }
  (*filterfile) >> nWriteFilters;
  for (int i = 0; i < nWriteFilters; i++) {
    std::cout << "\n------------- Write filter -----------------" << std::endl;
    std::string s;
    getline(*filterfile, s, ',');
    SmartPtr<std::istream> istrp(new std::istringstream(s));
    RbtTokenIterPtr ti(new RbtStringTokenIter(istrp, contextp));
    RbtFilterExpressionPtr filterExpr = p.Parse(ti, contextp);
    PrettyPrintVisitor visitor1(contextp);
    filterExpr->Accept(visitor1);
    writtingFilter.push_back(filterExpr);
  }
  maxnruns = 1000;
  std::cout << std::endl;
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtFilter::~RbtFilter() {
#ifdef _DEBUG
  std::cout << _CT << "  destructor" << std::endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

////////////////////////////////////////
// Public methods
////////////////

RbtModelPtr RbtFilter::GetReceptor() const { return m_spReceptor; }
RbtModelPtr RbtFilter::GetLigand() const { return m_spLigand; }

// Override RbtObserver pure virtual
// Notify observer that subject has changed
void RbtFilter::Update(RbtSubject *theChangedSubject) {
  RbtWorkSpace *pWorkSpace = GetWorkSpace();
  if (theChangedSubject == pWorkSpace) {
    // Check if receptor has been updated (model #0)
    if (pWorkSpace->GetNumModels() >= 1) {
      RbtModelPtr spReceptor = GetWorkSpace()->GetModel(0);
      if (spReceptor != m_spReceptor) {
#ifdef _DEBUG
        std::cout << "RbtBaseInterSF::Update(): Receptor has been updated"
                  << std::endl;
#endif //_DEBUG
        m_spReceptor = spReceptor;
        SetupReceptor();
      }
    }
    // Check if ligand has been updated (model #1)
    if (pWorkSpace->GetNumModels() >= 2) {
      RbtModelPtr spLigand = GetWorkSpace()->GetModel(1);
      if (spLigand != m_spLigand) {
#ifdef _DEBUG
        std::cout << "RbtBaseInterSF::Update(): Ligand has been updated"
                  << std::endl;
#endif //_DEBUG
        m_spLigand = spLigand;
        SetupLigand();
      }
    }
  }
}

// Called by Update when receptor is changed
void RbtFilter::SetupReceptor() {
  // read filter file
  ((RbtStringContextPtr)contextp)
      ->UpdateSite(m_spReceptor, GetWorkSpace()->GetDockingSite());
}

// Called by Update when ligand is changed
void RbtFilter::SetupLigand() {
  ((RbtStringContextPtr)contextp)->UpdateLigs(m_spLigand);
  filteridx = 0;
  nruns = 1;
  contextp->Assign("SCORE.NRUNS", nruns);
}

// Called by Update when either model has changed
void RbtFilter::SetupScore() {
  ((RbtStringContextPtr)contextp)
      ->UpdateScores(GetWorkSpace()->GetSF(), GetWorkSpace()->GetModel(1));
  // write down SCORE.NRUNS to make sure is getting the
  // right value
  //  std::exit(1);
}

// Finished with ligand?
bool RbtFilter::Terminate() {
  SetupScore();
  bool bTerm;
  if (nTermFilters > 0) {
    EvaluateVisitor visitor2(contextp);
    terminationFilters[filteridx]->Accept(visitor2);
    //    std::cout << filteridx << "\t"
    //         << terminationFilters[filteridx]->GetValue()
    //         << "\tnruns: " << nruns << std::endl;
    double val = terminationFilters[filteridx]->GetValue();
    if (val == STOP) {
      if (GetTrace() > 1) {
        std::cout << "Terminate with this ligand\n";
      }
      bTerm = true;
    } else if (val == NEXT) {
      // keep checking filters until no more or one of them is not passed
      // by this conformation. I do not check if it return stop
      // or cont, just next filter.
      while ((val == NEXT) && (filteridx < nTermFilters)) {
        filteridx++;
        if (filteridx < nTermFilters) {
          nruns = 0; // it should not stop because of NRUNS
          contextp->Assign("SCORE.NRUNS", nruns);
          terminationFilters[filteridx]->Accept(visitor2);
          val = terminationFilters[filteridx]->GetValue();
          if (GetTrace() > 1) {
            std::cout << "Go to next phase\n";
          }
        }
      }
      if (filteridx == nTermFilters) {
        if (GetTrace() > 1) {
          std::cout << "Terminate with this ligand\n";
        }
        bTerm = true; // output of filter is next phase, but there are no
                      // more phases
      } else {
        nruns = 1;
        contextp->Assign("SCORE.NRUNS", nruns);
        bTerm = false;
      }
    } else if (val == CONT) {
      // finish because it's reached limit runs
      if (nruns >= maxnruns) {
        if (GetTrace() > 1) {
          std::cout << "Terminate with this ligand\n";
        }
        bTerm = true;
      } else {
        nruns++;
        contextp->Assign("SCORE.NRUNS", nruns);
        if (GetTrace() > 1) {
          std::cout << "Continue in this phase\n";
        }
        bTerm = false;
      }
    } else {
      std::ostringstream error;
      error << "Wrong output: " << val;
      error << "; Termination filters should return 0,1 or -1" << std::ends;
      throw RbtError(_WHERE_, error.str());
    }
  } else
    bTerm = true; // if there are no Term filters, no need to write anything
                  // Simply terminate.
  return bTerm;
}

// Output conformation?
bool RbtFilter::Write() {

  bool bWrite = true;
  for (int i = 0; i < nWriteFilters; i++) {
    EvaluateVisitor visitor2(contextp);
    writtingFilter[i]->Accept(visitor2);
    //    std::cout << writtingFilter[i]->GetValue() << std::endl;
    double val = writtingFilter[i]->GetValue();
    if (val >= 0.0) {
      bWrite = false;
      break;
    }
    // else it evaluates next filter. If there are no more, all
    // have returned a negative value and it returns true
  }
  // if (bWrite)
  //  std::cout << "write this ligand\n";
  // else
  //  std::cout << "Do not write this ligand\n";
  return bWrite;
}
