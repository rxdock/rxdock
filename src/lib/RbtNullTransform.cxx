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

#include "RbtNullTransform.h"
#include "RbtPopulation.h"
#include "RbtWorkSpace.h"

// Static data member for class type
std::string RbtNullTransform::_CT("RbtNullTransform");

////////////////////////////////////////
// Constructors/destructors
RbtNullTransform::RbtNullTransform(const std::string &strName)
    : RbtBaseTransform(_CT, strName) {
#ifdef _DEBUG
  std::cout << _CT << " parameterised constructor" << std::endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtNullTransform::~RbtNullTransform() {
#ifdef _DEBUG
  std::cout << _CT << " destructor" << std::endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

////////////////////////////////////////
// Public methods
////////////////
void RbtNullTransform::Update(RbtSubject *theChangedSubject) {}

////////////////////////////////////////
// Private methods
////////////////
void RbtNullTransform::Execute() {
  int iTrace = GetTrace();
  // Trace level 1 = just output the current conformation (presumed to be the
  // best from the previous transform) to the history file
  if (iTrace == 5) {
    GetWorkSpace()->SaveHistory(true);
    double s = GetWorkSpace()->GetSF()->Score();
    std::cout << "SCORE = " << s << std::endl;
  }
  // Trace level 2 = output entire GA population to history file, unless pop
  // does not exist, in which case same as level 1
  else if (iTrace > 5) {
    RbtPopulationPtr pop = GetWorkSpace()->GetPopulation();
    if (pop.Ptr() == nullptr) {
      GetWorkSpace()->SaveHistory(true);
    } else {
      const RbtGenomeList &genList = pop->GetGenomeList();
      for (RbtGenomeListConstIter gIter = genList.begin();
           gIter != genList.end(); ++gIter) {
        (*gIter)->GetChrom()->SyncToModel();
        GetWorkSpace()->SaveHistory(true);
      }
    }
  }
}
