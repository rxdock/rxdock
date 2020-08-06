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

#include "rxdock/NullTransform.h"
#include "rxdock/Population.h"
#include "rxdock/WorkSpace.h"

#include <loguru.hpp>

using namespace rxdock;

// Static data member for class type
std::string NullTransform::_CT("NullTransform");

////////////////////////////////////////
// Constructors/destructors
NullTransform::NullTransform(const std::string &strName)
    : BaseTransform(_CT, strName) {
  LOG_F(2, "NullTransform parameterised constructor");
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

NullTransform::~NullTransform() {
  LOG_F(2, "NullTransform destructor");
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

////////////////////////////////////////
// Public methods
////////////////
void NullTransform::Update(Subject *theChangedSubject) {}

////////////////////////////////////////
// Private methods
////////////////
void NullTransform::Execute() {
  // Trace level 1 = just output the current conformation (presumed to be the
  // best from the previous transform) to the history file
  GetWorkSpace()->SaveHistory(true);
  LOG_F(1, "NullTransform::Execute: SCORE = {}",
        GetWorkSpace()->GetSF()->Score());
  // Trace level 2 = output entire GA population to history file, unless pop
  // does not exist, in which case same as level 1
  /* PopulationPtr pop = GetWorkSpace()->GetPopulation();
  if (pop.Ptr() == nullptr) {
    GetWorkSpace()->SaveHistory(true);
  } else {
    const GenomeList &genList = pop->GetGenomeList();
    for (GenomeListConstIter gIter = genList.begin(); gIter != genList.end();
         ++gIter) {
      (*gIter)->GetChrom()->SyncToModel();
      GetWorkSpace()->SaveHistory(true);
    }
  } */
}
