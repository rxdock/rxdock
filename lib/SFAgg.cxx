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

#include "rxdock/SFAgg.h"
#include "rxdock/Model.h"
#include "rxdock/WorkSpace.h"

#include <loguru.hpp>

#include <functional>

using namespace rxdock;

// Static data member for class type
const std::string SFAgg::_CT = "SFAgg";

////////////////////////////////////////
// Constructors/destructors
SFAgg::SFAgg(const std::string &strName)
    : BaseSF(_CT, strName), m_nNonHLigandAtoms(0) {
  LOG_F(2, "SFAgg parameterised constructor");
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

SFAgg::~SFAgg() {
  LOG_F(2, "SFAgg destructor");
  // Delete all our children
  // We need to iterate using a while loop because each deletion will
  // reduce the size of m_sf, hence conventional iterators would become
  // invalid
  LOG_F(1, "Deleting child scoring functions of {}", GetName());
  while (m_sf.size() > 0) {
    BaseSF *pSF = m_sf.back();
    // Assertion: parent of child is this object
    Assert<Assertion>(!SFAGG_CHECK || pSF->m_parent == this);
    LOG_F(1, "Deleting {} from {}", pSF->GetName(), GetName());
    delete pSF;
  }
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

////////////////////////////////////////
// Public methods
////////////////

// Returns all child component scores as a string-variant map
// Key = fully qualified component name, value = weighted score
//(for saving in a Model's data fields)
void SFAgg::ScoreMap(StringVariantMap &scoreMap) const {
  if (isEnabled()) {
    // First populate all the child entries in the scoring function map
    for (BaseSFListConstIter iter = m_sf.begin(); iter != m_sf.end(); iter++) {
      (*iter)->ScoreMap(scoreMap);
    }
    // Now we can pick up the raw score for this aggregate from the map entry
    // without calculating it directly,
    // as the child terms will have incremented the total as part of their
    // ScoreMap implementation
    std::string name = GetFullName();
    double rs = scoreMap[name];

    // Cascade to the parent of this aggregate
    AddToParentMapEntry(scoreMap, rs);

    // 7 Feb 2005 (DM, Enspiral Discovery) - normalise the raw aggregate score
    // by the number of heavy (non-H) ligand atoms
    if (m_nNonHLigandAtoms > 0) {
      double norm_s = rs / m_nNonHLigandAtoms;
      scoreMap[name + ".norm"] = norm_s;
      // Only record the number of heavy atoms for the root aggregate
      // (rxdock.score) root = has no parent
      if (!GetParentSF()) {
        scoreMap[name + ".heavy"] = m_nNonHLigandAtoms;
      }
    }
  }
}

// Aggregate handling methods
void SFAgg::Add(BaseSF *pSF) {
  // By first orphaning the scoring function to be added,
  // we handle attempts to readd existing children automatically,
  pSF->Orphan();
  pSF->m_parent = this;
  LOG_F(1, "SFAgg::Add: Adding {} to {}", pSF->GetName(), GetName());
  m_sf.push_back(pSF);
}

void SFAgg::Remove(BaseSF *pSF) {
  BaseSFListIter iter = std::find(m_sf.begin(), m_sf.end(), pSF);
  if (iter == m_sf.end()) {
    throw BadArgument(_WHERE_,
                      "Remove(): pSF is not a member of this aggregate");
  } else {
    // Assertion: parent of child is this object
    Assert<Assertion>(!SFAGG_CHECK || pSF->m_parent == this);
    LOG_F(2, "SFAgg::Remove: Removing {} from {}", pSF->GetName(), GetName());
    m_sf.erase(iter);
    pSF->m_parent = nullptr; // Nullify the parent pointer of the child that has
                             // been removed
  }
}

bool SFAgg::isAgg() const { return true; }
unsigned int SFAgg::GetNumSF() const { return m_sf.size(); }
BaseSF *SFAgg::GetSF(unsigned int iSF) const {
  if (iSF >= m_sf.size()) {
    throw BadArgument(_WHERE_, "GetSF(): iSF out of range");
  } else {
    return m_sf[iSF];
  }
}

// WorkSpace handling methods
// Register scoring function with a workspace
// Aggregate version registers all children, AND itself (new behaviour, 7 Feb
// 2005 (DM))
void SFAgg::Register(WorkSpace *pWorkSpace) {
  // 7 Feb 2005 (DM, Enspiral Discovery) - register the aggregate as well
  // as we need the number of heavy atoms in the ligand for calculating
  // the normalised scores
  BaseObject::Register(pWorkSpace);
  LOG_F(1, "SFAgg::Register: Registering child scoring functions of {}",
        GetName());
  for (BaseSFListIter iter = m_sf.begin(); iter != m_sf.end(); iter++) {
    (*iter)->Register(pWorkSpace);
  }
}

// Unregister with a workspace
// Aggregate version unregisters all children, AND itself (new behaviour, 7 Feb
// 2005 (DM))
void SFAgg::Unregister() {
  LOG_F(1, "SFAgg::Unregister: Unregistering child scoring functions of {}",
        GetName());
  for (BaseSFListIter iter = m_sf.begin(); iter != m_sf.end(); iter++) {
    (*iter)->Unregister();
  }
  BaseObject::Unregister(); // 7 Feb 2005 (DM) - unregister the agg, for the
                            // same reasons as above
}

// Override Observer pure virtual
// Notify observer that subject has changed
// 7 Feb 2005 (DM, Enspiral Discovery) - keep track of number of non-H ligand
// atoms so can report aggregate scores normalised by #heavy atoms
void SFAgg::Update(Subject *theChangedSubject) {
  WorkSpace *pWorkSpace = GetWorkSpace();
  if (theChangedSubject == pWorkSpace) {
    // Check if ligand has been updated (model #1)
    if (pWorkSpace->GetNumModels() >= 2) {
      ModelPtr spLigand = pWorkSpace->GetModel(1);
      if (spLigand.Ptr()) {
        m_nNonHLigandAtoms = GetNumAtomsWithPredicate(
            spLigand->GetAtomList(), std::not1(isAtomicNo_eq(1)));
      } else {
        m_nNonHLigandAtoms = 0;
      }
    }
  }
}

// Request Handling method
// Aggregates handle the request themselves first, then cascade to all children
void SFAgg::HandleRequest(RequestPtr spRequest) {
  BaseObject::HandleRequest(spRequest);
  for (BaseSFListIter iter = m_sf.begin(); iter != m_sf.end(); iter++) {
    (*iter)->HandleRequest(spRequest);
  }
}

// Virtual function for dumping scoring function details to an output stream
// Called by operator <<
void SFAgg::Print(std::ostream &s) const {
  // First print the parameters for this aggregate
  BaseObject::Print(s);
  // Now call Print() for each of the children
  for (BaseSFListConstIter iter = m_sf.begin(); iter != m_sf.end(); iter++) {
    (*iter)->Print(s);
  }
}

////////////////////////////////////////
// Private methods
/////////////////
// Raw score for an aggregate is the sum of the weighted scores of its children
double SFAgg::RawScore() const {
  double score(0.0);
  for (BaseSFListConstIter iter = m_sf.begin(); iter != m_sf.end(); iter++) {
    score += (*iter)->Score();
  }
  return score;
}
