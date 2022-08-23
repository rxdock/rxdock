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

#include "rxdock/BaseSF.h"
#include "rxdock/SFRequest.h"

#include <loguru.hpp>

using namespace rxdock;

// Static data members
const std::string BaseSF::_CT = "BaseSF";
const std::string BaseSF::_WEIGHT = "weight";
const std::string BaseSF::_RANGE = "range";
const std::string BaseSF::_SYSTEM_SF = "rxdock.score.system";
const std::string BaseSF::_INTRA_SF = "rxdock.score.intra";

////////////////////////////////////////
// Constructors/destructors
BaseSF::BaseSF(const std::string &strClass, const std::string &strName)
    : BaseObject(strClass, strName), m_parent(nullptr), m_weight(1.0),
      m_range(10.0) {
  LOG_F(2, "BaseSF parameterised constructor for {}", strClass);
  // Add parameters
  AddParameter(_WEIGHT, m_weight);
  AddParameter(_RANGE, m_range);
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

// Dummy default constructor for virtual base subclasses
// Should never get called
BaseSF::BaseSF() {
  LOG_F(WARNING, "BaseSF default constructor");
  //_RBTOBJECTCOUNTER_CONSTR_(_CT);
}

BaseSF::~BaseSF() {
  LOG_F(2, "BaseSF destructor");
  Orphan(); // Remove object from parent aggregate
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

////////////////////////////////////////
// Public methods
////////////////

// Fully qualified name, prefixed by all ancestors (e.g.
// rxdock.score.inter.hbond)
std::string BaseSF::GetFullName() const {
  // If we have a parent, prefix our short name with our parents full name,
  // else full name is just our short name
  return (m_parent) ? m_parent->GetFullName() + "." + GetName() : GetName();
}

double BaseSF::GetWeight() const { return m_weight; }
void BaseSF::SetWeight(double weight) { SetParameter(_WEIGHT, weight); }

double BaseSF::GetRange() const { return m_range; }
void BaseSF::SetRange(double range) { SetParameter(_RANGE, range); }

// Returns weighted score if scoring function is enabled, else returns zero
double BaseSF::Score() const {
  return isEnabled() ? GetWeight() * RawScore() : 0.0;
}

// Returns all child component scores as a string-variant map
// Key = fully qualified component name, value = weighted score
//(for saving in a Model's data fields)
void BaseSF::ScoreMap(StringVariantMap &scoreMap) const {
  if (isEnabled()) {
    // DM 17 Jan 2006.
    // New approach:
    // 1) We record the raw, unweighted score for this term
    // in the map
    double rs = RawScore();
    std::string name = GetFullName();
    scoreMap[name] = rs;
    // 2) We add the weighted score to the parent aggregate
    // entry. This gives us the opportunity to override
    // ScoreMap in order to divert scores away from their
    // natural parent entry. e.g. rxdock.score.inter.vdw could
    // record intra-receptor and intra-solvent contributions
    // under rxdock.score.system.vdw
    AddToParentMapEntry(scoreMap, rs);
  }
}

// Helper method for ScoreMap
void BaseSF::AddToParentMapEntry(StringVariantMap &scoreMap, double rs) const {
  if (m_parent) {
    double w = GetWeight();
    double s = w * rs;
    std::string parentName = m_parent->GetFullName();
    double parentScore = scoreMap[parentName];
    parentScore += s;
    scoreMap[parentName] = parentScore;
  }
}

// Aggregate handling (virtual) methods
// Base class throws an InvalidRequest error

void BaseSF::Add(BaseSF *) {
  throw InvalidRequest(_WHERE_,
                       "Add() invalid for non-aggregate scoring functions");
}
void BaseSF::Remove(BaseSF *) {
  throw InvalidRequest(_WHERE_,
                       "Remove() invalid for non-aggregate scoring functions");
}

BaseSF *BaseSF::GetSF(unsigned int iSF) const {
  throw InvalidRequest(_WHERE_,
                       "GetSF() invalid for non-aggregate scoring functions");
}
bool BaseSF::isAgg() const { return false; }
unsigned int BaseSF::GetNumSF() const { return 0; }

// Aggregate handling (concrete) methods
BaseSF *BaseSF::GetParentSF() const { return m_parent; }
// Force removal from the parent aggregate
void BaseSF::Orphan() {
  if (m_parent) {
    LOG_F(1, "BaseSF::Orphan: Removing {} from {}", GetName(),
          m_parent->GetName());
    m_parent->Remove(this);
  }
}

// DM 25 Oct 2000 - track changes to parameter values in local data members
// ParameterUpdated is invoked by ParamHandler::SetParameter
void BaseSF::ParameterUpdated(const std::string &strName) {
  if (strName == _WEIGHT) {
    m_weight = GetParameter(_WEIGHT);
  } else if (strName == _RANGE) {
    m_range = GetParameter(_RANGE);
  } else {
    BaseObject::ParameterUpdated(strName);
  }
}

void rxdock::to_json(json &j, const BaseSF &baseSF) {
  // parent pointer skipped
  j = json{{"weight", baseSF.m_weight}, {"range", baseSF.m_range}};
}

void rxdock::from_json(const json &j, BaseSF &baseSF) {
  // parent pointer skipped
  j.at("weight").get_to(baseSF.m_weight);
  j.at("range").get_to(baseSF.m_range);
}
