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

#include "rxdock/ConstSF.h"

#include <loguru.hpp>

using namespace rxdock;

// Static data member for class type
const std::string ConstSF::_CT = "ConstSF";
const std::string ConstSF::_SOLVENT_PENALTY = "solvent-penalty";

ConstSF::ConstSF(const std::string &strName)
    : BaseSF(_CT, strName), m_solventPenalty(0.5) {
  LOG_F(2, "ConstSF parameterised constructor");
  AddParameter(_SOLVENT_PENALTY, m_solventPenalty);
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

ConstSF::~ConstSF() {
  LOG_F(2, "ConstSF destructor");
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

void ConstSF::ScoreMap(StringVariantMap &scoreMap) const {
  if (isEnabled()) {
    // Divide the total raw score into "system" and "inter" components.
    double rs = InterScore();

    // First deal with the inter score which is stored in its natural location
    // in the map
    std::string name = GetFullName();
    scoreMap[name] = rs;
    AddToParentMapEntry(scoreMap, rs);

    // Now deal with the system raw score which needs to be stored in
    // rxdock.score.system.const
    double system_rs = SystemScore();
    std::string systemName = BaseSF::_SYSTEM_SF + "." + GetName();
    scoreMap[systemName] = system_rs;
    // increment the rxdock.score.system total
    double parentScore = scoreMap[BaseSF::_SYSTEM_SF];
    parentScore += system_rs * GetWeight();
    scoreMap[BaseSF::_SYSTEM_SF] = parentScore;
  }
}

double ConstSF::RawScore() const { return InterScore() + SystemScore(); }

double ConstSF::InterScore() const {
  // With current protocols, ligand should always be enabled, so RawScore is 1
  double retVal = GetLigand()->GetEnabled() ? 1.0 : 0.0;
  return retVal;
}

double ConstSF::SystemScore() const {
  double retVal = 0.0;
  ModelList solventList = GetSolvent();
  for (ModelListConstIter iter = solventList.begin(); iter != solventList.end();
       ++iter) {
    if ((*iter)->GetEnabled()) {
      retVal += m_solventPenalty;
    }
  }
  return retVal;
}

void ConstSF::ParameterUpdated(const std::string &strName) {
  if (strName == _SOLVENT_PENALTY) {
    m_solventPenalty = GetParameter(_SOLVENT_PENALTY);
  } else {
    BaseSF::ParameterUpdated(strName);
  }
}

void rxdock::to_json(json &j, const ConstSF &consf) {
  j = json{{"solv-penalty", consf.m_solventPenalty}};
}

void rxdock::from_json(const json &j, ConstSF &consf) {
  j.at("solv-penalty").get_to(consf.m_solventPenalty);
}