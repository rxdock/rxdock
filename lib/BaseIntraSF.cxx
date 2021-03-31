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

#include "rxdock/BaseIntraSF.h"
#include "rxdock/WorkSpace.h"

#include <loguru.hpp>

using namespace rxdock;

// Static data members
const std::string BaseIntraSF::_CT = "BaseIntraSF";

BaseIntraSF::BaseIntraSF() : m_zero(0.0) {
  LOG_F(2, "BaseIntraSF default constructor");
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

BaseIntraSF::~BaseIntraSF() {
  LOG_F(2, "BaseIntraSF destructor");
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

////////////////////////////////////////
// Public methods
////////////////

ModelPtr BaseIntraSF::GetLigand() const { return m_spLigand; }

// Override Observer pure virtual
// Notify observer that subject has changed
void BaseIntraSF::Update(Subject *theChangedSubject) {
  WorkSpace *pWorkSpace = GetWorkSpace();
  if (theChangedSubject == pWorkSpace) {
    // Check if ligand has been updated (model #1)
    if (pWorkSpace->GetNumModels() >= 2) {
      ModelPtr spLigand = GetWorkSpace()->GetModel(1);
      if (spLigand != m_spLigand) {
        LOG_F(1, "BaseIntraSF::Update(): Ligand has been updated");
        m_spLigand = spLigand;
        SetupScore();
        // Retain the zero-point offset from the ligand model data if present
        // Otherwise set the zero-point to the current score
        std::string name = GetFullName() + ".0";
        if (m_spLigand->isDataFieldPresent(name)) {
          m_zero = m_spLigand->GetDataValue(name);
        } else {
          m_zero = RawScore();
        }
      }
    }
  }
}

// Override BaseSF::ScoreMap to provide additional raw descriptors
void BaseIntraSF::ScoreMap(StringVariantMap &scoreMap) const {
  if (isEnabled()) {
    double rs =
        RawScore() -
        m_zero; // report the raw score relative to the zero-point offset
    std::string name = GetFullName();
    scoreMap[name] = rs;
    AddToParentMapEntry(scoreMap, rs);
    scoreMap[name + ".0"] = m_zero;
  }
}
