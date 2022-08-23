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

#include "rxdock/BaseUniMolTransform.h"
#include "rxdock/WorkSpace.h"

#include <loguru.hpp>

using namespace rxdock;

// Static data members
const std::string BaseUniMolTransform::_CT = "BaseUniMolTransform";

////////////////////////////////////////
// Constructors/destructors
BaseUniMolTransform::BaseUniMolTransform(const std::string &strClass,
                                         const std::string &strName)
    : BaseTransform(strClass, strName) {
  LOG_F(2, "BaseUniMolTransform parameterised constructor");
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

BaseUniMolTransform::~BaseUniMolTransform() {
  LOG_F(2, "BaseUniMolTransform destructor")
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

////////////////////////////////////////
// Public methods
////////////////

ModelPtr BaseUniMolTransform::GetLigand() const { return m_spLigand; }

// Override Observer pure virtual
// Notify observer that subject has changed
void BaseUniMolTransform::Update(Subject *theChangedSubject) {
  WorkSpace *pWorkSpace = GetWorkSpace();
  if (theChangedSubject == pWorkSpace) {
    // Check if ligand has been updated (model #1)
    if (pWorkSpace->GetNumModels() >= 2) {
      ModelPtr spLigand = GetWorkSpace()->GetModel(1);
      if (spLigand != m_spLigand) {
        LOG_F(1, "BaseUniMolTransform::Update: Ligand has been updated");
        m_spLigand = spLigand;
        SetupTransform();
      }
    }
  }
}

void rxdock::to_json(json &j, const BaseUniMolTransform &baseUniMolTrans) {
  j = json{{"ligand", *baseUniMolTrans.m_spLigand}};
}

void rxdock::from_json(const json &j, BaseUniMolTransform &baseUniMolTrans) {
  j.at("ligand").get_to(*baseUniMolTrans.m_spLigand);
}