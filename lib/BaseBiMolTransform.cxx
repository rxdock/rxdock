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

#include "rxdock/BaseBiMolTransform.h"
#include "rxdock/WorkSpace.h"

#include <loguru.hpp>

using namespace rxdock;

// Static data members
const std::string BaseBiMolTransform::_CT = "BaseBiMolTransform";

////////////////////////////////////////
// Constructors/destructors
BaseBiMolTransform::BaseBiMolTransform(const std::string &strClass,
                                       const std::string &strName)
    : BaseTransform(strClass, strName) {
  LOG_F(2, "BaseBiMolTransform parameterised constructor");
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

BaseBiMolTransform::~BaseBiMolTransform() {
  LOG_F(2, "BaseBiMolTransform destructor")
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

////////////////////////////////////////
// Public methods
////////////////

ModelPtr BaseBiMolTransform::GetReceptor() const { return m_spReceptor; }
ModelPtr BaseBiMolTransform::GetLigand() const { return m_spLigand; }
ModelList BaseBiMolTransform::GetSolvent() const { return m_solventList; }

// Override Observer pure virtual
// Notify observer that subject has changed
void BaseBiMolTransform::Update(Subject *theChangedSubject) {
  WorkSpace *pWorkSpace = GetWorkSpace();
  if (theChangedSubject == pWorkSpace) {
    int numModels = pWorkSpace->GetNumModels();
    // Check if receptor has been updated (model #0)
    if (numModels >= 1) {
      ModelPtr spReceptor = pWorkSpace->GetModel(0);
      if (spReceptor != m_spReceptor) {
        LOG_F(1, "BaseBiMolTransform::Update(): Receptor has been updated");
        m_spReceptor = spReceptor;
        SetupReceptor();
      }
    }
    // Check if ligand has been updated (model #1)
    if (numModels >= 2) {
      ModelPtr spLigand = pWorkSpace->GetModel(1);
      if (spLigand != m_spLigand) {
        LOG_F(1, "BaseBiMolTransform::Update(): Ligand has been updated");
        m_spLigand = spLigand;
        SetupLigand();
      }
    }
    // check if solvent has been updated (model #2 onwards)
    // At the moment, we only check whether the number of solvent
    // molecules has changed
    ModelList solventList;
    if (numModels >= 3) {
      solventList = pWorkSpace->GetModels(2);
    }
    if (solventList.size() != m_solventList.size()) {
      m_solventList = solventList;
      SetupSolvent();
    }
    SetupTransform();
  }
}

void rxdock::to_json(json &j, const BaseBiMolTransform &baseBiMolTrans) {
  json modelList;
  for (const auto &aIter : baseBiMolTrans.m_solventList) {
    json model = *aIter;
    modelList.push_back(model);
  }

  j = json{{"receptor", *baseBiMolTrans.m_spReceptor},
           {"ligand", *baseBiMolTrans.m_spLigand},
           {"model-list", modelList}};
}

void rxdock::from_json(const json &j, BaseBiMolTransform &baseBiMolTrans) {
  j.at("receptor").get_to(*baseBiMolTrans.m_spReceptor);
  j.at("ligand").get_to(*baseBiMolTrans.m_spLigand);

  for (auto &model : j.at("model-list")) {
    ModelPtr spModel = ModelPtr(new Model(model));
    baseBiMolTrans.m_solventList.push_back(spModel);
  }
}