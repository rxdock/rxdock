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

#include "BaseBiMolTransform.h"
#include "WorkSpace.h"

using namespace rxdock;

// Static data members
std::string BaseBiMolTransform::_CT("BaseBiMolTransform");

////////////////////////////////////////
// Constructors/destructors
BaseBiMolTransform::BaseBiMolTransform(const std::string &strClass,
                                       const std::string &strName)
    : BaseTransform(strClass, strName) {
#ifdef _DEBUG
  std::cout << _CT << " parameterised constructor" << std::endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

BaseBiMolTransform::~BaseBiMolTransform() {
#ifdef _DEBUG
  std::cout << _CT << " destructor" << std::endl;
#endif //_DEBUG
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
#ifdef _DEBUG
        std::cout << _CT << "::Update(): Receptor has been updated"
                  << std::endl;
#endif //_DEBUG
        m_spReceptor = spReceptor;
        SetupReceptor();
      }
    }
    // Check if ligand has been updated (model #1)
    if (numModels >= 2) {
      ModelPtr spLigand = pWorkSpace->GetModel(1);
      if (spLigand != m_spLigand) {
#ifdef _DEBUG
        std::cout << _CT << "::Update(): Ligand has been updated" << std::endl;
#endif //_DEBUG
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
