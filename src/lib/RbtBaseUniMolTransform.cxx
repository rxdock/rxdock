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

#include "RbtBaseUniMolTransform.h"
#include "RbtWorkSpace.h"

using namespace rxdock;

// Static data members
std::string RbtBaseUniMolTransform::_CT("RbtBaseUniMolTransform");

////////////////////////////////////////
// Constructors/destructors
RbtBaseUniMolTransform::RbtBaseUniMolTransform(const std::string &strClass,
                                               const std::string &strName)
    : RbtBaseTransform(strClass, strName) {
#ifdef _DEBUG
  std::cout << _CT << " parameterised constructor" << std::endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtBaseUniMolTransform::~RbtBaseUniMolTransform() {
#ifdef _DEBUG
  std::cout << _CT << " destructor" << std::endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

////////////////////////////////////////
// Public methods
////////////////

RbtModelPtr RbtBaseUniMolTransform::GetLigand() const { return m_spLigand; }

// Override RbtObserver pure virtual
// Notify observer that subject has changed
void RbtBaseUniMolTransform::Update(RbtSubject *theChangedSubject) {
  RbtWorkSpace *pWorkSpace = GetWorkSpace();
  if (theChangedSubject == pWorkSpace) {
    // Check if ligand has been updated (model #1)
    if (pWorkSpace->GetNumModels() >= 2) {
      RbtModelPtr spLigand = GetWorkSpace()->GetModel(1);
      if (spLigand != m_spLigand) {
#ifdef _DEBUG
        std::cout << _CT << "::Update(): Ligand has been updated" << std::endl;
#endif //_DEBUG
        m_spLigand = spLigand;
        SetupTransform();
      }
    }
  }
}
