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

#include "DihedralIntraSF.h"

using namespace rxdock;

// Static data members
std::string DihedralIntraSF::_CT("DihedralIntraSF");

DihedralIntraSF::DihedralIntraSF(const std::string &strName)
    : BaseSF(_CT, strName) {
#ifdef _DEBUG
  std::cout << _CT << " parameterised constructor" << std::endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

DihedralIntraSF::~DihedralIntraSF() {
  ClearModel();
#ifdef _DEBUG
  std::cout << _CT << " destructor" << std::endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

void DihedralIntraSF::SetupScore() {
  ClearModel();
  if (GetLigand().Null())
    return;
  if (GetLigand()->isFlexible()) {
    m_dihList = CreateDihedralList(GetLigand()->GetFlexBonds());
  }
}

double DihedralIntraSF::RawScore() const {
  double score = 0.0; // Total score
  for (DihedralListConstIter iter = m_dihList.begin(); iter != m_dihList.end();
       iter++) {
    score += (**iter)();
  }
  return score;
}

// Clear the dihedral list
// As we are not using smart pointers, there is some memory management to do
void DihedralIntraSF::ClearModel() {
  for (DihedralListIter iter = m_dihList.begin(); iter != m_dihList.end();
       iter++) {
    delete *iter;
  }
  m_dihList.clear();
}
