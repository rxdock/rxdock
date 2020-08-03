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

#include "DihedralTargetSF.h"

#include <loguru.hpp>

using namespace rxdock;

// Static data members
std::string DihedralTargetSF::_CT("DihedralTargetSF");

DihedralTargetSF::DihedralTargetSF(const std::string &strName)
    : BaseSF(_CT, strName) {
  LOG_F(2, "DihedralTargetSF parameterised constructor");
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

DihedralTargetSF::~DihedralTargetSF() {
  LOG_F(2, "DihedralTargetSF destructor");
  ClearReceptor();
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

void DihedralTargetSF::SetupReceptor() {
  ClearReceptor();
  if (GetReceptor().Null())
    return;
  if (GetReceptor()->isFlexible()) {
    m_dihList = CreateDihedralList(GetReceptor()->GetFlexBonds());
  }
}
void DihedralTargetSF::SetupLigand() {}
void DihedralTargetSF::SetupScore() {}

double DihedralTargetSF::RawScore() const {
  double score = 0.0; // Total score
  for (DihedralListConstIter iter = m_dihList.begin(); iter != m_dihList.end();
       iter++) {
    score += (**iter)();
  }
  return score;
}

// Clear the dihedral list
// As we are not using smart pointers, there is some memory management to do
void DihedralTargetSF::ClearReceptor() {
  for (DihedralListIter iter = m_dihList.begin(); iter != m_dihList.end();
       iter++) {
    delete *iter;
  }
  m_dihList.clear();
}
