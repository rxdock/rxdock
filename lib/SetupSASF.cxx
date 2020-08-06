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

#include "rxdock/SetupSASF.h"

#include <loguru.hpp>

using namespace rxdock;

std::string SetupSASF::_CT("SetupSASF");

SetupSASF::SetupSASF(const std::string &strName) : BaseSF(_CT, strName) {
  LOG_F(2, "SetupSASF parameterized constructor");
  Disable();
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

SetupSASF::~SetupSASF() {
  LOG_F(2, "SetupSASF destructor");
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

void SetupSASF::SetupReceptor() {
  LOG_F(2, "SetupSASF::SetupReceptor");
  // get receptor atoms
  // theReceptorList =
  // GetAtomListWithPredicate(GetReceptor()->GetAtomList(),
  // std::not1(isAtomicNo_eq(1)));
  theReceptorList = GetReceptor()->GetAtomList();
  LOG_F(1, "SetupSASF::SetupReceptor: #atoms = {}", theReceptorList.size());
}

void SetupSASF::SetupScore() { LOG_F(2, "SetupSASF::SetupScore"); }

void SetupSASF::SetupLigand() {
  LOG_F(2, "SetupSASF::SetupLigand");
  theLigandList = GetLigand()->GetAtomList();
  LOG_F(1, "SetupSASF::SetupLigand: #atoms= {}", theLigandList.size());
}

double SetupSASF::RawScore() const { return 0.0; }
