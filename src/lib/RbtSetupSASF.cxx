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

#include "RbtSetupSASF.h"

using namespace rxdock;

std::string RbtSetupSASF::_CT("RbtSetupSASF");

RbtSetupSASF::RbtSetupSASF(const std::string &strName)
    : RbtBaseSF(_CT, strName) {
#ifdef _DEBUG
  std::cout << _CT << "parameterized constructor      <--------" << std::endl;
#endif //_DEBUG
  Disable();
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtSetupSASF::~RbtSetupSASF() {
#ifdef _DEBUG
  std::cout << _CT << " destructor                    <--------" << std::endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

void RbtSetupSASF::SetupReceptor() {
#ifdef _DEBUG
  std::cout << _CT << " SetupReceptor                 <--------" << std::endl;
#endif //_DEBUG
       // get receptor atoms
       // theReceptorList =
       // GetAtomListWithPredicate(GetReceptor()->GetAtomList(),std::not1(isAtomicNo_eq(1)));
  theReceptorList = GetReceptor()->GetAtomList();
#ifdef _DEBUG
  std::cout << _CT << "::SetupReceptor(): #ATOMS = " << theReceptorList.size()
            << std::endl;
#endif //_DEBUG
}

void RbtSetupSASF::SetupScore() {
#ifdef _DEBUG
  std::cout << _CT << " SetupScore                    <--------" << std::endl;
#endif //_DEBUG
}

void RbtSetupSASF::SetupLigand() {
#ifdef _DEBUG
  std::cout << _CT << " SetupLigand                   <--------" << std::endl;
#endif //_DEBUG
  theLigandList = GetLigand()->GetAtomList();
#ifdef _DEBUG
  std::cout << _CT << "::SetupLigand(): #ATOMS = " << theLigandList.size()
            << std::endl;
#endif //_DEBUG
}

double RbtSetupSASF::RawScore() const { return 0.0; }
