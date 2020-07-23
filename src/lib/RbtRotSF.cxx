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

#include "RbtRotSF.h"

#include <functional>

using namespace rxdock;

// Static data member for class type
std::string RbtRotSF::_CT("RbtRotSF");
std::string RbtRotSF::_INCNH3("INCNH3");
std::string RbtRotSF::_INCOH("INCOH");

RbtRotSF::RbtRotSF(const std::string &strName)
    : RbtBaseSF(_CT, strName), nRot(0), bIncNH3(false), bIncOH(false) {
  // Add parameters
  AddParameter(_INCNH3, bIncNH3);
  AddParameter(_INCOH, bIncOH);
#ifdef _DEBUG
  std::cout << _CT << " parameterised constructor" << std::endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtRotSF::~RbtRotSF() {
#ifdef _DEBUG
  std::cout << _CT << " destructor" << std::endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

void RbtRotSF::SetupReceptor() {}

void RbtRotSF::SetupLigand() {
  if (GetLigand().Null()) {
    nRot = 0;
  } else {
    RbtBondList rotatableBondList =
        GetBondListWithPredicate(GetLigand()->GetBondList(), isBondRotatable());
    if (!bIncNH3) {
      rotatableBondList =
          GetBondListWithPredicate(rotatableBondList, std::not1(isBondToNH3()));
    }
    if (!bIncOH) {
      rotatableBondList =
          GetBondListWithPredicate(rotatableBondList, std::not1(isBondToOH()));
    }
    nRot = rotatableBondList.size();
  }
}

void RbtRotSF::SetupScore() {}

double RbtRotSF::RawScore() const { return double(nRot); }

void RbtRotSF::ParameterUpdated(const std::string &strName) {
  if (strName == _INCNH3) {
    bIncNH3 = GetParameter(_INCNH3);
  } else if (strName == _INCOH) {
    bIncOH = GetParameter(_INCOH);
  } else {
    RbtBaseSF::ParameterUpdated(strName);
  }
}
