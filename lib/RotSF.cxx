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

#include "rxdock/RotSF.h"

#include <loguru.hpp>

#include <functional>

using namespace rxdock;

// Static data member for class type
std::string RotSF::_CT("RotSF");
std::string RotSF::_INCNH3("INCNH3");
std::string RotSF::_INCOH("INCOH");

RotSF::RotSF(const std::string &strName)
    : BaseSF(_CT, strName), nRot(0), bIncNH3(false), bIncOH(false) {
  LOG_F(2, "RotSF parameterised constructor");
  // Add parameters
  AddParameter(_INCNH3, bIncNH3);
  AddParameter(_INCOH, bIncOH);
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RotSF::~RotSF() {
  LOG_F(2, "RotSF destructor");
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

void RotSF::SetupReceptor() {}

void RotSF::SetupLigand() {
  if (GetLigand().Null()) {
    nRot = 0;
  } else {
    BondList rotatableBondList =
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

void RotSF::SetupScore() {}

double RotSF::RawScore() const { return double(nRot); }

void RotSF::ParameterUpdated(const std::string &strName) {
  if (strName == _INCNH3) {
    bIncNH3 = GetParameter(_INCNH3);
  } else if (strName == _INCOH) {
    bIncOH = GetParameter(_INCOH);
  } else {
    BaseSF::ParameterUpdated(strName);
  }
}
