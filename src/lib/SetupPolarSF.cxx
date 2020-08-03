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

#include "SetupPolarSF.h"

#include <loguru.hpp>

#include <functional>
#include <iomanip>

using namespace rxdock;

// Static data members
std::string SetupPolarSF::_CT("SetupPolarSF");
std::string SetupPolarSF::_RADIUS("RADIUS");
std::string SetupPolarSF::_NORM("NORM");
std::string SetupPolarSF::_POWER("POWER");
std::string SetupPolarSF::_CHGFACTOR("CHGFACTOR");
std::string SetupPolarSF::_GUANFACTOR("GUANFACTOR");

// NB - Virtual base class constructor (BaseSF) gets called first,
// implicit constructor for BaseInterSF is called second
SetupPolarSF::SetupPolarSF(const std::string &strName) : BaseSF(_CT, strName) {
  LOG_F(2, "SetupPolarSF parameterised constructor");
  // Add parameters
  AddParameter(_RADIUS, 5.0);
  AddParameter(_NORM, 25);
  AddParameter(_POWER, 0.5);
  AddParameter(_CHGFACTOR, 0.5);
  AddParameter(_GUANFACTOR, 1.0);
  Disable();
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

SetupPolarSF::~SetupPolarSF() {
  LOG_F(2, "SetupPolarSF destructor");
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

void SetupPolarSF::SetupReceptor() {
  if (GetReceptor().Null())
    return;
  AtomList atomList = GetReceptor()->GetAtomList();
  AtomList heavyAtomList =
      GetAtomListWithPredicate(atomList, std::not1(isAtomicNo_eq(1)));
  SetupAtomList(atomList, heavyAtomList);
}

void SetupPolarSF::SetupLigand() {
  if (GetLigand().Null())
    return;
  AtomList atomList = GetLigand()->GetAtomList();
  AtomList emptyList;
  SetupAtomList(atomList, emptyList);
}

void SetupPolarSF::SetupSolvent() {
  ModelList solventList = GetSolvent();
  if (solventList.empty())
    return;
  // Concatenate all solvent atoms into a single list
  AtomList solventAtomList;
  for (ModelListConstIter iter = solventList.begin(); iter != solventList.end();
       ++iter) {
    AtomList atomList = (*iter)->GetAtomList();
    std::copy(atomList.begin(), atomList.end(),
              std::back_inserter(solventAtomList));
  }
  AtomList emptyList;
  SetupAtomList(solventAtomList, emptyList);
}

// Reusable method for setting up an arbitrary atom list
// If neigbourList is empty, then do not include the neighbour density factor
void SetupPolarSF::SetupAtomList(AtomList &atomList,
                                 const AtomList &neighbourList) {
  if (atomList.empty())
    return;

  bool bCalcNeighbourDensity = !neighbourList.empty();
  double radius = GetParameter(_RADIUS);
  double norm = GetParameter(_NORM);
  double power = GetParameter(_POWER);
  double chgFactor = GetParameter(_CHGFACTOR);
  double guanFactor = GetParameter(_GUANFACTOR);

  isAtomHBondAcceptor bIsHBA;
  isAtomHBondDonor bIsHBD;
  isAtomGuanidiniumCarbon bIsGuan;
  isAtomLipophilic bIsLipo;

  LOG_F(1,
        "SetupPolarSF: Atom\t\t, f(neighb), f(charge), isLipo,  isHBA,  isHBD");
  for (AtomListIter iter = atomList.begin(); iter != atomList.end(); iter++) {
    double fNeighb = 1.0;
    if (bCalcNeighbourDensity) {
      // Neighbour density function
      int nNeighb = GetNumAtomsWithPredicate(
          neighbourList, isAtomInsideSphere((*iter)->GetCoords(), radius));
      nNeighb--; // Don't count the atom itself!
      fNeighb = std::pow(nNeighb / norm, power);
    }
    // Charge factor (always >= 1, sign reflects whether charge is positive or
    // negative)
    double charge = (*iter)->GetGroupCharge();
    // DM 18 Jan 2001 - need to check for neutral HBond acceptors and set sign
    // to -1
    double sign = ((charge < 0.0) || (bIsHBA(*iter))) ? -1.0 : 1.0;
    charge = sign * (1.0 + (std::fabs(charge) * chgFactor));
    if (bIsGuan(*iter)) {
      charge *= guanFactor; // Adjustable weight for "IONIC" interactions with
                            // guanidinium carbons
    }
    // XB apply reweighting factor for receptor's polar atoms:
    // if (reweight)
    // double wxb = (*iter)->GetReweight();
    // (*iter)->SetUser1Value(fNeighb * charge * wxb);
    // LOG_F(1, "Atom: {} Polar_User1Value: {} weight: {}", (*iter)->GetName(),
    //       (*iter)->GetUser1Value(), (*iter)->GetReweight());
    // } else {
    (*iter)->SetUser1Value(fNeighb * charge);
    // }
    // XB END MODIFICATIONS
    (*iter)->SetUser1Flag(bIsLipo(*iter));
    LOG_F(1, "SetupPolarSF:{}\t,{:10.3f}{:10.3f},{:7d},{:7d},{:7d}",
          (*iter)->GetFullAtomName(), fNeighb, charge, (*iter)->GetUser1Flag(),
          bIsHBA(*iter), bIsHBD(*iter));
  }
}

void SetupPolarSF::SetupScore() {
  // No further setup required
}

double SetupPolarSF::RawScore() const { return 0.0; }
