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

#include "RbtSetupPolarSF.h"
#include <functional>
#include <iomanip>

using namespace rxdock;

// Static data members
std::string RbtSetupPolarSF::_CT("RbtSetupPolarSF");
std::string RbtSetupPolarSF::_RADIUS("RADIUS");
std::string RbtSetupPolarSF::_NORM("NORM");
std::string RbtSetupPolarSF::_POWER("POWER");
std::string RbtSetupPolarSF::_CHGFACTOR("CHGFACTOR");
std::string RbtSetupPolarSF::_GUANFACTOR("GUANFACTOR");

// NB - Virtual base class constructor (RbtBaseSF) gets called first,
// implicit constructor for RbtBaseInterSF is called second
RbtSetupPolarSF::RbtSetupPolarSF(const std::string &strName)
    : RbtBaseSF(_CT, strName) {
  // Add parameters
  AddParameter(_RADIUS, 5.0);
  AddParameter(_NORM, 25);
  AddParameter(_POWER, 0.5);
  AddParameter(_CHGFACTOR, 0.5);
  AddParameter(_GUANFACTOR, 1.0);
  Disable();
#ifdef _DEBUG
  std::cout << _CT << " parameterised constructor" << std::endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtSetupPolarSF::~RbtSetupPolarSF() {
#ifdef _DEBUG
  std::cout << _CT << " destructor" << std::endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

void RbtSetupPolarSF::SetupReceptor() {
  if (GetReceptor().Null())
    return;
  RbtAtomList atomList = GetReceptor()->GetAtomList();
  RbtAtomList heavyAtomList =
      GetAtomListWithPredicate(atomList, std::not1(isAtomicNo_eq(1)));
  int traceTriggerLevel = 1;
  SetupAtomList(atomList, heavyAtomList, traceTriggerLevel);
}

void RbtSetupPolarSF::SetupLigand() {
  if (GetLigand().Null())
    return;
  RbtAtomList atomList = GetLigand()->GetAtomList();
  RbtAtomList emptyList;
  int traceTriggerLevel = 2;
  SetupAtomList(atomList, emptyList, traceTriggerLevel);
}

void RbtSetupPolarSF::SetupSolvent() {
  RbtModelList solventList = GetSolvent();
  if (solventList.empty())
    return;
  // Concatenate all solvent atoms into a single list
  RbtAtomList solventAtomList;
  for (RbtModelListConstIter iter = solventList.begin();
       iter != solventList.end(); ++iter) {
    RbtAtomList atomList = (*iter)->GetAtomList();
    std::copy(atomList.begin(), atomList.end(),
              std::back_inserter(solventAtomList));
  }
  RbtAtomList emptyList;
  int traceTriggerLevel = 1;
  SetupAtomList(solventAtomList, emptyList, traceTriggerLevel);
}

// Reusable method for setting up an arbitrary atom list
// If neigbourList is empty, then do not include the neighbour density factor
void RbtSetupPolarSF::SetupAtomList(RbtAtomList &atomList,
                                    const RbtAtomList &neighbourList,
                                    int traceTriggerLevel) {
  if (atomList.empty())
    return;

  int iTrace = GetTrace();
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

  std::ios_base::fmtflags oldflags = std::cout.flags();
  std::cout.precision(3);
  std::cout.setf(std::ios_base::fixed, std::ios_base::floatfield);
  std::cout.setf(std::ios_base::right, std::ios_base::adjustfield);
  if (iTrace > traceTriggerLevel) {
    std::cout << std::endl
              << _CT
              << ": Atom\t\t, f(neighb), f(charge), isLipo,  isHBA,  isHBD"
              << std::endl;
  }

  for (RbtAtomListIter iter = atomList.begin(); iter != atomList.end();
       iter++) {
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
    // if (traceTriggerLevel == 1) {
    //  RbtDouble wxb = (*iter)->GetReweight();
    //  (*iter)->SetUser1Value(fNeighb*charge*wxb);
    // std::cout << "Atom: " << (*iter)->GetName() << " Polar_User1Value: "
    // <<
    // (*iter)->GetUser1Value() << " weigth: " << wxb << std::endl;
    //}else{
    (*iter)->SetUser1Value(fNeighb * charge);
    // }
    // XB END MODIFICATIONS
    (*iter)->SetUser1Flag(bIsLipo(*iter));
    if (iTrace > traceTriggerLevel) {
      //   std::cout << _CT << ": " << (*iter)->GetFullAtomName() << "\t," <<
      //   std::setw(10) << fNeighb << "," << std::setw(10) << charge
      //	   << "," << std::setw(7) << (*iter)->GetUser1Flag() << "," <<
      // std::setw(7) << bIsHBA(*iter) << "," << std::setw(7) << bIsHBD(*iter)
      // << std::endl;
    }
  }
  std::cout.flags(oldflags);
}

void RbtSetupPolarSF::SetupScore() {
  // No further setup required
}

double RbtSetupPolarSF::RawScore() const { return 0.0; }
