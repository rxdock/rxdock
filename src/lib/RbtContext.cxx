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

#include "RbtContext.h"
#include "RbtDebug.h"
#include <cstring>
#include <functional>

using namespace rxdock;

std::string RbtContext::_CT("RbtContext");
// std::string RbtStringContext::_CT("RbtStringContext");
// std::string RbtCellContext::_CT("RbtCellContext");
std::string RbtVble::_CT("RbtVble");

RbtContext::RbtContext() {
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
#ifdef _DEBUG
  std::cout << _CT << "base context const\n";
#endif //_DEBUG
}

RbtContext::RbtContext(const RbtContext &c) {
  _RBTOBJECTCOUNTER_COPYCONSTR_(_CT);
}

RbtContext::~RbtContext() { _RBTOBJECTCOUNTER_DESTR_(_CT); }

RbtStringContext::RbtStringContext() {
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
  throw RbtError(_WHERE_, "Can't construct string context by default");
}

RbtStringContext::RbtStringContext(const RbtStringContext &c)
    : RbtContext(c), vm(c.vm) {
  _RBTOBJECTCOUNTER_COPYCONSTR_(_CT);
}

RbtStringContext::~RbtStringContext() { _RBTOBJECTCOUNTER_DESTR_(_CT); }

RbtStringContext::RbtStringContext(SmartPtr<std::ifstream> ifile) {
  /*    RbtInt nvbles, nctes;
      (*ifile) >> nvbles >> nctes;
      std::string name;
      RbtInt i = 0;
      for (; i < nvbles ; i++)
      {
          (*ifile) >> name;
          vm[name] = new RbtVble(name);
      }
      for (; i < nctes + nvbles ; i++)
      {
          (*ifile) >> name;
          RbtDouble val = std::atof(name.c_str());
          vm[name] = new RbtVble(name, val);
      } */
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtCellContext::RbtCellContext() { _RBTOBJECTCOUNTER_CONSTR_(_CT); }

RbtCellContext::RbtCellContext(const RbtCellContext &c)
    : RbtContext(c), vm(c.vm) {
  _RBTOBJECTCOUNTER_COPYCONSTR_(_CT);
}

RbtCellContext::~RbtCellContext() { _RBTOBJECTCOUNTER_DESTR_(_CT); }

RbtCellContext::RbtCellContext(std::ifstream &ifile) {
  int nvbles, nctes;
  ifile >> nvbles >> nctes;
  ninputs = nvbles + nctes;
  std::string name;
  int i = 0, key;
  for (; i < nvbles; i++) {
    ifile >> key >> name;
    vm[key] = new RbtVble(name);
  }
  for (; i < nctes + nvbles; i++) {
    ifile >> key >> name;
    double val = std::atof(name.c_str());
    vm[key] = new RbtVble(name, val);
  }
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
#ifdef _DEBUG
  std::cout << "cell context const\n";
#endif //_DEBUG
}

/*
void RbtCellContext::Clear()
{
    for (RbtIntVbleMapIter it = vm.begin() ; it != vm.end() ; it++)
    {
        if (it->first >= ninputs)
            vm.erase(it);
    }
}*/

double RbtStringContext::Get(RbtModelPtr lig, std::string name) {
  if (name == "LIG_MW")
    return lig->GetTotalMass();
  if (name == "LIG_NATOMS")
    return lig->GetNumAtoms();
  int nAtoms = lig->GetNumAtoms();
  RbtAtomList atomList = lig->GetAtomList();
  RbtBondList bondList = lig->GetBondList();
  if (name == "LIG_NLIPOC")
    return GetNumAtomsWithPredicate(atomList, isAtomLipophilic());
  isHybridState_eq bIsArom(RbtAtom::AROM);
  if (name == "LIG_NAROMATOMS")
    return GetNumAtomsWithPredicate(atomList, bIsArom);
  if (name == "LIG_NHBD")
    return GetNumAtomsWithPredicate(atomList, isAtomHBondDonor());
  if (name == "LIG_NMETAL")
    return GetNumAtomsWithPredicate(atomList, isAtomMetal());
  if (name == "LIG_NGUAN")
    return GetNumAtomsWithPredicate(atomList, isAtomGuanidiniumCarbon());
  if (name == "LIG_NHBA")
    return GetNumAtomsWithPredicate(atomList, isAtomHBondAcceptor());
  if (name == "LIG_NROT") {
    RbtBondList rotatableBondList =
        GetBondListWithPredicate(bondList, isBondRotatable());
    rotatableBondList =
        GetBondListWithPredicate(rotatableBondList, std::not1(isBondToNH3()));
    return static_cast<double>(rotatableBondList.size());
  }
  if (name == "LIG_PERC_LIPOC") {
    unsigned int n = GetNumAtomsWithPredicate(atomList, isAtomLipophilic());
    return (100.0 * n / nAtoms);
  }
  if (name == "LIG_PERC_AROMATOMS") {
    unsigned int n = GetNumAtomsWithPredicate(atomList, bIsArom);
    return (100.0 * n / nAtoms);
  }
  if (name == "LIG_PERC_HBD") {
    unsigned int n = GetNumAtomsWithPredicate(atomList, isAtomHBondDonor());
    return (100.0 * n / nAtoms);
  }
  if (name == "LIG_PERC_METAL") {
    unsigned int n = GetNumAtomsWithPredicate(atomList, isAtomMetal());
    return (100.0 * n / nAtoms);
  }
  if (name == "LIG_PERC_GUAN") {
    unsigned int n =
        GetNumAtomsWithPredicate(atomList, isAtomGuanidiniumCarbon());
    return (100.0 * n / nAtoms);
  }
  if (name == "LIG_PERC_HBA") {
    unsigned int n = GetNumAtomsWithPredicate(atomList, isAtomHBondAcceptor());
    return (100.0 * n / nAtoms);
  }
  if (name == "LIG_NAROMRINGS") {
    RbtAtomListList ringLists = lig->GetRingAtomLists();
    int nAromRings(0); //# aromatic rings
    for (RbtAtomListListIter rIter = ringLists.begin();
         rIter != ringLists.end(); rIter++) {
      if (GetNumAtomsWithPredicate(*rIter, isPiAtom()) == (*rIter).size())
        nAromRings++;
    }
    return nAromRings;
  }
  double posChg(0.0);
  double negChg(0.0);
  for (RbtAtomListConstIter iter = atomList.begin(); iter != atomList.end();
       iter++) {
    double chg = (*iter)->GetGroupCharge();
    if (chg > 0.0)
      posChg += chg;
    else if (chg < 0.0)
      negChg += chg;
  }
  if (name == "LIG_POS_CHG")
    return posChg;
  if (name == "LIG_NEG_CHG")
    return negChg;
  if (name == "LIG_TOT_CHG")
    return (posChg + negChg);
  // FIXME: throw an exception here so this is never reached
  return 0.0;
}

double RbtStringContext::Get(RbtModelPtr spReceptor,
                             RbtDockingSitePtr spDockSite, std::string name) {
  if (name == "SITE_VOL")
    return spDockSite->GetVolume();
  double cavDist = 4.0;
  // Use a fixed definition of cavity atoms -
  // all those within 4A of docking volume
  double neighbR = 4.0; // Sphere radius for counting nearest neighbours
  double threshold = 15;
  // Definition of solvent exposed: neighbours < threshold
  RbtAtomList recepAtomList = spReceptor->GetAtomList();
  RbtAtomList cavAtomList =
      spDockSite->GetAtomList(recepAtomList, 0.0, cavDist);
  RbtAtomList exposedAtomList; // The list of exposed cavity atoms
  for (RbtAtomListConstIter iter = cavAtomList.begin();
       iter != cavAtomList.end(); iter++) {
    unsigned int nNeighb = GetNumAtomsWithPredicate(
        recepAtomList, isAtomInsideSphere((*iter)->GetCoords(), neighbR));
    nNeighb--;
    if (nNeighb < threshold) {
      // std::cout << (*iter)->GetFullAtomName() << "\t" << nNeighb <<
      // std::endl;
      exposedAtomList.push_back(*iter);
    }
  }
  size_t nAtoms = exposedAtomList.size();
  if (name == "SITE_NATOMS")
    return static_cast<double>(nAtoms);
  if (name == "SITE_NLIPOC")
    return GetNumAtomsWithPredicate(exposedAtomList, isAtomLipophilic());
  isHybridState_eq bIsArom(RbtAtom::AROM);
  if (name == "SITE_NAROMATOMS")
    return GetNumAtomsWithPredicate(exposedAtomList, bIsArom);
  if (name == "SITE_NHBD")
    return GetNumAtomsWithPredicate(exposedAtomList, isAtomHBondDonor());
  if (name == "SITE_NMETAL")
    return GetNumAtomsWithPredicate(exposedAtomList, isAtomMetal());
  if (name == "SITE_NGUAN")
    return GetNumAtomsWithPredicate(exposedAtomList, isAtomGuanidiniumCarbon());
  if (name == "SITE_NHBA")
    return GetNumAtomsWithPredicate(exposedAtomList, isAtomHBondAcceptor());
  if (name == "SITE_PERC_LIPOC") {
    unsigned int n =
        GetNumAtomsWithPredicate(exposedAtomList, isAtomLipophilic());
    return (100.0 * n / nAtoms);
  }
  if (name == "SITE_PERC_AROMATOMS") {
    unsigned int n = GetNumAtomsWithPredicate(exposedAtomList, bIsArom);
    return (100.0 * n / nAtoms);
  }
  if (name == "SITE_PERC_HBD") {
    unsigned int n =
        GetNumAtomsWithPredicate(exposedAtomList, isAtomHBondDonor());
    return (100.0 * n / nAtoms);
  }
  if (name == "SITE_PERC_METAL") {
    unsigned int n = GetNumAtomsWithPredicate(exposedAtomList, isAtomMetal());
    return (100.0 * n / nAtoms);
  }
  if (name == "SITE_PERC_GUAN") {
    unsigned int n =
        GetNumAtomsWithPredicate(exposedAtomList, isAtomGuanidiniumCarbon());
    return (100.0 * n / nAtoms);
  }
  if (name == "SITE_PERC_HBA") {
    unsigned int n =
        GetNumAtomsWithPredicate(exposedAtomList, isAtomHBondAcceptor());
    return (100.0 * n / nAtoms);
  }
  double posChg(0.0);
  double negChg(0.0);
  for (RbtAtomListConstIter iter = exposedAtomList.begin();
       iter != exposedAtomList.end(); iter++) {
    double chg = (*iter)->GetGroupCharge();
    if (chg > 0.0)
      posChg += chg;
    else if (chg < 0.0)
      negChg += chg;
  }
  if (name == "SITE_POS_CHG")
    return posChg;
  if (name == "SITE_NEG_CHG")
    return negChg;
  if (name == "SITE_TOT_CHG")
    return (posChg + negChg);
  // FIXME: throw an exception here so this is never reached
  return 0.0;
}

double RbtStringContext::Get(RbtBaseSF *spSF, std::string name,
                             RbtModelPtr lig) {
  RbtStringVariantMap scoreMap;
  spSF->ScoreMap(scoreMap);
  if (scoreMap[name].isEmpty())
    return vm[name]->GetValue(); // lig->GetDataValue(name);
  return scoreMap[name];
}

void RbtStringContext::UpdateLigs(RbtModelPtr lig) {
  for (RbtStringVbleMapIter it = vm.begin(); it != vm.end(); it++) {
    if ((*it).second->IsLig())
      (*it).second->SetValue(Get(lig, (*it).first));
  }
}

void RbtStringContext::UpdateSite(RbtModelPtr spReceptor,
                                  RbtDockingSitePtr spDockSite) {
  for (RbtStringVbleMapIter it = vm.begin(); it != vm.end(); it++) {
    if ((*it).second->IsSite())
      (*it).second->SetValue(Get(spReceptor, spDockSite, (*it).first));
  }
}

void RbtStringContext::UpdateScores(RbtBaseSF *spSF, RbtModelPtr lig) {
  RbtStringVariantMap scoreMap;
  spSF->ScoreMap(scoreMap);
  for (RbtStringVbleMapIter it = vm.begin(); it != vm.end(); it++) {
    if ((*it).second->IsScore()) {
      RbtVariant v = scoreMap[(*it).first];
      if (!v.isEmpty()) {
        (*it).second->SetValue(double(v));
      }
    }
  }
}
