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

#include "rxdock/SetupPMFSF.h"

#include <loguru.hpp>

#include <functional>

using namespace rxdock;

std::string SetupPMFSF::_CT("SetupPMFSF");

SetupPMFSF::SetupPMFSF(const std::string &strName) : BaseSF(_CT, strName) {
  LOG_F(2, "SetupPMFSF parameterised constructor");
  Disable();
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

SetupPMFSF::~SetupPMFSF() {
  LOG_F(2, "SetupPMFSF destructor");
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

void SetupPMFSF::SetupReceptor() {
  theReceptorList = GetAtomListWithPredicate(GetReceptor()->GetAtomList(),
                                             std::not1(isAtomicNo_eq(1)));
  SetupReceptorPMFTypes();
  for (long i = 0; i < theReceptorList.size(); i++) {
    LOG_F(1, "SetupPMFSF::SetupReceptor: {} {} type {} No {}", i,
          theReceptorList[i]->GetFullAtomName(),
          PMFType2Str(theReceptorList[i]->GetPMFType()),
          theReceptorList[i]->GetAtomicNo());
  }
}

void SetupPMFSF::SetupScore() {}

double SetupPMFSF::RawScore() const { return 0.0; }

void SetupPMFSF::SetupLigand() {
  theLigandList.clear();
  if (GetLigand().Null()) {
    LOG_F(WARNING, "SetupPMFSF::SetupLigand: Ligand is not defined");
    return;
  } else {
    // theLigandList = GetLigand()->GetAtomList();
    theLigandList = GetAtomListWithPredicate(GetLigand()->GetAtomList(),
                                             std::not1(isAtomicNo_eq(1)));
    SetupLigandPMFTypes();
    for (int i = 0; i < theLigandList.size(); i++) {
      LOG_F(1, "SetupPMFSF:SetupLigand: {} {} type {} No {}", i,
            theLigandList[i]->GetFullAtomName(),
            PMFType2Str(theLigandList[i]->GetPMFType()),
            theLigandList[i]->GetAtomicNo());
    }
  }
}

PMFType SetupPMFSF::GetPMFfor_rC(AtomPtr anAtom) {
  LOG_F(2, "SetupPMFSF::GetPMFfor_rC: {} Hybrid: {}", anAtom->GetFullAtomName(),
        ConvertHybridStateToString(anAtom->GetHybridState()));

  isAtomicNo_eq isO(8);
  isAtomicNo_eq isN(7);
  isAtomicNo_eq isS(16);

  AtomList obdAtomList =
      GetAtomListWithPredicate(GetBondedAtomList(anAtom), isO);
  AtomList nbdAtomList =
      GetAtomListWithPredicate(GetBondedAtomList(anAtom), isN);
  AtomList sbdAtomList =
      GetAtomListWithPredicate(GetBondedAtomList(anAtom), isS);
  int nBoundedHetero =
      obdAtomList.size() + nbdAtomList.size() + sbdAtomList.size();

  if (Atom::AROM == anAtom->GetHybridState()) { // if aromatic
    if (nBoundedHetero > 0)                     // polar
      return cP;
    else
      return cF;
  }
  // non-aromatic
  // check for charged [CN] neighbours
  AtomList bList = GetBondedAtomList(anAtom);
  for (AtomListIter bIter = bList.begin(); bIter != bList.end(); bIter++) {
    if (isO((*bIter)) &&
        (*bIter)->GetGroupCharge() < 0.0) // negatively charged oxigen neigbour
      return CO;                          // cee-oh
    if (isN((*bIter)) &&
        IsChargedNitrogen(
            (*bIter))) // check for guanadine or any other charge on N
      return CN;
  }
  // rest is non-aromatic and has no charged [ONS] neighbours
  // check for polarity
  if (nBoundedHetero > 0) // has at least 1 [ONS] neighbour
    return CP;            // polar
  else
    return CF; //
}

PMFType SetupPMFSF::GetPMFfor_rO(AtomPtr anAtom) {
  LOG_F(2, "SetupPMFSF::GetPMFfor_rO: {}", anAtom->GetFullAtomName());
  // check charge
  if (anAtom->GetGroupCharge() < 0.0)
    return OC;
  // check wether is it water
  isAtomicNo_eq bIsH(1);
  AtomList hbdAtomList =
      GetAtomListWithPredicate(GetBondedAtomList(anAtom), bIsH);
  if (2 == hbdAtomList.size()) // must be water
    return OW;
  else if (1 == hbdAtomList.size()) // if not water assume it is a H-bond donor
                                    // in an -OH group
    return OD;
  else
    return OA; // still no match? must be a backbone O or ASN, GLN O
}

PMFType SetupPMFSF::GetPMFfor_rN(AtomPtr anAtom) {
  LOG_F(2, "SetupPMFSF::GetPMFfor_rN: {} Hybrid: {}", anAtom->GetFullAtomName(),
        ConvertHybridStateToString(anAtom->GetHybridState()));
  // check for charge
  if (IsChargedNitrogen(anAtom))
    return NC;
  // check for donors
  isAtomHBondDonor isHBondDonor;
  AtomList donorList =
      GetAtomListWithPredicate(GetBondedAtomList(anAtom), isHBondDonor);
  if (donorList.size() > 0)
    return ND;
  else
    return NR;
}

// sulphur in MET (SA) and CYS (SD)
// it is not likely to find them elsewhere than in proteins
PMFType SetupPMFSF::GetPMFfor_rS(AtomPtr anAtom) {
  LOG_F(2, "SetupPMFSF::GetPMFfor_rS: {} Hybrid: {}", anAtom->GetFullAtomName(),
        ConvertHybridStateToString(anAtom->GetHybridState()));
  if (std::string::npos != anAtom->GetFullAtomName().find("MET"))
    return SA; // metionin H-bond acceptor
  else
    return SD; // cystein H-bond donor
}

void SetupPMFSF::SetupReceptorPMFTypes(void) {
  LOG_F(1, "SetupPMFSF receptor size {}", theReceptorList.size());
  for (unsigned long i = 0; i < theReceptorList.size(); i++) {
    switch (theReceptorList[i]->GetAtomicNo()) {
    case 1:
      theReceptorList[i]->SetPMFType(HH);
      break;
    case 6: // carbon
      theReceptorList[i]->SetPMFType(GetPMFfor_rC(theReceptorList[i]));
      break;
    case 7: // nitrogen
      theReceptorList[i]->SetPMFType(GetPMFfor_rN(theReceptorList[i]));
      break;
    case 8: // oxigen
      theReceptorList[i]->SetPMFType(GetPMFfor_rO(theReceptorList[i]));
      break;
    case 15:
      theReceptorList[i]->SetPMFType(P);
      break;
    case 16: // sulphur
      theReceptorList[i]->SetPMFType(GetPMFfor_rS(theReceptorList[i]));
      break;
    // these heteroatoms belongs to the receptor
    case 12:
    case 20:
      theReceptorList[i]->SetPMFType(Mg);
      break; // treat Ca as Mg
    default:
      // throw an exception but not bothered about it now
      // throw ();
      break;
    }
    LOG_F(2, "SetupPMFSF receptor type: {} for: {}",
          PMFType2Str(theReceptorList[i]->GetPMFType()),
          theReceptorList[i]->GetAtomicNo());
  }
}

bool SetupPMFSF::IsChargedNitrogen(AtomPtr anAtom) {
  LOG_F(2, "SetupPMFSF::IsChargedNitrogen");
  LOG_F(1, "amino acid {}", anAtom->GetSubunitName());
  isAtomicNo_eq bIsN(7);
  if (!bIsN(anAtom)) // is it N at all?
    return false;
  if (anAtom->GetGroupCharge() > 0.0) // some charge on the N itself
    return true;
  else {
    isAtomicNo_eq bIsH(1);
    AtomList hbdAtomList =
        GetAtomListWithPredicate(GetBondedAtomList(anAtom), bIsH);
    for (AtomListConstIter hIter = hbdAtomList.begin();
         hIter != hbdAtomList.end(); hIter++) {
      if ((*hIter)->GetGroupCharge() >
          0.0) // if there is charge on H -> guanadine
        return true;
    }
  }
  // othervise no charge
  return false;
}

// Set PMF type for carbons.
PMFType SetupPMFSF::GetPMFfor_lC(AtomPtr anAtom) {
  LOG_F(2, "SetupPMFSF::GetPMFfor_lC: {} Hybrid: {}", anAtom->GetFullAtomName(),
        ConvertHybridStateToString(anAtom->GetHybridState()));
  if (Atom::SP == anAtom->GetHybridState()) // sp has only one PMF type
    return C0;                              // C-zero and not cee-oh

  isAtomicNo_eq isH(1);
  isAtomicNo_eq isC(6);

  AtomList hbdAtomList =
      GetAtomListWithPredicate(GetBondedAtomList(anAtom), isH);
  AtomList cbdAtomList =
      GetAtomListWithPredicate(GetBondedAtomList(anAtom), isC);
  // bonded list will be needed anyway
  AtomList bList = GetBondedAtomList(anAtom);
  if (bList.size() == hbdAtomList.size() + cbdAtomList.size()) { // if non-polar
    switch (anAtom->GetHybridState()) {
    case Atom::SP2:
      return C3;
      break;
    case Atom::SP3:
      return CF;
      break;
    case Atom::AROM:
      return cF;
      break;
    default:
      return PMF_UNDEFINED;
      break;
    }
  } else { // considered as polar
    // there is a chance to bound to a charged O or N
    isAtomicNo_eq isO(8);
    isAtomicNo_eq isN(7);

    AtomList bList = GetBondedAtomList(anAtom);
    for (AtomListIter bIter = bList.begin(); bIter != bList.end(); bIter++) {
      if (isO((*bIter)) && (*bIter)->GetGroupCharge() <
                               0.0) // negatively charged oxigen neigbour
        return CO;                  // cee-oh
      else if (IsChargedNitrogen(anAtom))
        return CN;
    }
    // if no charged atom, check hybridisation
    switch (anAtom->GetHybridState()) {
    case Atom::SP2:
      return CW;
      break;
    case Atom::SP3:
      return CP;
      break;
    case Atom::AROM:
      return cP;
      break;
    default:
      return PMF_UNDEFINED;
      break;
    }
  }
}

// same for nitrogen
PMFType SetupPMFSF::GetPMFfor_lN(AtomPtr anAtom) {
  Atom::eHybridState theHybState = anAtom->GetHybridState();
  LOG_F(2, "SetupPMFSF::GetPMFfor_lN: {} Hybrid: {}", anAtom->GetFullAtomName(),
        ConvertHybridStateToString(theHybState));
  // sp has only one PMF type
  if (Atom::SP == theHybState)
    return N0; // N-zero and not en-oh
  // checking are there neighbours other than [CH]
  // the lists we get will be useful later as well
  // NOTE: excluding N as well since most of the ligand Ns will be typed
  // as NS instead of ND/NA/NR
  isAtomicNo_eq bIsH(1);
  isAtomicNo_eq bIsC(6);
  isAtomicNo_eq bIsN(7);
  AtomList hbdAtomList =
      GetAtomListWithPredicate(GetBondedAtomList(anAtom), bIsH);
  AtomList cbdAtomList =
      GetAtomListWithPredicate(GetBondedAtomList(anAtom), bIsC);
  AtomList nbdAtomList =
      GetAtomListWithPredicate(GetBondedAtomList(anAtom), bIsN);
  // if num_of_C + num_of_H < all_of_bounded than type is NS (there is something
  // else than C or H)
  if ((hbdAtomList.size() + cbdAtomList.size() + nbdAtomList.size() <
       GetBondedAtomList(anAtom).size()) &&
      Atom::AROM != theHybState)
    return NS;
  // if in aromatic (planar) ring:
  if (Atom::AROM == theHybState)
    return NR;
  // charged
  else if (IsChargedNitrogen(anAtom))
    return NC;
  // if planar but not in aromatic ring (sp2 or tri)
  // not bound to H but to 2 or 3 C
  if ((Atom::SP2 == theHybState || Atom::TRI == theHybState) &&
      (0 == hbdAtomList.size() && cbdAtomList.size() > 1))
    return NP;

  // if H-donor/acceptor outside of an aromatic ring
  isAtomHBondDonor isHBondDonor;
  isAtomHBondAcceptor isHBondAcceptor;

  AtomList donorList =
      GetAtomListWithPredicate(GetBondedAtomList(anAtom), isHBondDonor);
  if (donorList.size() > 0)
    return ND;
  else if (isHBondAcceptor(anAtom))
    return NA;
  // in ring like in FMN but not necessary aromatic and can be even SP3
  isAtomCyclic isCyclic;
  if (isCyclic(anAtom))
    return NR;

  // should not be reached
  return PMF_UNDEFINED;
}

PMFType SetupPMFSF::GetPMFfor_lO(AtomPtr anAtom) {
  Atom::eHybridState theHybState = anAtom->GetHybridState();
  LOG_F(2, "SetupPMFSF::GetPMFfor_lO: {} Hybrid: {}", anAtom->GetFullAtomName(),
        ConvertHybridStateToString(theHybState));
  isAtomHBondDonor isHBondDonor;
  isAtomHBondAcceptor isHBondAcceptor;
  isAtomCyclic isCyclic;
  isAtomicNo_eq bIsC(6);

  // first check negative charge
  if (anAtom->GetGroupCharge() < 0.0)
    return OC;
  // check for planar ring
  else if (Atom::AROM == theHybState)
    return OR;
  // first we have to check ether bonds since isAtomHBondAcceptor includes
  // ethers
  AtomList cbdAtomList =
      GetAtomListWithPredicate(GetBondedAtomList(anAtom), bIsC);
  if (cbdAtomList.size() > 1)
    return OE;
  // check for H-bond donor role
  AtomList donorList =
      GetAtomListWithPredicate(GetBondedAtomList(anAtom), isHBondDonor);
  if (donorList.size() > 0)
    return OD;
  // H-bond acceptor role
  else if (isHBondAcceptor(anAtom))
    return OA;
  // remainder
  return OS;
}

// PMF type for sulphur
PMFType SetupPMFSF::GetPMFfor_lS(AtomPtr anAtom) {
  LOG_F(2, "SetupPMFSF::GetPMFfor_lS: {} Hybrid: {}", anAtom->GetFullAtomName(),
        ConvertHybridStateToString(anAtom->GetHybridState()));
  isAtomHBondDonor isHBondDonor;
  AtomList donorList =
      GetAtomListWithPredicate(GetBondedAtomList(anAtom), isHBondDonor);
  if (donorList.size() > 0)
    return SD;
  else
    return SA;
}

void SetupPMFSF::SetupLigandPMFTypes(void) {
  for (unsigned int i = 0; i < theLigandList.size(); i++) {
    switch (theLigandList[i]->GetAtomicNo()) {
    case 1:
      theLigandList[i]->SetPMFType(HL);
      break;
    case 6: // carbon
      theLigandList[i]->SetPMFType(GetPMFfor_lC(theLigandList[i]));
      break;
    case 7: // nitrogen
      theLigandList[i]->SetPMFType(GetPMFfor_lN(theLigandList[i]));
      break;
    case 8: // oxigen
      theLigandList[i]->SetPMFType(GetPMFfor_lO(theLigandList[i]));
      break;
    case 9:
      theLigandList[i]->SetPMFType(F);
      break;
    case 12:
    case 20:
      theLigandList[i]->SetPMFType(Mg);
      break; // treat Ca as Mg
    case 15:
      theLigandList[i]->SetPMFType(P);
      break;
    case 16: // sulphur
      theLigandList[i]->SetPMFType(GetPMFfor_lS(theLigandList[i]));
      break;
    case 17:
      theLigandList[i]->SetPMFType(CL);
      break;
    case 23:
      theLigandList[i]->SetPMFType(V);
      break;
    case 25:
      theLigandList[i]->SetPMFType(Mn);
      break;
    case 26:
      theLigandList[i]->SetPMFType(Fe);
      break;
    case 30:
      theLigandList[i]->SetPMFType(Zn);
      break;
    case 35:
      theLigandList[i]->SetPMFType(Br);
      break;
    default:
      // throw an exception but not bothered about it now
      // throw ();
      break;
    }
  }
}
