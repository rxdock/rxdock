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

#include "rxdock/Bond.h"

#include <loguru.hpp>

using namespace rxdock;

// Default constructor
// Doesn't make much sense to have a 'default' bond
// but allow it for the time being.
// Default is for atom 0 bonded to itself!
Bond::Bond()
    : m_nBondId(0), m_spAtom1(), m_spAtom2(), m_nFormalBondOrder(1),
      m_dPartialBondOrder(1.0), m_bCyclic(false), m_bSelected(false) {
  _RBTOBJECTCOUNTER_CONSTR_("Bond");
}

// Parameterised constructor
Bond::Bond(int nBondId, AtomPtr &spAtom1, AtomPtr &spAtom2,
           int nFormalBondOrder)
    : m_nBondId(nBondId), m_spAtom1(spAtom1), m_spAtom2(spAtom2),
      m_nFormalBondOrder(nFormalBondOrder),
      m_dPartialBondOrder(nFormalBondOrder), m_bCyclic(false),
      m_bSelected(false) {
  // DM 04 Dec 1998  Register the bond back with the constituent atoms
  bool bOK1 = m_spAtom1->AddBond(this);
  bool bOK2 = m_spAtom2->AddBond(this);
  if (!bOK1) {
    LOG_F(ERROR, "Failed to add bond {} to atom {}", m_nBondId,
          m_spAtom1->GetAtomId());
  }
  if (!bOK2) {
    LOG_F(ERROR, "Failed to add bond {} to atom {}", m_nBondId,
          m_spAtom2->GetAtomId());
  }

  _RBTOBJECTCOUNTER_CONSTR_("Bond");
}

Bond::Bond(json j) {
  j.get_to(*this);
  _RBTOBJECTCOUNTER_CONSTR_("Bond");
}

// Default destructor
Bond::~Bond() {
  // DM 04 Dec 1998  Unregister the bond with the constituent atoms
  bool bOK1 = m_spAtom1->RemoveBond(this);
  bool bOK2 = m_spAtom2->RemoveBond(this);
  if (!bOK1) {
    LOG_F(ERROR, "Failed to remove bond {} from atom {}", m_nBondId,
          m_spAtom1->GetAtomId());
  }
  if (!bOK2) {
    LOG_F(ERROR, "Failed to remove bond {} from atom {}", m_nBondId,
          m_spAtom2->GetAtomId());
  }

  _RBTOBJECTCOUNTER_DESTR_("Bond");
}

// Copy constructor
Bond::Bond(const Bond &bond) {
  m_nBondId = bond.m_nBondId;
  m_spAtom1 = bond.m_spAtom1;
  m_spAtom2 = bond.m_spAtom2;
  m_nFormalBondOrder = bond.m_nFormalBondOrder;
  m_dPartialBondOrder = bond.m_dPartialBondOrder;
  // Copy cyclic flag as the copied bond still points to the same atom
  m_bCyclic = bond.m_bCyclic;
  // But reset selection flag to false
  SetSelectionFlag(false);
  _RBTOBJECTCOUNTER_COPYCONSTR_("Bond");
}

// Copy assignment
Bond &Bond::operator=(const Bond &bond) {
  if (this != &bond) {
    m_nBondId = bond.m_nBondId;
    m_spAtom1 = bond.m_spAtom1;
    m_spAtom2 = bond.m_spAtom2;
    m_nFormalBondOrder = bond.m_nFormalBondOrder;
    m_dPartialBondOrder = bond.m_dPartialBondOrder;
    // Copy cyclic flag as the copied bond still points to the same atom
    m_bCyclic = bond.m_bCyclic;
    // But reset selection flag to false
    SetSelectionFlag(false);
  }
  return *this;
}

///////////////////////////////////////////////
// Stream functions
///////////////////////////////////////////////

// Insertion operator (primarily for debugging)
std::ostream &rxdock::operator<<(std::ostream &s, const Bond &bond) {
  AtomPtr spAtom1 = bond.m_spAtom1;
  AtomPtr spAtom2 = bond.m_spAtom2;

  return s << "Bond ID=" << bond.m_nBondId << ", Atoms " << spAtom1->GetAtomId()
           << "-" << spAtom2->GetAtomId()
           << ", Bond order=" << bond.m_nFormalBondOrder;
  //	   << std::endl << *spAtom1 << std::endl << *spAtom2;
}

void rxdock::to_json(json &j, const Bond &bond) {
  j = json{{"bond-id", bond.m_nBondId},
           {"atom1", *bond.m_spAtom1},
           {"atom2", *bond.m_spAtom2},
           {"formal-bond-order", bond.m_nFormalBondOrder},
           {"partial-bond-order", bond.m_dPartialBondOrder},
           {"cyclic", bond.m_bCyclic},
           {"selected", bond.m_bSelected}};
}

void rxdock::from_json(const json &j, Bond &bond) {
  j.at("bond-id").get_to(bond.m_nBondId);
  bond.m_spAtom1 = new Atom(j.at("atom1"));
  bond.m_spAtom2 = new Atom(j.at("atom2"));
  j.at("formal-bond-order").get_to(bond.m_nFormalBondOrder);
  j.at("partial-bond-order").get_to(bond.m_dPartialBondOrder);
  j.at("cyclic").get_to(bond.m_bCyclic);
  j.at("selected").get_to(bond.m_bSelected);
}

////////////////////////////////////////
// Public accessor functions
///////////////////////////

////////////////////////////////////////
// Public methods
////////////////

// Returns bond length
double Bond::Length() const { return BondLength(m_spAtom1, m_spAtom2); }

//////////////////////////////////////////
// Non-member functions (in rxdock namespace)
//////////////////////////////////////////

// Is bond rotatable ?
bool rxdock::isBondRotatable::operator()(Bond *pBond) const {
  // Useful predicates
  isFFType_eq bIsC_SP2("C_SP2");
  isFFType_eq bIsC_SP("C_SP");
  isFFType_eq bIsN_TRI("N_TRI");
  isFFType_eq bIsN_SP2P("N_SP2+");
  isFFType_eq bIsM_SP("N_SP");

  // Eliminate cyclic bonds
  if (pBond->GetCyclicFlag())
    return false;

  // Eliminate non-single bonds
  if (pBond->GetFormalBondOrder() > 1)
    return false;

  // Allow single bonds between two pendant ring systems
  AtomPtr spAtom1(pBond->GetAtom1Ptr());
  AtomPtr spAtom2(pBond->GetAtom2Ptr());
  if (spAtom1->GetCyclicFlag() && spAtom2->GetCyclicFlag())
    return true;

  // Eliminate Csp2-Ntri bonds (e.g. amides, guanidines)
  // and Csp2-Nsp2+ bonds (e.g. aromatic nitro groups) as these are partially
  // conjugated if ( (bIsC_SP2(spAtom1) && (bIsN_TRI(spAtom2) ||
  // bIsN_SP2P(spAtom2))) ||
  //     (bIsC_SP2(spAtom2) && (bIsN_TRI(spAtom1) || bIsN_SP2P(spAtom1))) )
  //  return false;
  // DM 21 Nov 2002 - consider these cases separately in a bit more detail
  // Amides, where the C_SP2 atom is specifically NOT in a ring
  if ((bIsC_SP2(spAtom1) && !spAtom1->GetCyclicFlag() && bIsN_TRI(spAtom2)) ||
      (bIsC_SP2(spAtom2) && !spAtom2->GetCyclicFlag() && bIsN_TRI(spAtom1)))
    return false;
  // Nitro groups
  if ((bIsC_SP2(spAtom1) && bIsN_SP2P(spAtom2)) ||
      (bIsC_SP2(spAtom2) && bIsN_SP2P(spAtom1)))
    return false;
  // Terminal Csp2-NH2 groups where the N is trigonal (mainly anilines)
  if ((bIsC_SP2(spAtom1) && bIsN_TRI(spAtom2) &&
       (spAtom2->GetCoordinationNumber(1) == 2)) ||
      (bIsC_SP2(spAtom2) && bIsN_TRI(spAtom1) &&
       (spAtom1->GetCoordinationNumber(1) == 2)))
    return false;

  // Eliminate terminal atoms (making only one bond)
  if ((spAtom1->GetNumBonds() == 1) || (spAtom2->GetNumBonds() == 1))
    return false;

  // Eliminate bonds by bond type (exclude any bond without at least one atom
  // which is SP3) i.e. we count SP3-SP3, SP3-SP2, SP3-AROM, SP3-TRI as
  // rotatable Atom::eHybridState eState1(spAtom1->GetHybridState());
  // Atom::eHybridState eState2(spAtom2->GetHybridState());
  // if (!( (eState1 == Atom::SP3) || (eState2 == Atom::SP3) ))
  //  return false;
  //
  // DM 05 May 1999 - also allow bonds to O_TRI and S_TRI
  // if ( (spAtom1->GetHybridState() != Atom::SP3) &&
  //     (spAtom2->GetHybridState() != Atom::SP3) &&
  //     !bIsO_TRI(spAtom1) && !bIsS_TRI(spAtom1) &&
  //     !bIsO_TRI(spAtom2) && !bIsS_TRI(spAtom2) )
  //  return false;

  // Eliminate terminal CH3 and NH3 groups (i.e. anything bonded to 3 explicit
  // or implicit hydrogens) DM 8 Feb 2000 - only eliminate CH3
  if (((spAtom1->GetAtomicNo() == 6) &&
       (spAtom1->GetNumImplicitHydrogens() +
            spAtom1->GetCoordinationNumber(1) ==
        3)) ||
      ((spAtom2->GetAtomicNo() == 6) &&
       (spAtom2->GetNumImplicitHydrogens() +
            spAtom2->GetCoordinationNumber(1) ==
        3)))
    return false;

  // DM 21 Nov 2002 - eliminate -C#N (-C triple-bond N)
  if ((bIsC_SP(spAtom1) && (spAtom1->GetCoordinationNumber("N_SP") > 0)) ||
      (bIsC_SP(spAtom2) && (spAtom2->GetCoordinationNumber("N_SP") > 0)))
    return false;

  // If we get this far, the bond is rotatable
  return true;
}

// DM 24 Sep 2001
// Is bond to a terminal NH3+ group?
// Use to filter rotable bond list for scoring function purposes
bool rxdock::isBondToNH3::operator()(Bond *pBond) const {
  AtomPtr spAtom1(pBond->GetAtom1Ptr());
  AtomPtr spAtom2(pBond->GetAtom2Ptr());
  if (((spAtom1->GetAtomicNo() == 7) &&
       (spAtom1->GetNumImplicitHydrogens() +
            spAtom1->GetCoordinationNumber(1) ==
        3)) ||
      ((spAtom2->GetAtomicNo() == 7) &&
       (spAtom2->GetNumImplicitHydrogens() +
            spAtom2->GetCoordinationNumber(1) ==
        3)))
    return true;
  else
    return false;
}

// Is bond to a terminal OH group?
bool rxdock::isBondToOH::operator()(Bond *pBond) const {
  AtomPtr spAtom1(pBond->GetAtom1Ptr());
  AtomPtr spAtom2(pBond->GetAtom2Ptr());
  if (((spAtom1->GetAtomicNo() == 8) &&
       (spAtom1->GetNumImplicitHydrogens() +
            spAtom1->GetCoordinationNumber(1) ==
        1)) ||
      ((spAtom2->GetAtomicNo() == 8) &&
       (spAtom2->GetNumImplicitHydrogens() +
            spAtom2->GetCoordinationNumber(1) ==
        1)))
    return true;
  else
    return false;
}

// DM 7 June 1999
// Is bond an amide bond?
bool rxdock::isBondAmide::operator()(Bond *pBond) const {
  // Useful predicates
  // isFFType_eq bIsC_SP2("C_SP2");
  // isFFType_eq bIsN_TRI("N_TRI");
  isAtomicNo_eq isC(6);
  isAtomicNo_eq isN(7);
  isAtomicNo_eq isO(8);
  isHybridState_eq isSP2(Atom::SP2);
  isHybridState_eq isTRI(Atom::TRI);

  // 07 Nov 2003 (DM) - return false for all cyclic bonds
  // even if they are "amide-like" (e.g. in thymine)
  if (pBond->GetCyclicFlag())
    return false;

  AtomPtr spAtom1(pBond->GetAtom1Ptr());
  AtomPtr spAtom2(pBond->GetAtom2Ptr());

  if (isC(spAtom1) && isSP2(spAtom1) && isN(spAtom2) && isTRI(spAtom2)) {
    AtomList oxygens =
        GetAtomListWithPredicate(GetBondedAtomList(spAtom1), isO);
    int nOSP2 = GetNumAtomsWithPredicate(oxygens, isSP2);
    LOG_F(1, "Amide check on {} - {} #OSP2 = {}", spAtom1->GetFullAtomName(),
          spAtom2->GetFullAtomName(), nOSP2);
    return (nOSP2 == 1);
  } else if (isC(spAtom2) && isSP2(spAtom2) && isN(spAtom1) && isTRI(spAtom1)) {
    AtomList oxygens =
        GetAtomListWithPredicate(GetBondedAtomList(spAtom2), isO);
    int nOSP2 = GetNumAtomsWithPredicate(oxygens, isSP2);
    LOG_F(1, "Amide check on {} - {} #OSP2 = {}", spAtom2->GetFullAtomName(),
          spAtom1->GetFullAtomName(), nOSP2);
    return (nOSP2 == 1);
  } else
    return false;
}

////////////////////////////////////////////
// Bond list functions (implemented as STL algorithms)
////////////////////////////////////////////

// Selected bonds
void rxdock::SetBondSelectionFlagsInList(BondList &bondList, bool bSelected) {
  for (BondListIter iter = bondList.begin(); iter != bondList.end(); iter++)
    (*iter)->SetSelectionFlag(bSelected);
}

// Cyclic bonds
void rxdock::SetBondCyclicFlagsInList(BondList &bondList, bool bCyclic) {
  for (BondListIter iter = bondList.begin(); iter != bondList.end(); iter++)
    (*iter)->SetCyclicFlag(bCyclic);
}
