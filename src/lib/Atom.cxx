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

#include "Atom.h"
#include "Bond.h"
#include "Model.h"
#include "PseudoAtom.h"

#include <loguru.hpp>

#include <functional>

using namespace rxdock;

///////////////////////////////////////////////
// Constructors / destructors
///////////////////////////////////////////////

// Default constructor
// Initialise 2-D attributes to suitable defaults
// Initialise 3-D attributes to zero/null
Atom::Atom()
    : m_nAtomicNo(6), m_nAtomId(0), m_strAtomName("C"), m_strSubunitId("1"),
      m_strSubunitName("RES"), m_strSegmentName("SEG1"),
      m_eState(UNDEFINED), // DM 8 Dec 1998 Changed from SP3 to UNDEFINED
      m_nHydrogens(0), m_nFormalCharge(0), m_pModel(nullptr), m_bCyclic(false),
      m_bSelected(false), m_bUser1(false),
      m_dUser1(1.0), // DM 27 Jul 2000 - initialise user values to 1 as they are
                     // commonly used as weightings
      m_dUser2(1.0), m_nPMFType(PMF_UNDEFINED),
      m_triposType(TriposAtomType::UNDEFINED), m_coord(0.0, 0.0, 0.0),
      m_dPartialCharge(0.0), m_dGroupCharge(0.0), m_dAtomicMass(0.0),
      m_dVdwRadius(0.0), m_strFFType("") {
  _RBTOBJECTCOUNTER_CONSTR_("Atom");
}

// Constructor supplying all 2-D parameters
// Initialise 3-D attributes to zero/null
Atom::Atom(int nAtomId, int nAtomicNo /*= 6*/,
           std::string strAtomName /*= "C"*/,
           std::string strSubunitId /*= "1"*/,
           std::string strSubunitName /*= "RES"*/,
           std::string strSegmentName /*= "SEG1"*/,
           eHybridState eState /*= UNDEFINED*/, unsigned int nHydrogens /*= 0*/,
           int nFormalCharge /*= 0.0*/
           )
    : m_nAtomicNo(nAtomicNo), m_nAtomId(nAtomId), m_strAtomName(strAtomName),
      m_strSubunitId(strSubunitId), m_strSubunitName(strSubunitName),
      m_strSegmentName(strSegmentName), m_eState(eState),
      m_nHydrogens(nHydrogens), m_nFormalCharge(nFormalCharge),
      m_pModel(nullptr), m_bCyclic(false), m_bSelected(false), m_bUser1(false),
      m_dUser1(1.0), // DM 27 Jul 2000 - initialise user values to 1 as they are
                     // commonly used as weightings
      m_dUser2(1.0), m_nPMFType(PMF_UNDEFINED),
      m_triposType(TriposAtomType::UNDEFINED), m_coord(0.0, 0.0, 0.0),
      m_dPartialCharge(0.0), m_dGroupCharge(0.0), m_dAtomicMass(0.0),
      m_dVdwRadius(0.0), m_strFFType("") {
  _RBTOBJECTCOUNTER_CONSTR_("Atom");
}

// Default destructor
Atom::~Atom() { _RBTOBJECTCOUNTER_DESTR_("Atom"); }

// Copy constructor
Atom::Atom(const Atom &atom) {
  m_nAtomId = atom.m_nAtomId;
  m_nAtomicNo = atom.m_nAtomicNo;
  m_strAtomName = atom.m_strAtomName;
  m_strSubunitId = atom.m_strSubunitId;
  m_strSubunitName = atom.m_strSubunitName;
  m_strSegmentName = atom.m_strSegmentName;
  m_eState = atom.m_eState;
  m_nHydrogens = atom.m_nHydrogens;
  m_nFormalCharge = atom.m_nFormalCharge;
  m_coord = atom.m_coord;
  m_dPartialCharge = atom.m_dPartialCharge;
  m_dGroupCharge = atom.m_dGroupCharge;
  m_dAtomicMass = atom.m_dAtomicMass;
  m_dVdwRadius = atom.m_dVdwRadius;
  m_strFFType = atom.m_strFFType;
  m_savedCoords = atom.m_savedCoords;
  m_bUser1 = atom.m_bUser1;
  m_dUser1 = atom.m_dUser1;
  m_dUser2 = atom.m_dUser2;
  m_nPMFType = atom.m_nPMFType;
  m_triposType = atom.m_triposType;
  // Copied atoms no longer belong to the model so set to nullptr here
  SetModelPtr(nullptr);
  // Copied atoms no longer belong to the bonds so erase the bond map
  ClearBondMap();
  // Set the cyclic flag to false as the atom isn't bonded to anything
  SetCyclicFlag(false);
  // Reset the selection flag also
  SetSelectionFlag(false);
  _RBTOBJECTCOUNTER_COPYCONSTR_("Atom");
}

// Copy assignment
Atom &Atom::operator=(const Atom &atom) {
  if (this != &atom) { // beware of self-assignment
    m_nAtomId = atom.m_nAtomId;
    m_nAtomicNo = atom.m_nAtomicNo;
    m_strAtomName = atom.m_strAtomName;
    m_strSubunitId = atom.m_strSubunitId;
    m_strSubunitName = atom.m_strSubunitName;
    m_strSegmentName = atom.m_strSegmentName;
    m_eState = atom.m_eState;
    m_nHydrogens = atom.m_nHydrogens;
    m_nFormalCharge = atom.m_nFormalCharge;
    m_coord = atom.m_coord;
    m_dPartialCharge = atom.m_dPartialCharge;
    m_dGroupCharge = atom.m_dGroupCharge;
    m_dAtomicMass = atom.m_dAtomicMass;
    m_dVdwRadius = atom.m_dVdwRadius;
    m_strFFType = atom.m_strFFType;
    m_savedCoords = atom.m_savedCoords;
    m_bUser1 = atom.m_bUser1;
    m_dUser1 = atom.m_dUser1;
    m_dUser2 = atom.m_dUser2;
    m_nPMFType = atom.m_nPMFType;
    m_triposType = atom.m_triposType;
    // Copied atoms no longer belong to the model so set to nullptr here
    SetModelPtr(nullptr);
    // Copied atoms no longer belong to the bonds so erase the bond map
    ClearBondMap();
    // Set the cyclic flag to false as the atom isn't bonded to anything
    SetCyclicFlag(false);
    // Reset the selection flag also
    SetSelectionFlag(false);
  }
  return *this;
}

///////////////////////////////////////////////
// Stream functions
///////////////////////////////////////////////

// Insertion operator (primarily for debugging)
std::ostream &rxdock::operator<<(std::ostream &s, const Atom &atom) {
  return atom.Print(s);
}

// Virtual function for dumping atom details to an output stream
// Derived classes (e.g. pseudoatom) can override if required
std::ostream &Atom::Print(std::ostream &s) const {
  // Get owning model name (if any)
  std::string strModelName = "Orphan";
  if (m_pModel != nullptr)
    strModelName = m_pModel->GetName();

  TriposAtomType triposType;

  return s << "(" << strModelName << ")"
           << ", ID=" << m_nAtomId << ", AtNo=" << m_nAtomicNo << ", "
           << m_strSegmentName << "," << m_strSubunitName << ","
           << m_strSubunitId << "," << m_strAtomName
           << ", Hybrid=" << ConvertHybridStateToString(m_eState)
           << ", #Hyd=" << m_nHydrogens << ", GrpChg=" << m_dGroupCharge
           << ", Tripos="
           << triposType.Type2Str(GetTriposType())
           //	   << ", FFType=" << m_strFFType
           << ", Cyclic=" << m_bCyclic;
}

///////////////////////////////////////////////
// Public accessor functions
// 2-D attributes
///////////////////////////////////////////////

bool Atom::GetEnabled() const {
  return (m_pModel) ? m_pModel->GetEnabled() : true;
}

// Returns composite of segment name, subunit id and name, and atom name
std::string Atom::GetFullAtomName() const {
  return m_strSegmentName + ":" + m_strSubunitName + "_" + m_strSubunitId +
         ":" + m_strAtomName;
}

// DM 04 Dec 1998  Add functions to handle bond map
// Add a bond to the bond map - returns true if OK, false if this atom is not a
// member of the bond or if bond already added
bool Atom::AddBond(Bond *pBond) {
  // Case 1 - this atom is the first atom in the bond
  //(key=pBond, value = true) as this atom is atom 1
  if (pBond->GetAtom1Ptr().Ptr() == this) {
    std::pair<BondMapIter, bool> p =
        m_bondMap.insert(std::make_pair(pBond, true));
    return p.second; // p.second is true if bond inserted, false if bond was
                     // already present
  }
  // Case 2 - this atom is the second atom in the bond
  //(key=pBond, value = false) as this atom is atom 2
  else if (pBond->GetAtom2Ptr().Ptr() == this) {
    std::pair<BondMapIter, bool> p =
        m_bondMap.insert(std::make_pair(pBond, false));
    return p.second; // p.second is true if bond inserted, false if bond was
                     // already present
  } else
    return false;
}

// Remove a bond from the bond map - returns true if OK, false if bond not
// present in map
bool Atom::RemoveBond(Bond *pBond) {
  if (m_bondMap.erase(pBond) > 0) // erase returns the number of element deleted
    return true;
  else
    return false;
}

// Returns number of cyclic bonds in map
unsigned int Atom::GetNumCyclicBonds() const {
  unsigned int nCount(0);
  for (BondMapConstIter iter = m_bondMap.begin(); iter != m_bondMap.end();
       iter++)
    if ((*iter).first->GetCyclicFlag())
      nCount++;
  return nCount;
}

// Returns bond map with only cyclic bonds included
BondMap Atom::GetCyclicBondMap() const {
  BondMap cyclicBondMap;
  for (BondMapConstIter iter = m_bondMap.begin(); iter != m_bondMap.end();
       iter++)
    if ((*iter).first->GetCyclicFlag())
      cyclicBondMap.insert(*iter);
  return cyclicBondMap;
}

// Returns total formal bond order around atom
int Atom::GetTotalFormalBondOrder() const {
  int nTot(0);
  for (BondMapConstIter iter = m_bondMap.begin(); iter != m_bondMap.end();
       iter++)
    nTot += (*iter).first->GetFormalBondOrder();
  return nTot;
}

// Returns max formal bond order for all bonds to atom
int Atom::GetMaxFormalBondOrder() const {
  int nMax(0);
  for (BondMapConstIter iter = m_bondMap.begin(); iter != m_bondMap.end();
       iter++)
    nMax = std::max(nMax, (*iter).first->GetFormalBondOrder());
  return nMax;
}

// Clears the bond map - should only need to be called by the copy constructors,
// hence private
void Atom::ClearBondMap() { m_bondMap.clear(); }

///////////////////////////////////////////////
// Other public methods
///////////////////////////////////////////////

// DM 08 Feb 1999 - all saved coords are now saved in a
// std::map<UInt,Coord> map key=0 is reserved for the default SaveCoords
// and RevertCoords
void Atom::SaveCoords(unsigned int coordNum) {
  m_savedCoords[coordNum] = m_coord;
}

void Atom::RevertCoords(unsigned int coordNum) {
  UIntCoordMapConstIter iter = m_savedCoords.find(coordNum);
  if (iter != m_savedCoords.end()) {
    m_coord = (*iter).second;
  } else
    throw InvalidRequest(_WHERE_,
                         "RevertCoords failed on atom " + GetFullAtomName());
}

// DM 04 Dec 1998  Now we have the bond map, we can easily provide coordination
// numbers This version returns the total number of coordinated atoms (includes
// implicit hydrogens) Equivalent to GetNumBonds()+GetNumImplicitHydrogens()
unsigned int Atom::GetCoordinationNumber() const {
  return m_bondMap.size() + m_nHydrogens;
}

// This version returns the number of coordinated atoms of a given element
// Note: for hydrogens, the implicit atoms are excluded
// so GetCoordinationNumber(1) will return the number of explicit hydrogens
unsigned int Atom::GetCoordinationNumber(int nAtomicNo) const {
  return GetNumAtomsWithPredicate(GetBondedAtomList(m_bondMap),
                                  isAtomicNo_eq(nAtomicNo));
}

// This version returns the number of coordinated atoms of a given force field
// type
unsigned int Atom::GetCoordinationNumber(const std::string &strFFType) const {
  return GetNumAtomsWithPredicate(GetBondedAtomList(m_bondMap),
                                  isFFType_eq(strFFType));
}

// This version returns the number of coordinated atoms of a given hybridisation
// state
unsigned int Atom::GetCoordinationNumber(eHybridState e) const {
  return GetNumAtomsWithPredicate(GetBondedAtomList(m_bondMap),
                                  isHybridState_eq(e));
}

//////////////////////////////////////////
// Non-member functions (in rxdock namespace)
//////////////////////////////////////////

// DM 23 Apr 1999 - provide custom comparison function
// to bond map, so that map is sorted by bond ID, not by pointer
bool rxdock::BondPCmp_BondId::operator()(Bond *pBond1, Bond *pBond2) const {
  return pBond1->GetBondId() < pBond2->GetBondId();
}

////////////////////////////////////////////
// Conversion functions
///////////////////////

// Converts hybridisation state enum to a string
std::string rxdock::ConvertHybridStateToString(Atom::eHybridState eState) {
  switch (eState) {
  case Atom::SP:
    return "SP";
  case Atom::SP2:
    return "SP2";
  case Atom::SP3:
    return "SP3";
  case Atom::AROM:
    return "AROM";
  case Atom::TRI:
    return "TRI";
  default:
    return "";
  }
}

// Convert formal charge to a string (e.g. +, -, ++, --, +3, -3 etc)
std::string rxdock::ConvertFormalChargeToString(int nCharge) {
  switch (nCharge) {
  case 0:
    return "";
  case 1:
    return "+";
  case -1:
    return "-";
  case 2:
    return "++";
  case -2:
    return "--";
  default:
    // For higher charges, return as +3,-3 etc
    std::ostringstream ostr;
    ostr.setf(std::ios_base::showpos);
    ostr << nCharge;
    std::string strCharge(ostr.str());
    return strCharge;
  }
}

////////////////////////////////////////////
// Predicates (for use by STL algorithms)
////////////////////////////////////////////

// Is atom a H-Bond Acceptor ?
// Checks atomic number for O(8), N(7) or S(16)
// DM 15 Feb 1999 - screens out positively charged atoms
// Checks nitrogen hybridisation state and coordination number for available
// lone pairs
bool rxdock::isAtomHBondAcceptor::operator()(const Atom *pAtom) const {
  // Screen out positively charged atoms
  // if (spAtom->GetFormalCharge() > 0.0)
  //  return false;

  const int nAtomicNo = pAtom->GetAtomicNo();
  const Atom::eHybridState hybrid = pAtom->GetHybridState();
  const unsigned int nCoord = pAtom->GetCoordinationNumber();
  switch (nAtomicNo) {
  case 7: // N
    // Only include nitrogens which have available lone pairs
    switch (hybrid) {
    case Atom::SP3:
      return (nCoord < 4);
    case Atom::SP2:
    case Atom::TRI:
    case Atom::AROM:
      return (nCoord < 3);
    case Atom::SP:
      return (nCoord < 2);
    default:
      return false;
    }
  case 8:
    switch (hybrid) {
    case Atom::SP3:
      return (nCoord < 3);
    case Atom::SP2:
    case Atom::TRI:
    case Atom::AROM:
      // Remove terminal oxygens bonded to nitrogen (e.g. Nitro oxygens)
      return ((nCoord > 1) || (pAtom->GetCoordinationNumber(7) == 0));
    default:
      return false;
    }
  case 16:
    // 16 May 2003 (DM) Sulphur - only include terminal S_SP2, e.g. in C=S
    // Previously, all S were included
    switch (hybrid) {
    case Atom::SP2:
      return (nCoord < 2);
    default:
      return false;
    }
  default:
    return false;
  }
}

// Is atom a H-Bond Donor ?
// Checks 1) is it a hydrogen, 2) does it make exactly one bond, 3) is the bond
// to O, N or S ?
bool rxdock::isAtomHBondDonor::operator()(const Atom *pAtom) const {
  return ((pAtom->GetAtomicNo() == 1) && (pAtom->GetNumBonds() == 1) &&
          ((pAtom->GetCoordinationNumber(8) == 1) ||
           (pAtom->GetCoordinationNumber(7) == 1) ||
           (pAtom->GetCoordinationNumber(16) == 1)));
}

// Is atom planar ?
// Checks if 1) atom makes 2 bonds (in which case must be planar) or
// 2) hybridisation state is SP2, AROM or TRI
bool rxdock::isAtomPlanar::operator()(const Atom *pAtom) const {
  const Atom::eHybridState hybrid = pAtom->GetHybridState();
  // 17 Dec 1998 (DM) GetCoordinationNumber works for extended atoms
  // as it includes any implicit hydrogens
  //  return ( (spAtom->GetNumBonds() == 2) ||
  return ((pAtom->GetCoordinationNumber() == 2) || (hybrid == Atom::SP2) ||
          (hybrid == Atom::AROM) || (hybrid == Atom::TRI));
}

// Is atom a pi-atom ?
// SP2,TRI or AROM
bool rxdock::isPiAtom::operator()(const Atom *pAtom) const {
  const Atom::eHybridState hybrid = pAtom->GetHybridState();
  return ((hybrid == Atom::SP2) || (hybrid == Atom::AROM) ||
          (hybrid == Atom::TRI));
}

// DM 29 Jul 1999
// Is atom within a given distance of any coord in the coord list
bool rxdock::isAtomNearCoordList::operator()(const Atom *pAtom) const {
  // DM 27 Oct 2000 - GetCoords now returns by reference
  const Coord &atomCoord = pAtom->GetCoords();
  for (CoordListConstIter iter = cl.begin(); iter != cl.end(); iter++) {
    if (Length2(atomCoord - (*iter)) <= r2)
      return true;
  }
  return false;
}

// Is atom2 1-2 connected (i.e. bonded) to atom1 ?
rxdock::isAtom_12Connected::isAtom_12Connected(Atom *pAtom) : pAtom1(pAtom) {
  bondedAtomList1 = GetBondedAtomList(pAtom1);
}

bool rxdock::isAtom_12Connected::operator()(Atom *pAtom2) const {
  // Check if atom1 and atom2 are the same
  if (pAtom1 == pAtom2)
    return false;

  // Checks if atom2 appears in the bonded atom list for atom1.
  // STL find_if will do the trick
  // Need to use std::bind to convert binary predicate to unary
  return std::find_if(bondedAtomList1.begin(), bondedAtomList1.end(),
                      std::bind(isAtomPtr_eq(), std::placeholders::_1,
                                pAtom2)) != bondedAtomList1.end();
}

// Is atom2 1-3 connected (i.e. via a bond angle) to atom1 ?
rxdock::isAtom_13Connected::isAtom_13Connected(Atom *pAtom) : pAtom1(pAtom) {
  bondedAtomList1 = GetBondedAtomList(pAtom1);
}

bool rxdock::isAtom_13Connected::operator()(Atom *pAtom2) const {
  // Check if atom1 and atom2 are the same
  if (pAtom1 == pAtom2)
    return false;

  AtomList bondedAtomList2 = GetBondedAtomList(pAtom2);
  // Check if any atoms appear in both bonded atom lists
  // STL find_first_of will do the trick.
  return std::find_first_of(bondedAtomList1.begin(), bondedAtomList1.end(),
                            bondedAtomList2.begin(), bondedAtomList2.end(),
                            isAtomPtr_eq()) != bondedAtomList1.end();
}

// Is atom's subunit an RNA-type (i.e. A,G,C or U)
// DM 4 Jan 1999 - check for 3 letter residue names also
// Note: strictly speaking these are used for DNA, not RNA, but they often crop
// up in RNA
bool rxdock::isAtomRNA::operator()(const Atom *pAtom) const {
  std::string strSubunitName = pAtom->GetSubunitName();
  return ((strSubunitName == "A") || (strSubunitName == "ADE") ||
          (strSubunitName == "G") || (strSubunitName == "GUA") ||
          (strSubunitName == "C") || (strSubunitName == "CYT") ||
          (strSubunitName == "U") || (strSubunitName == "URI"));
}

// DM 21 Jul 1999 Is atom lipophilic ?
// DM 16 May 2003 Total rewrite to be much more comprehensive
bool rxdock::isAtomLipophilic::operator()(Atom *pAtom) const {
  // Remove the obvious polar stuff
  if (isIonic(pAtom) || isHBD(pAtom) || isHBA(pAtom) || isMetal(pAtom) ||
      isO(pAtom) || isN(pAtom))
    return false;

  const int nAtomicNo = pAtom->GetAtomicNo();
  AtomList bondedAtomList;
  AtomListConstIter iter;

  switch (nAtomicNo) {
  case 1: // H
    return true;
  case 6:  // C
  case 16: // S
    if (isSP3(pAtom))
      return true;
    bondedAtomList = GetBondedAtomList(pAtom);
    // Remove C,S bonded to O_SP2 (gets rid of carbonyl, sulphone, sulphoxide,
    // sulphonamide etc)
    for (iter = bondedAtomList.begin(); iter != bondedAtomList.end(); iter++) {
      // Check this way so that predicate works for PSF/CRD files (formal bond
      // order not set)
      if (isO(*iter) && isSP2(*iter))
        return false;
    }
    return true;
  case 9: // F
    return false;
  case 17: // Cl
  case 35: // Br
  case 53: // I
    return true;
  default:
    return false;
  }
}

// DM 24 Jan 2001
// Checks for common metal ions by atomic number (Na,Mg,K->Zn)
rxdock::isAtomMetal::isAtomMetal() {}
bool rxdock::isAtomMetal::operator()(const Atom *pAtom) const {
  int atNo = pAtom->GetAtomicNo();
  return (atNo == 11) || (atNo == 12) || ((atNo >= 19) && (atNo <= 30));
}

// Select the intramolecular flexible interactions to this atom (sets
// SelectionFlag = true)
void rxdock::SelectFlexAtoms::operator()(Atom *pAtom) {
  Model *pModel = pAtom->GetModelPtr();
  if (pModel) {
    pModel->SelectFlexAtoms(pAtom);
  }
}

////////////////////////////////////////////
// Atom list functions (implemented as STL algorithms)
////////////////////////////////////////////

void rxdock::SetAtomSelectionFlagsInList(AtomList &atomList, bool bSelected) {
  SelectAtom select(bSelected);
  std::for_each(atomList.begin(), atomList.end(), select);
}

void rxdock::InvertAtomSelectionFlags(AtomList &atomList) {
  InvertSelectAtom invert;
  std::for_each(atomList.begin(), atomList.end(), invert);
}

void rxdock::SetAtomCyclicFlagsInList(AtomList &atomList, bool bCyclic) {
  CyclicAtom cyclic(bCyclic);
  std::for_each(atomList.begin(), atomList.end(), cyclic);
}

// BondedAtoms
// Helper function to return the atom pointer for the "other" atom in the bond
// i.e. return (bondPair.second) ? bondPair.first->GetAtom2Ptr() :
// bondPair.first->GetAtom1Ptr();
AtomPtr rxdock::GetBondedAtomPtr(std::pair<Bond *, bool> bondBoolPair) {
  // If this atom is first in the bond (bondPair.second == true), get atom 2
  // ptr, else get atom 1 ptr
  return (bondBoolPair.second) ? (bondBoolPair.first)->GetAtom2Ptr()
                               : (bondBoolPair.first)->GetAtom1Ptr();
}

unsigned int rxdock::GetNumBondedAtoms(const BondMap &bondMap) {
  return bondMap.size(); // This function is not really needed, but provide for
                         // consistency
}

AtomList rxdock::GetBondedAtomList(const BondMap &bondMap) {
  AtomList atomList;
  std::transform(
      bondMap.begin(), bondMap.end(), std::back_inserter(atomList),
      std::function<AtomPtr(std::pair<Bond *, bool>)>(&GetBondedAtomPtr));
  // DM 23 Apr 1999 - sort the bonded atom list by atom ID so we return the
  // atoms in the same order every time The problem with the bond map is that
  // the key is Bond* which is not guaranteed to increase in order of
  // construction.
  std::sort(atomList.begin(), atomList.end(), AtomPtrCmp_AtomId());
  return atomList;
}

// If spAtom is a regular Atom, these two functions just behave like the two
// above If spAtom can be dynamically_casted to an PseudoAtom, these return
// the constituent atom list for the pseudoatom
unsigned int rxdock::GetNumBondedAtoms(const Atom *pAtom) {
  // Attempt the dynamic cast (on the regular pointer)
  PseudoAtom *pPseudoAtom =
      dynamic_cast<PseudoAtom *>(const_cast<Atom *>(pAtom));
  if (pPseudoAtom != nullptr)
    return pPseudoAtom->GetNumAtoms();
  else
    return GetNumBondedAtoms(pAtom->GetBondMap());
}

AtomList rxdock::GetBondedAtomList(const Atom *pAtom) {
  // Attempt the dynamic cast (on the regular pointer)
  PseudoAtom *pPseudoAtom =
      dynamic_cast<PseudoAtom *>(const_cast<Atom *>(pAtom));
  if (pPseudoAtom != nullptr)
    return pPseudoAtom->GetAtomList();
  else
    return GetBondedAtomList(pAtom->GetBondMap());
}

// Binary

// DM 1 Feb 1999
// Attempts to match atoms in atomList2 with those in atomList1, where the match
// is performed using isAtom_eq (tests subunit ID, subunit name and atom name
// for equality) and NOT isAtomPtr_eq (which tests the Atom* pointers
// themselves
// for equality, i.e. same objects) Returns list of Atom smart pointers to
// atoms in atomList1 for which a match is found.
unsigned int rxdock::GetNumMatchingAtoms(const AtomList &atomList1,
                                         const AtomList &atomList2) {
  unsigned int nCount(0);
  for (AtomListConstIter iter = atomList2.begin(); iter != atomList2.end();
       iter++) {
    nCount += GetNumAtomsWithPredicate(
        atomList1, std::bind(isAtom_eq(), std::placeholders::_1, *iter));
  }
  return nCount;
}

AtomList rxdock::GetMatchingAtomList(const AtomList &atomList1,
                                     const AtomList &atomList2) {
  AtomList matchingAtomList;
  for (AtomListConstIter iter = atomList2.begin(); iter != atomList2.end();
       iter++) {
    AtomList tempList = GetAtomListWithPredicate(
        atomList1, std::bind(isAtom_eq(), std::placeholders::_1, *iter));
    std::copy(tempList.begin(), tempList.end(),
              std::back_inserter(matchingAtomList));
  }
  return matchingAtomList;
}

// DM 30 Mar 1999
// These versions match atoms against a full atom name specifier of the form
// given by GetFullAtomName(): <segment name>:<subunit name>_<subunit id>:<atom
// name>  e.g. A:U_23:O4 All fields are optional, and if missing, will match all
// atoms So, for example, these are all valid full names: N7     (equiv. to
// ::N7)
// - matches N7's in all subunits, all segments A::O4                   -
// matches O4's in all subunits in segment A U:O4 (equiv. to :U:O4)  - matches
// O4's in subunits named U in all segments U_23:O4                 - matches
// O4's in subunit U23 in all segments _23:                    - matches all
// atoms in
// subunit ID=23 in all segments
unsigned int rxdock::GetNumMatchingAtoms(const AtomList &atomList,
                                         const std::string &strFullName) {
  return GetMatchingAtomList(atomList, strFullName).size();
}

AtomList rxdock::GetMatchingAtomList(const AtomList &atomList,
                                     const std::string &strFullName) {
  AtomList matchingAtoms = atomList;

  // Split the name into it's constituent subunit and atom names
  std::vector<std::string> componentList =
      ConvertDelimitedStringToList(strFullName, ":");
  unsigned int idx(0); // Index into constituent list
  // Switch on how many constituent names are specified
  switch (componentList.size()) {

  case 3: // we have segment, subunit and atom names
    if (!componentList[idx].empty()) {
      matchingAtoms = GetAtomListWithPredicate(
          matchingAtoms, isSegmentName_eq(componentList[idx]));
      LOG_F(1, "Matching segment name={}, #atoms={}", componentList[idx],
            matchingAtoms.size());
    }
    idx++;

  case 2: // Fall-through!! We have subunit and atom names
    // Check if the subunit name also specifies a subunit ID
    if (!componentList[idx].empty()) {
      std::vector<std::string> subunitList =
          ConvertDelimitedStringToList(componentList[idx], "_");
      switch (subunitList.size()) {
      case 2: // Subunit ID and name are specified
        if (!subunitList[1].empty()) {
          matchingAtoms = GetAtomListWithPredicate(
              matchingAtoms, isSubunitId_eq(subunitList[1]));
          LOG_F(1, "Matching subunit id={}, #atoms={}", subunitList[1],
                matchingAtoms.size());
        }
      case 1: // Fall-through!! Only the subunit name is specified
        if (!subunitList[0].empty()) {
          matchingAtoms = GetAtomListWithPredicate(
              matchingAtoms, isSubunitName_eq(subunitList[0]));
          LOG_F(1, "Matching subunit name={}, #atoms={}", subunitList[0],
                matchingAtoms.size());
        }
        break;
      default:
        LOG_F(1, "Invalid subunit string in {}", strFullName);
        break;
      }
    }
    idx++;

  case 1: // Fall-through!! We only have an atom name specifier
    if (!componentList[idx].empty()) {
      matchingAtoms = GetAtomListWithPredicate(
          matchingAtoms, isAtomName_eq(componentList[idx]));
      LOG_F(1, "Matching atom name={}, #atoms={}", componentList[idx],
            matchingAtoms.size());
    }
    break;

  default:
    LOG_IF_F(WARNING, !strFullName.empty(), "Too many colons (:) in {}",
             strFullName);
    break;
  }

  return matchingAtoms;
}

// DM 15 Apr 1999 - as above, but match against a list of full atom name
// specifiers Returns total list (i.e. all matches OR'd). Does not remove
// duplicates.
unsigned int
rxdock::GetNumMatchingAtoms(const AtomList &atomList,
                            const std::vector<std::string> &fullNameList) {
  return GetMatchingAtomList(atomList, fullNameList).size();
}

AtomList
rxdock::GetMatchingAtomList(const AtomList &atomList,
                            const std::vector<std::string> &fullNameList) {
  AtomList matchingAtoms;
  for (std::vector<std::string>::const_iterator iter = fullNameList.begin();
       iter != fullNameList.end(); iter++) {
    // Find the atoms which match the specifier string
    AtomList selectedAtoms = GetMatchingAtomList(atomList, *iter);
    // Add to the total list of matching atoms
    std::copy(selectedAtoms.begin(), selectedAtoms.end(),
              std::back_inserter(matchingAtoms));
  }
  return matchingAtoms;
}

/////////////////////////
// DM 07 Jan 1999
// Actions on atom lists
/////////////////////////

// Save coords by number for all atoms in the list
void rxdock::SaveAtomCoords(const AtomList &atomList, unsigned int coordNum) {
  for (AtomListConstIter iter = atomList.begin(); iter != atomList.end();
       iter++) {
    AtomPtr spAtom(*iter);
    spAtom->SaveCoords(coordNum);
  }
}

// Revert to numbered coords for all atoms in the list
void rxdock::RevertAtomCoords(const AtomList &atomList, unsigned int coordNum) {
  for (AtomListConstIter iter = atomList.begin(); iter != atomList.end();
       iter++) {
    AtomPtr spAtom(*iter);
    spAtom->RevertCoords(coordNum);
  }
}

// Returns total atomic mass (molecular weight) for all atoms in the list
double rxdock::GetTotalAtomicMass(const AtomList &atomList) {
  return std::accumulate(atomList.begin(), atomList.end(), 0.0,
                         AccumAtomicMass);
}

// Returns center of mass of atoms in the list
Coord rxdock::GetCenterOfAtomicMass(const AtomList &atomList) {
  // Accumulate sum of mass*coord
  Coord com = std::accumulate(atomList.begin(), atomList.end(), Coord(),
                              AccumMassWeightedCoords);
  // Divide by total mass
  com /= GetTotalAtomicMass(atomList);
  return com;
}

// DM 20 May 1999 - returns the coords of all atoms in the list
CoordList rxdock::GetCoordList(const AtomList &atomList) {
  CoordList coordList;
  std::transform(atomList.begin(), atomList.end(),
                 std::back_inserter(coordList), ExtractAtomCoord);
  return coordList;
}

// DM 09 Aug 2001 - returns coordList via argument
void rxdock::GetCoordList(const AtomList &atomList, CoordList &coordList) {
  coordList.clear();
  std::transform(atomList.begin(), atomList.end(),
                 std::back_inserter(coordList), ExtractAtomCoord);
}

// Streams an atom list in Quanta CSD file format (for easy color coding of
// selected atoms in Quanta) nFormat  = 0 => Receptor atom format: "zone A:1 #
// pick O5T = col 2" nFormat != 0 => Ligand atom format "pick N1 .and. segm H =
// col 2" std::ostream should have been opened before calling QuantaCSDFileDump
// DM 30 Jul 1999 - added segment name to receptor atom format
void rxdock::PrintQuantaCSDFormat(const AtomList &atomList, std::ostream &s,
                                  int nColor, int nFormat) {
  for (AtomListConstIter iter = atomList.begin(); iter != atomList.end();
       iter++) {
    AtomPtr spAtom(*iter);
    switch (nFormat) {
    case 0:
      s << "zone " << spAtom->GetSegmentName() << ":" << spAtom->GetSubunitId()
        << " # pick " << spAtom->GetName() << " = col " << nColor << std::endl;
      break;
    default:
      s << "pick " << spAtom->GetName() << " .and. segm "
        << spAtom->GetSegmentName() << " = col " << nColor << std::endl;
    }
  }
}

// DM 25 Feb 1999
// Modified DM 6 Apr 1999 to operate using the GroupCharge attribute
// Scans an atom list and sets the group charge to zero on any
// cationic-ionic pairs which are 1-2 or 1-3 connected
//(e.g. on -OOC.CR.NH3+ would zero the C and N)
void rxdock::RemoveZwitterions(AtomList &atomList) {
  AtomList cationicAtomList =
      GetAtomListWithPredicate(atomList, isAtomCationic());
  AtomList anionicAtomList =
      GetAtomListWithPredicate(atomList, isAtomAnionic());

  // Progressively partition the cationic atom list into those that are 1-3
  // connected or less to one of the anionic atoms, and those that aren't
  // c14begin points to the first cationic atom which is 1-4 connected or above
  // to all of the anionic atoms All cationic atoms up to c14begin can be
  // removed
  AtomListIter c14begin = cationicAtomList.begin();
  for (AtomListConstIter aIter = anionicAtomList.begin();
       aIter != anionicAtomList.end(); aIter++) {
    // Need to partition in two stages: firstly 1-2 connected, then 1-3
    // connected
    c14begin = std::stable_partition(c14begin, cationicAtomList.end(),
                                     isAtom_12Connected(*aIter));
    c14begin = std::stable_partition(c14begin, cationicAtomList.end(),
                                     isAtom_13Connected(*aIter));
  }

  // Now repeat the other way round, checking the anionic atom list against ALL
  // the cationic atoms (including those we know are due to be removed) We need
  // to do it in this way to cope with situations like (O-)-(N+)-(O-) We
  // couldn't simply delete the first N,O pair we find as this would leave the
  // second O with nothing to match.
  //
  AtomListIter a14begin = anionicAtomList.begin();
  for (AtomListConstIter cIter = cationicAtomList.begin();
       cIter != cationicAtomList.end(); cIter++) {
    // Need to partition in two stages: firstly 1-2 connected, then 1-3
    // connected
    a14begin = std::stable_partition(a14begin, anionicAtomList.end(),
                                     isAtom_12Connected(*cIter));
    a14begin = std::stable_partition(a14begin, anionicAtomList.end(),
                                     isAtom_13Connected(*cIter));
  }

  // Now we can safely neutralise the 1-2 and 1-3 atoms
  for (AtomListIter iter = cationicAtomList.begin(); iter != c14begin; iter++) {
    (*iter)->SetGroupCharge(0.0);
    LOG_F(1, "Removing zwitterions: neutralising cation {}",
          (*iter)->GetFullAtomName());
  }

  for (AtomListIter iter = anionicAtomList.begin(); iter != a14begin; iter++) {
    (*iter)->SetGroupCharge(0.0);
    LOG_F(1, "Removing zwitterions: neutralising anion {}",
          (*iter)->GetFullAtomName());
  }
}
