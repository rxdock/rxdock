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

// Handles all the properties for a single atom

#ifndef _RBTATOM_H_
#define _RBTATOM_H_

#include "rxdock/Config.h"
#include "rxdock/Coord.h"
#include "rxdock/PMF.h"
#include "rxdock/Quat.h"
#include "rxdock/TriposAtomType.h"

#include <list>

namespace rxdock {

class Model; // Forward declaration
class Bond;  // Forward declaration

// DM 23 Apr 1999 - provide custom comparison function
// to bond map, so that map is sorted by bond ID, not by pointer
class BondPCmp_BondId {
public:
  bool operator()(Bond *pBond1, Bond *pBond2) const;
};

// Map of bonds the atom is bonded to
// Key is regular Bond*, value is Bool
// Bool is true if this atom is atom1 in the bond, false if atom2
typedef std::map<Bond *, bool, BondPCmp_BondId> BondMap;
typedef BondMap::iterator BondMapIter;
typedef BondMap::const_iterator BondMapConstIter;

class Atom {
public:
  ///////////////////////////////////////////////
  // Enums
  ///////////////////////////////////////////////

  // Hybridisation state
  enum eHybridState {
    UNDEFINED = 0,
    SP = 1,   // SP hybridised (linear geom)
    SP2 = 2,  // SP2 hybridised (planar geom)
    SP3 = 3,  // SP3 hybridised (tetrahedral geom)
    AROM = 4, // Aromatic SP2 (planar geom)
    TRI = 5,  // Trigonal nitrogen (planar geom)
  };

  ///////////////////////////////////////////////
  // Constructors / destructors
  ///////////////////////////////////////////////

  // Default constructor
  RBTDLL_EXPORT Atom();

  // Constructor supplying all 2-D parameters
  RBTDLL_EXPORT
  explicit Atom(int nAtomId, int nAtomicNo = 6, std::string strAtomName = "C",
                std::string strSubunitId = "1",
                std::string strSubunitName = "RES",
                std::string strSegmentName = "SEG1",
                eHybridState eState =
                    UNDEFINED, // DM 8 Dec 1998 Changed from SP3 to UNDEFINED
                unsigned int nHydrogens = 0, int nFormalCharge = 0);

  // Default destructor
  virtual ~Atom();

  // Copy constructor
  Atom(const Atom &atom);

  // Copy assignment
  Atom &operator=(const Atom &atom);

  ///////////////////////////////////////////////
  // Stream functions
  ///////////////////////////////////////////////

  // Insertion operator (primarily for debugging)
  // Note: needs to be friend so can access atom private data
  // without using the accessor functions
  RBTDLL_EXPORT friend std::ostream &operator<<(std::ostream &s,
                                                const Atom &atom);

  // Virtual function for dumping atom details to an output stream
  // Derived classes (e.g. pseudoatom) can override if required
  virtual std::ostream &Print(std::ostream &s) const;

  ///////////////////////////////////////////////
  // Public accessor functions
  // 2-D attributes
  ///////////////////////////////////////////////

  // DM 7 June 2006
  // atom enabled state is controlled by the parent Model
  // enabled state
  RBTDLL_EXPORT bool GetEnabled() const;

  // AtomId
  int GetAtomId() const { return m_nAtomId; }
  void SetAtomId(const int nAtomId) { m_nAtomId = nAtomId; }

  // AtomicNo
  int GetAtomicNo() const { return m_nAtomicNo; }
  void SetAtomicNo(const int nAtomicNo) { m_nAtomicNo = nAtomicNo; }

  // AtomName
  std::string GetName() const { return m_strAtomName; }
  void SetAtomName(const std::string &strAtomName) {
    m_strAtomName = strAtomName;
  }
  RBTDLL_EXPORT std::string
  GetFullAtomName() const; // Returns composite of segment name,
                           // subunit id and name, and atom name

  // SubunitId
  std::string GetSubunitId() const { return m_strSubunitId; }
  void SetSubunitId(const std::string &strSubunitId) {
    m_strSubunitId = strSubunitId;
  }

  // SubunitName
  std::string GetSubunitName() const { return m_strSubunitName; }
  void SetSubunitName(const std::string &strSubunitName) {
    m_strSubunitName = strSubunitName;
  }

  // SegmentName
  std::string GetSegmentName() const { return m_strSegmentName; }
  void SetSegmentName(const std::string &strSegmentName) {
    m_strSegmentName = strSegmentName;
  }

  // Hybridisation state
  eHybridState GetHybridState() const { return m_eState; }
  void SetHybridState(const eHybridState eState) { m_eState = eState; }

  // Attached hydrogens
  unsigned int GetNumImplicitHydrogens() const { return m_nHydrogens; }
  void SetNumImplicitHydrogens(const unsigned int nHydrogens) {
    m_nHydrogens = nHydrogens;
  }

  // FormalCharge (DM 24 Mar 1999 - changed from double to int)
  int GetFormalCharge() const { return m_nFormalCharge; }
  void SetFormalCharge(const int nFormalCharge) {
    m_nFormalCharge = nFormalCharge;
  }

  // CyclicFlag - flag to indicate atom is in a ring (set by Model::FindRing)
  bool GetCyclicFlag() const { return m_bCyclic; }
  void SetCyclicFlag(bool bCyclic = true) { m_bCyclic = bCyclic; }

  // SelectionFlag - general purpose flag can be set/cleared by various search
  // algorithms (e.g. FindRings)
  bool GetSelectionFlag() const { return m_bSelected; }
  void SetSelectionFlag(bool bSelected = true) { m_bSelected = bSelected; }
  void InvertSelectionFlag() { m_bSelected = !m_bSelected; }

  // DM 21 Jul 1999 - user-defined flag
  bool GetUser1Flag() const { return m_bUser1; }
  void SetUser1Flag(bool bUser1 = true) { m_bUser1 = bUser1; }

  // DM 29 Jan 2000 - user-defined double
  double GetUser1Value() const { return m_dUser1; }
  void SetUser1Value(double dUser1 = 0.0) { m_dUser1 = dUser1; }

  // DM 27 Jul 2000 - user-defined double#2
  double GetUser2Value() const { return m_dUser2; }
  void SetUser2Value(double dUser2 = 0.0) { m_dUser2 = dUser2; }

  // ModelPtr - pointer to the parent model
  Model *GetModelPtr() const { return m_pModel; }
  void SetModelPtr(Model *pModel = nullptr) { m_pModel = pModel; }

  // DM 04 Dec 1998  Add functions to handle bond map
  // Returns number of bonds in map
  unsigned int GetNumBonds() const { return m_bondMap.size(); }
  // Returns the bond map
  // DM 26 Oct 2000 - return by reference
  // BondMap GetBondMap() const {return m_bondMap;}
  const BondMap &GetBondMap() const { return m_bondMap; }
  // Add a bond to the bond map - returns true if OK, false if this atom is not
  // a member of the bond or if bond already added
  bool AddBond(Bond *pBond);
  // Remove a bond from the bond map - returns true if OK, false if bond not
  // present in map
  bool RemoveBond(Bond *pBond);
  // Returns number of cyclic bonds in map
  unsigned int GetNumCyclicBonds() const;
  // Returns bond map with only cyclic bonds included
  BondMap GetCyclicBondMap() const;
  // Returns total formal bond order around atom
  int GetTotalFormalBondOrder() const;
  // Returns max formal bond order for all bonds to atom
  int GetMaxFormalBondOrder() const;

  ///////////////////////////////////////////////
  // Public accessor functions
  // 3-D attributes
  ///////////////////////////////////////////////

  // Coords
  // DM 8 Dec 1998 - make virtual so we can override in the pseudoatom class
  // DM 11 Jul 2000 - remove overhead of virtual methods
  // Pseudo atoms must now refresh their coords each time the constituent atoms
  // move DM 27 Oct 2000 - return by reference
  const Coord &GetCoords() const { return m_coord; }
  double GetX() const { return m_coord.xyz(0); }
  double GetY() const { return m_coord.xyz(1); }
  double GetZ() const { return m_coord.xyz(2); }

  void SetCoords(const Coord &coord) { m_coord = coord; }
  void SetCoords(const double x, const double y, const double z) {
    m_coord.xyz(0) = x;
    m_coord.xyz(1) = y;
    m_coord.xyz(2) = z;
  }
  void SetX(const double x) { m_coord.xyz(0) = x; }
  void SetY(const double y) { m_coord.xyz(1) = y; }
  void SetZ(const double z) { m_coord.xyz(2) = z; }

  // PartialCharge
  double GetPartialCharge() const { return m_dPartialCharge; }
  void SetPartialCharge(const double dPartialCharge) {
    m_dPartialCharge = dPartialCharge;
  }

  // GroupCharge (added DM 24 Mar 1999, for ionic interaction group charges)
  double GetGroupCharge() const { return m_dGroupCharge; }
  void SetGroupCharge(const double dGroupCharge) {
    m_dGroupCharge = dGroupCharge;
  }

  // AtomicMass
  double GetAtomicMass() const { return m_dAtomicMass; }
  void SetAtomicMass(const double dAtomicMass) { m_dAtomicMass = dAtomicMass; }

  // VdwRadius
  double GetVdwRadius() const { return m_dVdwRadius; }
  void SetVdwRadius(const double dVdwRadius) { m_dVdwRadius = dVdwRadius; }

  // AtomType
  std::string GetFFType() const { return m_strFFType; }
  void SetFFType(const std::string &strFFType) { m_strFFType = strFFType; }
  PMFType GetPMFType() const { return m_nPMFType; }
  void SetPMFType(PMFType aType) { m_nPMFType = aType; }
  TriposAtomType::eType GetTriposType() const { return m_triposType; }
  void SetTriposType(TriposAtomType::eType aType) { m_triposType = aType; }

  // XB
  // reweighting factor
  // double GetReweight() const { return m_dReweight; }
  // void SetReweight(const double dReweight) { m_dReweight = dReweight; }
  // XB END MODIFICATIONS

  ///////////////////////////////////////////////
  // Other public methods
  ///////////////////////////////////////////////
  // DM 11 Jan 1999 - give atoms the ability to save their current coords
  // and revert to the saved coords. Useful for Monte Carlo simulations.
  // DM 08 Feb 1999 - all saved coords are now saved in a
  // std::map<UInt,Coord> map key=0 is reserved for the default SaveCoords
  // and RevertCoords
  void SaveCoords(unsigned int coordNum = 0);
  void RevertCoords(unsigned int coordNum = 0);

  // Translate - translate coordinates by the supplied vector
  void Translate(const Vector &vector) { m_coord += vector; }

  void Translate(const double vx, const double vy, const double vz) {
    m_coord += Coord(vx, vy, vz);
  }

  // DM 07 Jan 1999 - rotate coordinates using the supplied quaternion
  void RotateUsingQuat(const Quat &q) { m_coord = q.Rotate(m_coord); }

  // DM 04 Dec 1998  Now we have the bond map, we can easily provide
  // coordination numbers This version returns the total number of coordinated
  // atoms (includes implicit hydrogens) Equivalent to
  // GetNumBonds()+GetNumImplicitHydrogens()
  unsigned int GetCoordinationNumber() const;
  // This version returns the number of coordinated atoms of a given element
  // Note: for hydrogens, the implicit atoms are excluded
  // so GetCoordinationNumber(1) will return the number of explicit hydrogens
  unsigned int GetCoordinationNumber(int nAtomicNo) const;
  // This version returns the number of coordinated atoms of a given force field
  // type
  unsigned int GetCoordinationNumber(const std::string &strFFType) const;
  // This version returns the number of coordinated atoms of a given
  // hybridisation state
  unsigned int GetCoordinationNumber(eHybridState e) const;

protected:
private:
  // Private methods
  // Clears the bond map - should only need to be called by the copy
  // constructors, hence private
  void ClearBondMap();

private:
  // Private data

  // These can be considered as 2-D params (i.e. define the chemistry and
  // topology)
  int m_nAtomicNo;              // Atomic number
  int m_nAtomId;                // Atom ID
  std::string m_strAtomName;    // Atom name
  std::string m_strSubunitId;   // Subunit(residue) ID (note: string not int)
  std::string m_strSubunitName; // Subunit(residue) name
  std::string m_strSegmentName; // Segment(molecule) name
  eHybridState m_eState;        // Hybridisation state
  unsigned int m_nHydrogens;    // Number of attached (implicit) hydrogens
  int m_nFormalCharge; // Formal charge (DM 24 Mar 1999 - changed from double
                       // to int)
  Model *m_pModel;     // Regular pointer to parent model
  BondMap m_bondMap;   // Map of bonds this atom is bonded to
  bool m_bCyclic;      // Is the atom in a ring ?
  bool m_bSelected;    // Can be set/cleared by various search algorithms (e.g.
                       // FindRings)
  bool m_bUser1;       // User flag 1
  double m_dUser1;     // User value 1
  double m_dUser2;     // User value 2
  PMFType m_nPMFType;  // PMF atom type: default
  TriposAtomType::eType m_triposType; // DM 14 May 2002: Tripos 5.2 atom type

  // These can be considered as 3-D params (i.e. the extra info required for 3-D
  // calculations)
  Coord m_coord;           // cartesian coords
  double m_dPartialCharge; // partial charge
  double m_dGroupCharge;   // interaction group charge (added DM 24 Mar 1999)
  double m_dAtomicMass;    // atomic mass
  double m_dVdwRadius;     // atomic mass
  std::string m_strFFType; // force field atom type
  // double m_dReweight;      // XB reweighting factor

  UIntCoordMap m_savedCoords; // DM 08 Feb 1999 - now store all saved coords
                              // in a std::map<UInt,Coord>
};

// Useful typedefs
typedef SmartPtr<Atom> AtomPtr;        // Smart pointer
typedef std::vector<AtomPtr> AtomList; // Vector of smart pointers
typedef AtomList::iterator AtomListIter;
typedef AtomList::const_iterator AtomListConstIter;

typedef std::vector<Atom *> AtomRList; // Vector of regular pointers
typedef AtomRList::iterator AtomRListIter;
typedef AtomRList::const_iterator AtomRListConstIter;

typedef std::vector<AtomList>
    AtomListList; // A vector of atom vectors (e.g. for storing ring systems)
typedef AtomListList::iterator AtomListListIter;
typedef AtomListList::const_iterator AtomListListConstIter;

typedef std::vector<AtomRList>
    AtomRListList; // A vector of atom vectors (regular pointer version)
typedef AtomRListList::iterator AtomRListListIter;
typedef AtomRListList::const_iterator AtomRListListConstIter;

typedef std::list<AtomPtr> AtomTrueList; // A real list of smart pointers
typedef AtomTrueList::iterator AtomTrueListIter;
typedef AtomTrueList::const_iterator AtomTrueListConstIter;

///////////////////////////////////////////////
// Non-member functions (in rxdock namespace)
//////////////////////////////////////////

std::ostream &operator<<(std::ostream &s, const Atom &atom);

////////////////////////////////////////////
// Calculation functions
///////////////////////

// DM 07 Jul 2000 - converted to regular Atom* versions
// Returns the "bond" length in Angstroms between 2 atoms
// Note: atoms do not have to be bonded
inline double BondLength(Atom *pAtom1, Atom *pAtom2) {
  return Length(pAtom1->GetCoords(), pAtom2->GetCoords());
}

// Returns the "bond" angle in degrees between 3 atoms
// Note: atoms do not have to be bonded
inline double BondAngle(Atom *pAtom1, Atom *pAtom2, Atom *pAtom3) {
  return Angle(pAtom1->GetCoords(), pAtom2->GetCoords(), pAtom3->GetCoords());
}

// Returns the "bond" dihedral in degrees between 4 atoms
// Note: atoms do not have to be bonded
inline double BondDihedral(Atom *pAtom1, Atom *pAtom2, Atom *pAtom3,
                           Atom *pAtom4) {
  return Dihedral(pAtom1->GetCoords(), pAtom2->GetCoords(), pAtom3->GetCoords(),
                  pAtom4->GetCoords());
}

////////////////////////////////////////////
// Conversion functions
///////////////////////

// Converts hybridisation state enum to a string
std::string ConvertHybridStateToString(Atom::eHybridState eState);
// Convert formal charge to a string (e.g. +, -, +2, -2 etc)
std::string ConvertFormalChargeToString(int nCharge);

////////////////////////////////////////////
// Predicates (for use by STL algorithms)
////////////////////////////////////////////
// DM 30 Apr 2002 - convert to regular Atom* functions
// More universal as AtomPtr parameters will be automatically degraded to
// regular pointers

// Is atom enabled ?
class isAtomEnabled : public std::function<bool(Atom *)> {
public:
  explicit isAtomEnabled() = default;
  bool operator()(Atom *pAtom) const { return pAtom->GetEnabled(); }
};

// Is atom selected ?
class isAtomSelected : public std::function<bool(Atom *)> {
public:
  explicit isAtomSelected() = default;
  bool operator()(Atom *pAtom) const { return pAtom->GetSelectionFlag(); }
};

// Is atom cyclic ?
class isAtomCyclic : public std::function<bool(Atom *)> {
public:
  explicit isAtomCyclic() = default;
  bool operator()(const Atom *pAtom) const { return pAtom->GetCyclicFlag(); }
};

// DM 9 Feb 2000 - is atom bridgehead? Checks for >2 cyclic bonds
class isAtomBridgehead : public std::function<bool(Atom *)> {
public:
  explicit isAtomBridgehead() = default;
  bool operator()(const Atom *pAtom) const {
    return (pAtom->GetNumCyclicBonds() > 2);
  }
};

// Is atom a H-Bond Acceptor ?
// Checks atomic number for O(8), N(7) or S(16)
class isAtomHBondAcceptor : public std::function<bool(Atom *)> {
public:
  explicit isAtomHBondAcceptor() = default;
  RBTDLL_EXPORT bool operator()(const Atom *) const;
};

// Is atom a H-Bond Donor ?
// Checks 1) is it a hydrogen, 2) does it make exactly one bond, 3) is the bond
// to O, N or S ?
class isAtomHBondDonor : public std::function<bool(Atom *)> {
public:
  explicit isAtomHBondDonor() = default;
  RBTDLL_EXPORT bool operator()(const Atom *) const;
};

// Is atom formally charged ?
// Checks if formal charge is != zero
class isAtomCharged : public std::function<bool(Atom *)> {
public:
  explicit isAtomCharged() = default;
  bool operator()(const Atom *pAtom) const {
    return (pAtom->GetFormalCharge() != 0);
  }
};

// Is atom formally positively charged ?
// Checks if formal charge is > zero
class isAtomPosCharged : public std::function<bool(Atom *)> {
public:
  explicit isAtomPosCharged() = default;
  bool operator()(const Atom *pAtom) const {
    return (pAtom->GetFormalCharge() > 0);
  }
};

// Is atom formally negatively charged ?
// Checks if formal charge is < zero
class isAtomNegCharged : public std::function<bool(Atom *)> {
public:
  explicit isAtomNegCharged() = default;
  bool operator()(const Atom *pAtom) const {
    return (pAtom->GetFormalCharge() < 0);
  }
};

// Does atom have implicit hydrogens
class isAtomExtended : public std::function<bool(Atom *)> {
public:
  explicit isAtomExtended() = default;
  bool operator()(const Atom *pAtom) const {
    return (pAtom->GetNumImplicitHydrogens() > 0);
  }
};

// Is atom planar ?
// Checks if 1) atom makes 2 bonds (in which case must be planar) or
// 2) hybridisation state is SP2, AROM or TRI
class isAtomPlanar : public std::function<bool(Atom *)> {
public:
  explicit isAtomPlanar() = default;
  bool operator()(const Atom *) const;
};

// Is atom a pi-atom ?
// SP2,TRI or AROM, or special case of OSP3
class isPiAtom : public std::function<bool(Atom *)> {
public:
  explicit isPiAtom() = default;
  RBTDLL_EXPORT bool operator()(const Atom *) const;
};

// Is atomic number equal to n ?
class isAtomicNo_eq : public std::function<bool(Atom *)> {
  int n;

public:
  explicit isAtomicNo_eq(int nn) : n(nn) {}
  bool operator()(const Atom *pAtom) const { return pAtom->GetAtomicNo() == n; }
};

// Is Force field type equal to s ?
class isFFType_eq : public std::function<bool(Atom *)> {
  std::string s;

public:
  explicit isFFType_eq(std::string ss) : s(std::move(ss)) {}
  bool operator()(const Atom *pAtom) const { return pAtom->GetFFType() == s; }
};

// Is Atom name equal to s ?
class isAtomName_eq : public std::function<bool(Atom *)> {
  std::string s;

public:
  explicit isAtomName_eq(std::string ss) : s(std::move(ss)) {}
  bool operator()(const Atom *pAtom) const { return pAtom->GetName() == s; }
};

// Is Subunit name equal to s ?
class isSubunitName_eq : public std::function<bool(Atom *)> {
  std::string s;

public:
  explicit isSubunitName_eq(std::string ss) : s(std::move(ss)) {}
  bool operator()(const Atom *pAtom) const {
    return pAtom->GetSubunitName() == s;
  }
};

// Is Subunit ID equal to s ?
class isSubunitId_eq : public std::function<bool(Atom *)> {
  std::string s;

public:
  explicit isSubunitId_eq(std::string ss) : s(std::move(ss)) {}
  bool operator()(const Atom *pAtom) const {
    return pAtom->GetSubunitId() == s;
  }
};

// Is Segment name equal to s ?
class isSegmentName_eq : public std::function<bool(Atom *)> {
  std::string s;

public:
  explicit isSegmentName_eq(std::string ss) : s(std::move(ss)) {}
  bool operator()(const Atom *pAtom) const {
    return pAtom->GetSegmentName() == s;
  }
};

// Is hybridisation state equal to e ? (DM 26 Mar 1999)
class isHybridState_eq : public std::function<bool(Atom *)> {
  Atom::eHybridState e;

public:
  explicit isHybridState_eq(Atom::eHybridState ee) : e(ee) {}
  bool operator()(const Atom *pAtom) const {
    return pAtom->GetHybridState() == e;
  }
};

// DM 6 Jan 1999
// Is atom2 equal to atom1 (checks if underlying regular pointers match)
// Note: this is a binary rather than unary predicate
// DM 1 Feb 1999 - renamed from isAtom_eq to isAtomPtr_eq
class isAtomPtr_eq : public std::function<bool(Atom *, Atom *)> {
public:
  explicit isAtomPtr_eq() = default;
  bool operator()(const Atom *pAtom1, const Atom *pAtom2) const {
    return pAtom1 == pAtom2;
  }
};

// DM 1 Feb 1999
// Is atom2 equal to atom1 (checks if subunit name, subunit ID and atom name
// match) Note: this is a binary rather than unary predicate
class isAtom_eq : public std::function<bool(Atom *, Atom *)> {
public:
  explicit isAtom_eq() = default;
  bool operator()(const Atom *pAtom1, const Atom *pAtom2) const {
    return ((pAtom1->GetSubunitId() == pAtom2->GetSubunitId()) &&
            (pAtom1->GetSubunitName() == pAtom2->GetSubunitName()) &&
            (pAtom1->GetName() == pAtom2->GetName()));
  }
};

// DM 22 Apr 2002
// Is atom2 ID equal to atom1 ID
// Note: this is a binary rather than unary predicate
class isAtomId_eq : public std::function<bool(Atom *, Atom *)> {
public:
  explicit isAtomId_eq() = default;
  bool operator()(const Atom *pAtom1, const Atom *pAtom2) const {
    return (pAtom1->GetAtomId() == pAtom2->GetAtomId());
  }
};

// DM 12 Apr 1999
// Is atom inside a sphere defined by radius=r, center=c
class isAtomInsideSphere : public std::function<bool(Atom *)> {
  const Coord &c; // center of sphere
  double r2;      // radius squared (to avoid taking square roots)
public:
  explicit isAtomInsideSphere(const Coord &cc, double rr)
      : c(cc), r2(rr * rr) {}
  bool operator()(const Atom *pAtom) const {
    return Length2(pAtom->GetCoords() - c) <= r2;
  }
};

// DM 12 Apr 1999
// Is atom within a cuboid defined by cmin,cmax coords?
class isAtomInsideCuboid : public std::function<bool(Atom *)> {
  const Coord &cmin;
  const Coord &cmax;

public:
  explicit isAtomInsideCuboid(const Coord &ccmin, const Coord &ccmax)
      : cmin(ccmin), cmax(ccmax) {}
  bool operator()(const Atom *pAtom) const {
    // DM 27 Oct 2000 - GetCoords now returns by reference
    const Coord &c = pAtom->GetCoords();
    return (c >= cmin) && (c <= cmax);
  }
};

// DM 29 Jul 1999
// Is atom within a given distance of any coord in the coord list
class isAtomNearCoordList : public std::function<bool(Atom *)> {
  const CoordList &cl;
  double r2; // radius squared (to avoid taking square roots)
public:
  explicit isAtomNearCoordList(const CoordList &ccl, double rr)
      : cl(ccl), r2(rr * rr) {}
  bool operator()(const Atom *pAtom) const;
};

// Is atom2 1-2 connected (i.e. bonded) to atom1 ?
class isAtom_12Connected : public std::function<bool(Atom *)> {
  Atom *pAtom1;
  AtomList bondedAtomList1;

public:
  explicit isAtom_12Connected(Atom *spAtom);
  bool operator()(Atom *spAtom2) const;
};

// Is atom2 1-3 connected (i.e. via a bond angle) to atom1 ?
class isAtom_13Connected : public std::function<bool(Atom *)> {
  Atom *pAtom1;
  AtomList bondedAtomList1;

public:
  explicit isAtom_13Connected(Atom *pAtom);
  bool operator()(Atom *pAtom2) const;
};

// Is atom's subunit an RNA-type (i.e. A,G,C or U)
class isAtomRNA : public std::function<bool(Atom *)> {
public:
  explicit isAtomRNA() = default;
  bool operator()(const Atom *) const;
};

// Is atom classified as an RNA receptor "anionic" atom
// For our purposes, these are N7,O2,O4,O6,O1P and O2P
// DM 22 Dec 1998 - center phosphate interactions on the P atom, not on O1P and
// O2P to avoid double counting of favourable interactions class
// isAtomRNAAnionic : public std::function<bool(Atom *)>  {
//  isAtomRNA bIsAtomRNA;
// public:
//  explicit isAtomRNAAnionic() {}
//  Bool operator() (AtomPtr) const;
//};

// DM 22 Dec 1998 - need ligand anionic atoms so we can include repulsive
// interactions with the RNA
//
// Is atom classified as a ligand "anionic" atom
// Includes COO-, PO3-
// class isAtomLigandAnionic : public std::function<bool(Atom *)>  {
// public:
//  explicit isAtomLigandAnionic() {}
//  Bool operator() (AtomPtr) const;
//};

// DM 22 Dec 1998 - renamed from isAtomRNACationic
//
// Is atom classified as a ligand "cationic" atom
// The most important category is charged Nitrogens
// class isAtomLigandCationic : public std::function<bool(Atom *)>  {
// public:
//  explicit isAtomLigandCationic() {}
//  Bool operator() (AtomPtr) const;
//};

//////////////////////////////
// DM 30 Mar 1999
// Now that the MolecularFileSources define the GroupCharge attribute in a
// file-format independent manner, we don't need separate predicates for
// isRNAAnionic, isLigandAnionic etc We just need isIonic, isCationic, isAnionic

// Is atom defined as an ionic interaction center ?
// Checks if group charge is != zero
class isAtomIonic : public std::function<bool(Atom *)> {
public:
  explicit isAtomIonic() = default;
  bool operator()(const Atom *pAtom) const {
    return (std::fabs(pAtom->GetGroupCharge()) > 0.001);
  }
};

// Is atom defined as a cationic interaction center ?
// Checks if group charge is > zero
class isAtomCationic : public std::function<bool(Atom *)> {
public:
  explicit isAtomCationic() = default;
  bool operator()(const Atom *pAtom) const {
    return (pAtom->GetGroupCharge() > 0.001);
  }
};

// Is atom defined as an anionic interaction center ?
// Checks if group charge is < zero
class isAtomAnionic : public std::function<bool(Atom *)> {
public:
  explicit isAtomAnionic() = default;
  bool operator()(const Atom *pAtom) const {
    return (pAtom->GetGroupCharge() < -0.001);
  }
};

// Is atom the central carbon in a guanidinium group ?
// Checks for cationic sp2/arom carbon
// DM 27 Jan 2000 - also check is acyclic
class isAtomGuanidiniumCarbon : public std::function<bool(Atom *)> {
  isAtomicNo_eq bIsCarbon;
  isAtomCationic bIsCationic;
  isPiAtom bIsPi;
  isAtomCyclic bIsCyclic;

public:
  explicit isAtomGuanidiniumCarbon() : bIsCarbon(6) {}
  bool operator()(const Atom *pAtom) const {
    return bIsCationic(pAtom) && bIsCarbon(pAtom) && bIsPi(pAtom) &&
           !bIsCyclic(pAtom);
  }
};

// DM 24 Jan 2001
// Checks for common metal ions by atomic number
// DM 19 Oct 2001 - extend to all common metals (Na,Mg,K->Zn)
class isAtomMetal : public std::function<bool(Atom *)> {
public:
  RBTDLL_EXPORT explicit isAtomMetal();
  RBTDLL_EXPORT bool operator()(const Atom *pAtom) const;
};

// DM 21 Jul 1999 Is atom lipophilic ?
// DM 16 May 2003 Total rewrite to be much more comprehensive
class isAtomLipophilic : public std::function<bool(Atom *)> {
  isAtomIonic isIonic;
  isAtomHBondDonor isHBD;
  isAtomHBondAcceptor isHBA;
  isAtomMetal isMetal;
  isHybridState_eq isSP3;
  isHybridState_eq isSP2;
  isAtomicNo_eq isO;
  isAtomicNo_eq isN;

public:
  explicit isAtomLipophilic()
      : isSP3(Atom::SP3), isSP2(Atom::SP2), isO(8), isN(7) {}
  RBTDLL_EXPORT bool operator()(Atom *pAtom) const;
};

//////////////////////////////

// DM 6 April 1999
// Is coordination number of atom equal to n
// First constructor checks total coordination number
// Other constructors check atomic number, force field type and hybridisation
// state coordination numbers
class isCoordinationNumber_eq : public std::function<bool(Atom *)> {
  enum {
    TOTAL,
    ATNO,
    FFTYPE,
    HYBRID
  } eCNType;      // Type of coordination number to check
  unsigned int n; // Coordination number value to check
  int atNo;
  std::string ffType;
  Atom::eHybridState hybrid;

public:
  // Total coordination number
  explicit isCoordinationNumber_eq(unsigned int nn)
      : eCNType(TOTAL), n(nn), atNo(0), ffType(""), hybrid(Atom::UNDEFINED) {}
  // Atomic number coordination number
  explicit isCoordinationNumber_eq(unsigned int nn, int nAt)
      : eCNType(ATNO), n(nn), atNo(nAt), ffType(""), hybrid(Atom::UNDEFINED) {}
  // Force field type coordination number
  explicit isCoordinationNumber_eq(unsigned int nn, std::string strType)
      : eCNType(FFTYPE), n(nn), atNo(0), ffType(std::move(strType)),
        hybrid(Atom::UNDEFINED) {}
  // Hybridisation state coordination number
  explicit isCoordinationNumber_eq(unsigned int nn, Atom::eHybridState eState)
      : eCNType(HYBRID), n(nn), atNo(0), ffType(""), hybrid(eState) {}
  bool operator()(Atom *pAtom) const {
    switch (eCNType) {
    case TOTAL:
      return pAtom->GetCoordinationNumber() == n;
    case ATNO:
      return pAtom->GetCoordinationNumber(atNo) == n;
    case FFTYPE:
      return pAtom->GetCoordinationNumber(ffType) == n;
    case HYBRID:
      return pAtom->GetCoordinationNumber(hybrid) == n;
    default:
      return false;
    }
  }
};

////////////////////////////////////////////////////////
// Comparison functions for sorting Atom* containers
// For use by STL sort algorithms
////////////////////////////////////////////////////////

// Less than operator for sorting Atom*s by atom ID
class AtomPtrCmp_AtomId {
public:
  bool operator()(Atom *pAtom1, Atom *pAtom2) const {
    return pAtom1->GetAtomId() < pAtom2->GetAtomId();
  }
};

// Less than operator for sorting Atom*s by atomic number
class AtomPtrCmp_AtomicNo {
public:
  bool operator()(Atom *pAtom1, Atom *pAtom2) const {
    return pAtom1->GetAtomicNo() < pAtom2->GetAtomicNo();
  }
};

// Less than operator for sorting Atom*s by atom name
class AtomPtrCmp_AtomName {
public:
  bool operator()(Atom *pAtom1, Atom *pAtom2) const {
    return pAtom1->GetName() < pAtom2->GetName();
  }
};

// Less than operator for sorting AtomPtrs by pointer value
class AtomPtrCmp_Ptr {
public:
  bool operator()(Atom *pAtom1, Atom *pAtom2) const { return pAtom1 < pAtom2; }
};

////////////////////////////////////////////
// Functions objects for performing actions on atoms
// For use by STL algorithms
////////////////////////////////////////////

// DM 28 Jul 1999 - extract atom coord (for use by std::transform)
// DM 27 Oct 2000 - return by reference
inline const Coord &ExtractAtomCoord(Atom *pAtom) { return pAtom->GetCoords(); }

// DM 09 Nov 1999 - Accumulate atomic mass (for use by std::accumulate)
inline double AccumAtomicMass(double val, Atom *pAtom) {
  return val + pAtom->GetAtomicMass();
}

// DM 11 Nov 1999 - Accumulate mass weighted coords (for use by std::accumulate)
inline Coord AccumMassWeightedCoords(const Coord &val, Atom *pAtom) {
  return val + pAtom->GetAtomicMass() * pAtom->GetCoords();
}

// Translate an atom by the supplied vector
class TranslateAtom {
  const Vector v;

public:
  explicit TranslateAtom(const Vector &vv) : v(vv) {}
  void operator()(Atom *pAtom) { pAtom->Translate(v); }
};

// Translate an atom by the supplied vector, but only if the selection flag is
// true
class TranslateAtomIfSelected {
  const Vector v;

public:
  explicit TranslateAtomIfSelected(const Vector &vv) : v(vv) {}
  void operator()(Atom *pAtom) {
    if (pAtom->GetSelectionFlag())
      pAtom->Translate(v);
  }
};

// Rotate an atom by the supplied quaternion
class RotateAtomUsingQuat {
  const Quat q;

public:
  explicit RotateAtomUsingQuat(const Quat &qq) : q(qq) {}
  void operator()(Atom *pAtom) { pAtom->RotateUsingQuat(q); }
};

// Rotate an atom by the supplied quaternion, but only if the selection flag is
// true
class RotateAtomUsingQuatIfSelected {
  const Quat q;

public:
  explicit RotateAtomUsingQuatIfSelected(const Quat &qq) : q(qq) {}
  void operator()(Atom *pAtom) {
    if (pAtom->GetSelectionFlag())
      pAtom->RotateUsingQuat(q);
  }
};

// Select/deselect the atom
class SelectAtom {
  bool b;

public:
  explicit SelectAtom(bool bb) : b(bb) {}
  void operator()(Atom *pAtom) { pAtom->SetSelectionFlag(b); }
};

// Invert the atom selection flag
class InvertSelectAtom {
public:
  explicit InvertSelectAtom() = default;
  void operator()(Atom *pAtom) { pAtom->InvertSelectionFlag(); }
};

// Select the intramolecular flexible interactions to this atom (sets
// SelectionFlag = true)
class SelectFlexAtoms {
public:
  explicit SelectFlexAtoms() = default;
  void operator()(Atom *pAtom);
};

// Set/unset the cyclic flag
class CyclicAtom {
  bool b;

public:
  explicit CyclicAtom(bool bb) : b(bb) {}
  void operator()(Atom *pAtom) { pAtom->SetCyclicFlag(b); }
};

////////////////////////////////////////////
// Atom list functions (implemented as STL algorithms)
////////////////////////////////////////////

// Unary

// Generic template version of GetNumAtoms, passing in your own predicate
template <class Predicate>
unsigned int GetNumAtomsWithPredicate(const AtomList &atomList,
                                      const Predicate &pred) {
  return std::count_if(atomList.begin(), atomList.end(), pred);
}

// Generic template version of GetAtomList, passing in your own predicate
template <class Predicate>
AtomList GetAtomListWithPredicate(const AtomList &atomList,
                                  const Predicate &pred) {
  AtomList newAtomList;
  std::copy_if(atomList.begin(), atomList.end(),
               std::back_inserter(newAtomList), pred);
  return newAtomList;
}

// Generic template version of FindAtomInList, passing in your own predicate
template <class Predicate>
AtomListIter FindAtomInList(AtomList &atomList, const Predicate &pred) {
  return std::find_if(atomList.begin(), atomList.end(), pred);
}

// Selected atoms
void SetAtomSelectionFlagsInList(AtomList &atomList, bool bSelected = true);
void InvertAtomSelectionFlags(AtomList &atomList); // DM 08 Jan 1999

inline unsigned int GetNumSelectedAtomsInList(const AtomList &atomList) {
  return GetNumAtomsWithPredicate(atomList, isAtomSelected());
}

inline AtomList GetSelectedAtomsFromList(const AtomList &atomList) {
  return GetAtomListWithPredicate(atomList, isAtomSelected());
}

// Cyclic atoms
void SetAtomCyclicFlagsInList(AtomList &atomList, bool bCyclic = true);

inline unsigned int GetNumCyclicAtomsInList(const AtomList &atomList) {
  return GetNumAtomsWithPredicate(atomList, isAtomCyclic());
}
inline AtomList GetCyclicAtomsFromList(const AtomList &atomList) {
  return GetAtomListWithPredicate(atomList, isAtomCyclic());
}

// Hydrogen bond acceptor atoms
inline unsigned int GetNumHBondAcceptorAtomsInList(const AtomList &atomList) {
  return GetNumAtomsWithPredicate(atomList, isAtomHBondAcceptor());
}
inline AtomList GetHBondAcceptorAtomsFromList(const AtomList &atomList) {
  return GetAtomListWithPredicate(atomList, isAtomHBondAcceptor());
}

// Hydrogen bond donor atoms
inline unsigned int GetNumHBondDonorAtomsInList(const AtomList &atomList) {
  return GetNumAtomsWithPredicate(atomList, isAtomHBondDonor());
}
inline AtomList GetHBondDonorAtomsFromList(const AtomList &atomList) {
  return GetAtomListWithPredicate(atomList, isAtomHBondDonor());
}

//(Formally) charged atoms
inline unsigned int GetNumChargedAtomsInList(const AtomList &atomList) {
  return GetNumAtomsWithPredicate(atomList, isAtomCharged());
}
inline AtomList GetChargedAtomsFromList(const AtomList &atomList) {
  return GetAtomListWithPredicate(atomList, isAtomCharged());
}

// Planar atoms
inline unsigned int GetNumPlanarAtomsInList(const AtomList &atomList) {
  return GetNumAtomsWithPredicate(atomList, isAtomPlanar());
}
inline AtomList GetPlanarAtomsFromList(const AtomList &atomList) {
  return GetAtomListWithPredicate(atomList, isAtomPlanar());
}

// RNA anionic atoms
// inline UInt GetNumRNAAnionicAtoms(const AtomList& atomList) {
//  return
//  GetNumAtomsWithPredicate(atomList,isAtomRNAAnionic());
//}
// inline AtomList GetRNAAnionicAtomList(const AtomList& atomList) {
//  return
//  GetAtomListWithPredicate(atomList,isAtomRNAAnionic());
//}

// Ligand anionic atoms
// DM 22 Dec 1998
// inline UInt GetNumLigandAnionicAtoms(const AtomList& atomList) {
//  return
//  GetNumAtomsWithPredicate(atomList,isAtomLigandAnionic());
//}
// inline AtomList GetLigandAnionicAtomList(const AtomList& atomList) {
//  return
//  GetAtomListWithPredicate(atomList,isAtomLigandAnionic());
//}

// Ligand cationic atoms
// DM 22 Dec 1998 -renamed from Get..RNACationic..
// inline UInt GetNumLigandCationicAtoms(const AtomList& atomList) {
//  return
//  GetNumAtomsWithPredicate(atomList,isAtomLigandCationic());
//}
// inline AtomList GetLigandCationicAtomList(const AtomList& atomList) {
//  return
//  GetAtomListWithPredicate(atomList,isAtomLigandCationic());
//}

// Guanidinium carbons
// inline UInt GetNumGuanidiniumCarbonAtoms(const AtomList& atomList) {
//  return
//  GetNumAtomsWithPredicate(atomList,isAtomGuanidiniumCarbon());
//}
// inline AtomList GetGuanidiniumCarbonAtomList(const AtomList& atomList)
// {
//  return
//  GetAtomListWithPredicate(atomList,isAtomGuanidiniumCarbon());
//}

// BondedAtoms
// Helper function to return the atom pointer for the "other" atom in the bond
// i.e. return (bondPair.second) ? bondPair.first->GetAtom2Ptr() :
// bondPair.first->GetAtom1Ptr();
AtomPtr GetBondedAtomPtr(std::pair<Bond *, bool> bondBoolPair);

// This pair of functions is a little different to the above
// They convert a BondMap to an AtomList. The BondMap is assumed to
// contain all bonds to a particular atom GetBondedAtomPtr is used to get the
// "other" atom in the pair.
unsigned int GetNumBondedAtoms(const BondMap &bondMap);
AtomList GetBondedAtomList(const BondMap &bondMap);

// If spAtom is a regular Atom, these two functions just behave like the two
// above If spAtom can be dynamically_casted to an PseudoAtom, these return
// the constituent atom list for the pseudoatom
unsigned int GetNumBondedAtoms(const Atom *pAtom);
RBTDLL_EXPORT AtomList GetBondedAtomList(const Atom *pAtom);

// Binary

// Atoms with atomic no = nAtomicNo
inline unsigned int GetNumAtomsWithAtomicNo_eq(const AtomList &atomList,
                                               int nAtomicNo) {
  return GetNumAtomsWithPredicate(atomList, isAtomicNo_eq(nAtomicNo));
}
inline AtomList GetAtomListWithAtomicNo_eq(const AtomList &atomList,
                                           int nAtomicNo) {
  return GetAtomListWithPredicate(atomList, isAtomicNo_eq(nAtomicNo));
}

// Atoms with FFType = strFFType
inline unsigned int GetNumAtomsWithFFType_eq(const AtomList &atomList,
                                             std::string strFFType) {
  return GetNumAtomsWithPredicate(atomList, isFFType_eq(std::move(strFFType)));
}
inline AtomList GetAtomListWithFFType_eq(const AtomList &atomList,
                                         std::string strFFType) {
  return GetAtomListWithPredicate(atomList, isFFType_eq(std::move(strFFType)));
}

// Atoms with AtomName = strAtomName
inline unsigned int GetNumAtomsWithAtomName_eq(const AtomList &atomList,
                                               std::string strAtomName) {
  return GetNumAtomsWithPredicate(atomList,
                                  isAtomName_eq(std::move(strAtomName)));
}
inline AtomList GetAtomListWithAtomName_eq(const AtomList &atomList,
                                           std::string strAtomName) {
  return GetAtomListWithPredicate(atomList,
                                  isAtomName_eq(std::move(strAtomName)));
}

// DM 1 Feb 1999
// Attempts to match atoms in atomList2 with those in atomList1, where the match
// is performed using isAtom_eq (tests subunit ID, subunit name and atom name
// for equality) and NOT isAtomPtr_eq (which tests the Atom* pointers
// themselves
// for equality, i.e. same objects) Returns list of Atom smart pointers to
// atoms in atomList1 for which a match is found.
unsigned int GetNumMatchingAtoms(const AtomList &atomList1,
                                 const AtomList &atomList2);
AtomList GetMatchingAtomList(const AtomList &atomList1,
                             const AtomList &atomList2);

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
unsigned int GetNumMatchingAtoms(const AtomList &atomList,
                                 const std::string &strFullName);
AtomList GetMatchingAtomList(const AtomList &atomList,
                             const std::string &strFullName);

// DM 15 Apr 1999 - as above, but match against a list of full atom name
// specifiers Returns total list (i.e. all matches OR'd). Does not remove
// duplicates.
unsigned int GetNumMatchingAtoms(const AtomList &atomList,
                                 const std::vector<std::string> &fullNameList);
AtomList GetMatchingAtomList(const AtomList &atomList,
                             const std::vector<std::string> &fullNameList);

/////////////////////////
// DM 07 Jan 1999
// Actions on atom lists
/////////////////////////
// Translate all atoms in the list by the supplied vector
inline void TranslateAtoms(const AtomList &atomList, const Vector &v) {
  std::for_each(atomList.begin(), atomList.end(), TranslateAtom(v));
}

// Translate all selected atoms in the list by the supplied vector
inline void TranslateSelectedAtoms(const AtomList &atomList, const Vector &v) {
  std::for_each(atomList.begin(), atomList.end(), TranslateAtomIfSelected(v));
}

// Rotate all atoms in the list by the supplied quaternion
inline void RotateAtomsUsingQuat(const AtomList &atomList, const Quat &q) {
  std::for_each(atomList.begin(), atomList.end(), RotateAtomUsingQuat(q));
}

// Rotate all selected atoms in the list by the supplied quaternion
inline void RotateSelectedAtomsUsingQuat(const AtomList &atomList,
                                         const Quat &q) {
  std::for_each(atomList.begin(), atomList.end(),
                RotateAtomUsingQuatIfSelected(q));
}

// Save coords by number for all atoms in the list
void SaveAtomCoords(const AtomList &atomList, unsigned int coordNum = 0);
// Revert to numbered coords for all atoms in the list
void RevertAtomCoords(const AtomList &atomList, unsigned int coordNum = 0);

// DM 9 Nov 1999
// Returns total atomic mass (molecular weight) for all atoms in the list
double GetTotalAtomicMass(const AtomList &atomList);

// Returns center of mass of atoms in the list
Coord GetCenterOfAtomicMass(const AtomList &atomList);

// DM 20 May 1999 - returns the coords of all atoms in the list
CoordList GetCoordList(const AtomList &atomList);
// DM 09 Aug 2001 - returns coordList via argument
RBTDLL_EXPORT void GetCoordList(const AtomList &atomList, CoordList &coordList);

// Streams an atom list in Quanta CSD file format (for easy color coding of
// selected atoms in Quanta) nFormat  = 0 => Receptor atom format: "zone 1 #
// pick O5T = col 2" nFormat != 0 => Ligand atom format "pick N1 .and. segm H =
// col 2" std::ostream should have been opened before calling
// PrintQuantaCSDFormat
void PrintQuantaCSDFormat(const AtomList &atomList, std::ostream &s,
                          int nColor = 2, int nFormat = 0);

// DM 25 Feb 1999
// Modified DM 6 Apr 1999 to operate using the GroupCharge attribute
// Scans an atom list and sets the group charge to zero on any
// cationic-ionic pairs which are 1-2 or 1-3 connected
//(e.g. on -OOC.CR.NH3+ would zero the C and N)
void RemoveZwitterions(AtomList &atomList);

} // namespace rxdock

#endif //_RBTATOM_H_
