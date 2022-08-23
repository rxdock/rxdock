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

// Grid for indexing interactions involving more than 2 atoms

#ifndef _RBTINTERACTIONGRID_H_
#define _RBTINTERACTIONGRID_H_

#include "rxdock/Atom.h"
#include "rxdock/BaseGrid.h"

#include <nlohmann/json.hpp>

using json = nlohmann::json;

namespace rxdock {

// simple container for up to 3 atoms, to hold one half of an interaction
// i.e. receptor atoms or ligand atoms.
class InteractionCenter {
public:
  enum eLP { NONE, PLANE, LONEPAIR };
  InteractionCenter(Atom *pAtom1 = nullptr, Atom *pAtom2 = nullptr,
                    Atom *pAtom3 = nullptr, eLP LP = NONE)
      : m_pAtom1(pAtom1), m_pAtom2(pAtom2), m_pAtom3(pAtom3), m_LP(LP) {
    _RBTOBJECTCOUNTER_CONSTR_("InteractionCenter");
  }
  ~InteractionCenter() { _RBTOBJECTCOUNTER_DESTR_("InteractionCenter"); }

  friend void to_json(json &j, const InteractionCenter &interactionCenter);
  friend void from_json(const json &j, InteractionCenter &interactionCenter);

  Atom *GetAtom1Ptr() const { return m_pAtom1; }
  Atom *GetAtom2Ptr() const { return m_pAtom2; }
  Atom *GetAtom3Ptr() const { return m_pAtom3; }
  eLP LP() const { return m_LP; }
  // Returns list of constituent atoms
  //(deconvolutes pseudoatoms into their constituent Atom* lists)
  AtomRList GetAtomList() const;
  bool isSelected() const;

private:
  // Could be a useful general function
  // If pAtom is a pseudo atom, then pushes all the constituent atoms onto list
  // else, push pAtom onto the list
  void AccumulateAtomList(const Atom *pAtom, AtomRList &atomList) const;
  Atom *m_pAtom1;
  Atom *m_pAtom2;
  Atom *m_pAtom3;
  eLP m_LP;
};

void to_json(json &j, const InteractionCenter &interactionCenter);
void from_json(const json &j, InteractionCenter &interactionCenter);

typedef std::vector<InteractionCenter *>
    InteractionCenterList; // Vector of regular pointers
typedef InteractionCenterList::iterator InteractionCenterListIter;
typedef InteractionCenterList::const_iterator InteractionCenterListConstIter;

// Less than operator for sorting InteractionCenter* by pointer value
class InteractionCenterCmp {
public:
  bool operator()(const InteractionCenter *pIC1,
                  const InteractionCenter *pIC2) const {
    return pIC1 < pIC2;
  }
};
// Is interaction center selected ?
class isInteractionCenterSelected
    : public std::function<bool(InteractionCenter *)> {
public:
  explicit isInteractionCenterSelected() {}
  bool operator()(const InteractionCenter *pIC) const {
    return pIC->isSelected();
  }
};

// Is the distance between interaction centers less than a given value ?
// Function checks d**2 to save performing a sqrt
class isInteractionD_lt : public std::function<bool(InteractionCenter *)> {
  double d_sq;
  Atom *a;

public:
  explicit isInteractionD_lt(const InteractionCenter *pIC1, double dd)
      : d_sq(dd * dd), a(pIC1->GetAtom1Ptr()) {}
  bool operator()(const InteractionCenter *pIC2) const {
    return Length2(pIC2->GetAtom1Ptr()->GetCoords(), a->GetCoords()) < d_sq;
  }
};

// Select/deselect the interaction center (selects all constituent atoms)
class SelectInteractionCenter {
  bool b;

public:
  explicit SelectInteractionCenter(bool bb) : b(bb) {}
  void operator()(InteractionCenter *pIC);
};

// A map of interaction centers indexed by unsigned int
// Used to store the receptor atom lists at each grid point
// DM 11 Jul 2000 - use map of regular Atom* list (not AtomPtr smart
// pointer list)

// DM 3 Nov 2000 - replace map by vector for faster lookup
// typedef std::map<UInt,InteractionCenterList> InteractionListMap;
typedef std::vector<InteractionCenterList> InteractionListMap;

typedef InteractionListMap::iterator InteractionListMapIter;
typedef InteractionListMap::const_iterator InteractionListMapConstIter;

class InteractionGrid : public BaseGrid {
public:
  // Class type string
  static const std::string _CT;
  ////////////////////////////////////////
  // Constructors/destructors
  // Construct a NXxNYxNZ grid running from gridMin at gridStep resolution
  InteractionGrid(const Coord &gridMin, const Coord &gridStep, unsigned int NX,
                  unsigned int NY, unsigned int NZ, unsigned int NPad = 0);

  // Constructor reading all params from binary stream
  InteractionGrid(json j);

  ~InteractionGrid(); // Default destructor

  /*friend void to_json(json &j, const InteractionGrid &interactionGrid);
  friend void from_json(const json &j, InteractionGrid &interactionGrid);*/

  // Copy constructor
  InteractionGrid(const InteractionGrid &);
  // Copy constructor taking base class argument
  InteractionGrid(const BaseGrid &);
  // Copy assignment
  InteractionGrid &operator=(const InteractionGrid &);
  // Copy assignment taking base class argument
  InteractionGrid &operator=(const BaseGrid &);

  ////////////////////////////////////////
  // Virtual functions for reading/writing grid data to streams in
  // text and binary format
  // Subclasses should provide their own private OwnPrint
  // method to handle subclass data members, and override the public
  // Print method
  virtual void Print(std::ostream &ostr) const; // Text output

  ////////////////////////////////////////
  // Public methods
  ////////////////

  /////////////////////////
  // Get attribute functions
  /////////////////////////
  const InteractionCenterList &GetInteractionList(unsigned int iXYZ) const;
  const InteractionCenterList &GetInteractionList(const Coord &c) const;

  /////////////////////////
  // Set attribute functions
  /////////////////////////
  void SetInteractionLists(InteractionCenter *pIntn, double radius);
  void ClearInteractionLists();
  void UniqueInteractionLists();

protected:
  ////////////////////////////////////////
  // Protected methods
  ///////////////////
  // Protected method for writing data members for this class to text stream
  void OwnPrint(std::ostream &ostr) const;

private:
  ////////////////////////////////////////
  // Private methods
  /////////////////
  InteractionGrid(); // Disable default constructor

  // Helper function called by copy constructor and assignment operator
  void CopyGrid(const InteractionGrid &);
  // DM 3 Nov 2000 - create InteractionListMap of the appropriate size
  void CreateMap();

protected:
  ////////////////////////////////////////
  // Protected data
  ////////////////

private:
  ////////////////////////////////////////
  // Private data
  //////////////
  InteractionListMap m_intnMap; // Used to store the interaction center lists
                                // at each grid point
  const InteractionCenterList m_emptyList; // Dummy list used by GetAtomList
};
/*
void to_json(json &j, const InteractionGrid &interactionGrid);
void from_json(const json &j, InteractionGrid &interactionGrid);*/

// Useful typedefs
typedef SmartPtr<InteractionGrid> InteractionGridPtr; // Smart pointer

} // namespace rxdock

#endif //_RBTINTERACTIONGRID_H_
