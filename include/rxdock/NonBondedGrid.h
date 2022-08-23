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

// Grid for indexing lipophilic interactions

#ifndef _RBTNONBONDEDGRID_H_
#define _RBTNONBONDEDGRID_H_

#include "rxdock/Atom.h"
#include "rxdock/BaseGrid.h"

#include <nlohmann/json.hpp>

using json = nlohmann::json;

namespace rxdock {

// A map of atom vectors indexed by unsigned int
// Used to store the receptor atom lists at each grid point
// DM 11 Jul 2000 - use map of regular Atom* list (not AtomPtr smart
// pointer list)

// DM 3 Nov 2000 - replace map by vector for faster lookup
// typedef std::map<UInt,AtomRList> AtomListMap;
typedef std::vector<AtomRList> AtomListMap;

typedef AtomListMap::iterator AtomListMapIter;
typedef AtomListMap::const_iterator AtomListMapConstIter;

class NonBondedGrid : public BaseGrid {
public:
  // Class type string
  static const std::string _CT;
  ////////////////////////////////////////
  // Constructors/destructors
  // Construct a NXxNYxNZ grid running from gridMin at gridStep resolution
  NonBondedGrid(const Coord &gridMin, const Coord &gridStep, unsigned int NX,
                unsigned int NY, unsigned int NZ, unsigned int NPad = 0);

  // Constructor reading all params from binary stream
  NonBondedGrid(json j);

  ~NonBondedGrid(); // Default destructor

  /*friend void to_json(json &j, const NonBondedGrid &nonBondedGrid);
  friend void from_json(const json &j, NonBondedGrid &nonBondedGrid);*/

  // Copy constructor
  NonBondedGrid(const NonBondedGrid &);
  // Copy constructor taking base class argument
  NonBondedGrid(const BaseGrid &);
  // Copy assignment
  NonBondedGrid &operator=(const NonBondedGrid &);
  // Copy assignment taking base class argument
  NonBondedGrid &operator=(const BaseGrid &);

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
  // AtomList GetAtomList(UInt iXYZ) const;
  // AtomList GetAtomList(const Coord& c) const;
  const AtomRList &GetAtomList(unsigned int iXYZ) const;
  const AtomRList &GetAtomList(const Coord &c) const;

  /////////////////////////
  // Set attribute functions
  /////////////////////////
  void SetAtomLists(Atom *pAtom, double radius);
  void ClearAtomLists();
  void UniqueAtomLists();

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
  NonBondedGrid(); // Disable default constructor

  // Helper function called by copy constructor and assignment operator
  void CopyGrid(const NonBondedGrid &);
  // DM 6 Nov 2000 - create AtomListMap of the appropriate size
  void CreateMap();

protected:
  ////////////////////////////////////////
  // Protected data
  ////////////////

private:
  ////////////////////////////////////////
  // Private data
  //////////////
  AtomListMap
      m_atomMap; // Used to store the receptor atom lists at each grid point
  const AtomRList m_emptyList; // Dummy list used by GetAtomList
};

/*void to_json(json &j, const NonBondedGrid &nonBondedGrid);
void from_json(const json &j, NonBondedGrid &nonBondedGrid);*/

// Useful typedefs
typedef SmartPtr<NonBondedGrid> NonBondedGridPtr; // Smart pointer

} // namespace rxdock

#endif //_RBTNONBONDEDGRID_H_
