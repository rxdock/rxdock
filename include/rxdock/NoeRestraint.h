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

// Miscellaneous helper classes/structs for handling NOE interactions

#ifndef _RBTNOERESTRAINT_H_
#define _RBTNOERESTRAINT_H_

#include "rxdock/Atom.h"

namespace rxdock {

// NOE Restraint type
// UNDEFINED: undefined restraint type
// OR: The distance to each of the atoms in the list is calculated separately,
//    and any of these can satisfy the restraint.
//    In practice, this means only the shortest distance is scored
// MEAN: The distance to the mean coords of the atoms in the list is calculated,
//      and used to evaluate the restraint
// OR: The distance to each of the atoms in the list is calculated separately,
//    and all of these should satisfy the restraint.
//    In practice, this means only the longest distance is scored
enum eNoeType { NOE_UNDEFINED, NOE_OR, NOE_MEAN, NOE_AND };

// Returns a string representing the NOE restraint type
std::string ConvertNoeRestraintTypeToString(const eNoeType type);

///////////////////////////////////////////////////////////////////////
// Simple struct for holding one end of an NOE restraint
// Has an atom name string list and a restraint type
class NoeEndNames {
public:
  NoeEndNames() {}
  bool isOK() const;
  // Public data
  std::vector<std::string> names;
  eNoeType type;
};

// Insertion operator for the above
std::ostream &operator<<(std::ostream &s, const NoeEndNames &n);

///////////////////////////////////////////////////////////////////////
// As above, but with a real atom list in place of the string list
class NoeEndAtoms {
public:
  // Default constructor
  NoeEndAtoms() {}
  // Constructor accepting an NoeEndNames and an atom list to match against
  NoeEndAtoms(const NoeEndNames &n, const AtomList &atomList);
  bool isOK() const;
  // Public data
  AtomList atoms;
  eNoeType type;
};

// Insertion operator for the above
std::ostream &operator<<(std::ostream &s, const NoeEndAtoms &n);

///////////////////////////////////////////////////////////////////////
// Simple struct for holding a complete NOE restraint definition
// Has two NoeEndNames, and max distance
class NoeRestraintNames {
public:
  NoeRestraintNames() : maxDist(5.0) {}
  bool isOK() const;
  // Public data
  NoeEndNames from;
  NoeEndNames to;
  double maxDist;
};

// Insertion operator for the above
std::ostream &operator<<(std::ostream &s, const NoeRestraintNames &noe);

///////////////////////////////////////////////////////////////////////
// As above, but with real NoeEndAtoms in place of NoeEndNames
class NoeRestraintAtoms {
public:
  // Default constructor
  NoeRestraintAtoms() : maxDist(5.0) {}
  // Constructor accepting an NoeRestraintNames and an atom list to match
  // against
  NoeRestraintAtoms(const NoeRestraintNames &n, const AtomList &atomList);
  bool isOK() const; // true if at least one atom at each end of the NOE and
                     // NOE type is not undefined
  bool isSimple()
      const; // true if exactly one atom at each end of the NOE (unambiguous)
  // Public data
  NoeEndAtoms from;
  NoeEndAtoms to;
  // Double minDist;
  double maxDist;
};

// Insertion operator for the above
std::ostream &operator<<(std::ostream &s, const NoeRestraintAtoms &noe);

///////////////////////////////////////////////////////////////////////
// Simple struct for holding a complete STD restraint definition
// Has one NoeEndNames, and max distance (from receptor)
class StdRestraintNames {
public:
  StdRestraintNames() : maxDist(5.0) {}
  bool isOK() const;
  // Public data
  NoeEndNames from;
  double maxDist;
};

// Insertion operator for the above
std::ostream &operator<<(std::ostream &s, const StdRestraintNames &std);

///////////////////////////////////////////////////////////////////////
// As above, but with real NoeEndAtoms in place of NoeEndNames
class StdRestraintAtoms {
public:
  // Default constructor
  StdRestraintAtoms() : maxDist(5.0) {}
  // Constructor accepting an StdRestraintNames and an atom list to match
  // against
  StdRestraintAtoms(const StdRestraintNames &n, const AtomList &atomList);
  bool isOK() const; // true if at least one atom in the definition and type
                     // is not undefined
  bool isSimple()
      const; // true if exactly one atom at in the STD definition (unambiguous)
  // Public data
  NoeEndAtoms from;
  double maxDist;
};

// Insertion operator for the above
std::ostream &operator<<(std::ostream &s, const StdRestraintAtoms &std);

// Useful typedefs
typedef std::vector<NoeRestraintNames> NoeRestraintNamesList;
typedef NoeRestraintNamesList::iterator NoeRestraintNamesListIter;
typedef NoeRestraintNamesList::const_iterator NoeRestraintNamesListConstIter;

typedef std::vector<NoeRestraintAtoms> NoeRestraintAtomsList;
typedef NoeRestraintAtomsList::iterator NoeRestraintAtomsListIter;
typedef NoeRestraintAtomsList::const_iterator NoeRestraintAtomsListConstIter;

typedef std::vector<StdRestraintNames> StdRestraintNamesList;
typedef StdRestraintNamesList::iterator StdRestraintNamesListIter;
typedef StdRestraintNamesList::const_iterator StdRestraintNamesListConstIter;

typedef std::vector<StdRestraintAtoms> StdRestraintAtomsList;
typedef StdRestraintAtomsList::iterator StdRestraintAtomsListIter;
typedef StdRestraintAtomsList::const_iterator StdRestraintAtomsListConstIter;

} // namespace rxdock

#endif //_RBTNOERESTRAINT_H_
