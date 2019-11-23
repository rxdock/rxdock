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

#include "RbtAtom.h"

namespace Rbt {
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
} // namespace Rbt

///////////////////////////////////////////////////////////////////////
// Simple struct for holding one end of an NOE restraint
// Has an atom name string list and a restraint type
class RbtNoeEndNames {
public:
  RbtNoeEndNames() {}
  bool isOK() const;
  // Public data
  std::vector<std::string> names;
  Rbt::eNoeType type;
};

// Insertion operator for the above
std::ostream &operator<<(std::ostream &s, const RbtNoeEndNames &n);

///////////////////////////////////////////////////////////////////////
// As above, but with a real atom list in place of the string list
class RbtNoeEndAtoms {
public:
  // Default constructor
  RbtNoeEndAtoms() {}
  // Constructor accepting an RbtNoeEndNames and an atom list to match against
  RbtNoeEndAtoms(const RbtNoeEndNames &n, const RbtAtomList &atomList);
  bool isOK() const;
  // Public data
  RbtAtomList atoms;
  Rbt::eNoeType type;
};

// Insertion operator for the above
std::ostream &operator<<(std::ostream &s, const RbtNoeEndAtoms &n);

///////////////////////////////////////////////////////////////////////
// Simple struct for holding a complete NOE restraint definition
// Has two RbtNoeEndNames, and max distance
class RbtNoeRestraintNames {
public:
  RbtNoeRestraintNames() : maxDist(5.0) {}
  bool isOK() const;
  // Public data
  RbtNoeEndNames from;
  RbtNoeEndNames to;
  double maxDist;
};

// Insertion operator for the above
std::ostream &operator<<(std::ostream &s, const RbtNoeRestraintNames &noe);

///////////////////////////////////////////////////////////////////////
// As above, but with real RbtNoeEndAtoms in place of RbtNoeEndNames
class RbtNoeRestraintAtoms {
public:
  // Default constructor
  RbtNoeRestraintAtoms() : maxDist(5.0) {}
  // Constructor accepting an RbtNoeRestraintNames and an atom list to match
  // against
  RbtNoeRestraintAtoms(const RbtNoeRestraintNames &n,
                       const RbtAtomList &atomList);
  bool isOK() const; // true if at least one atom at each end of the NOE and
                     // NOE type is not undefined
  bool isSimple()
      const; // true if exactly one atom at each end of the NOE (unambiguous)
  // Public data
  RbtNoeEndAtoms from;
  RbtNoeEndAtoms to;
  // RbtDouble minDist;
  double maxDist;
};

// Insertion operator for the above
std::ostream &operator<<(std::ostream &s, const RbtNoeRestraintAtoms &noe);

///////////////////////////////////////////////////////////////////////
// Simple struct for holding a complete STD restraint definition
// Has one RbtNoeEndNames, and max distance (from receptor)
class RbtStdRestraintNames {
public:
  RbtStdRestraintNames() : maxDist(5.0) {}
  bool isOK() const;
  // Public data
  RbtNoeEndNames from;
  double maxDist;
};

// Insertion operator for the above
std::ostream &operator<<(std::ostream &s, const RbtStdRestraintNames &std);

///////////////////////////////////////////////////////////////////////
// As above, but with real RbtNoeEndAtoms in place of RbtNoeEndNames
class RbtStdRestraintAtoms {
public:
  // Default constructor
  RbtStdRestraintAtoms() : maxDist(5.0) {}
  // Constructor accepting an RbtStdRestraintNames and an atom list to match
  // against
  RbtStdRestraintAtoms(const RbtStdRestraintNames &n,
                       const RbtAtomList &atomList);
  bool isOK() const; // true if at least one atom in the definition and type
                     // is not undefined
  bool isSimple()
      const; // true if exactly one atom at in the STD definition (unambiguous)
  // Public data
  RbtNoeEndAtoms from;
  double maxDist;
};

// Insertion operator for the above
std::ostream &operator<<(std::ostream &s, const RbtStdRestraintAtoms &std);

// Useful typedefs
typedef std::vector<RbtNoeRestraintNames> RbtNoeRestraintNamesList;
typedef RbtNoeRestraintNamesList::iterator RbtNoeRestraintNamesListIter;
typedef RbtNoeRestraintNamesList::const_iterator
    RbtNoeRestraintNamesListConstIter;

typedef std::vector<RbtNoeRestraintAtoms> RbtNoeRestraintAtomsList;
typedef RbtNoeRestraintAtomsList::iterator RbtNoeRestraintAtomsListIter;
typedef RbtNoeRestraintAtomsList::const_iterator
    RbtNoeRestraintAtomsListConstIter;

typedef std::vector<RbtStdRestraintNames> RbtStdRestraintNamesList;
typedef RbtStdRestraintNamesList::iterator RbtStdRestraintNamesListIter;
typedef RbtStdRestraintNamesList::const_iterator
    RbtStdRestraintNamesListConstIter;

typedef std::vector<RbtStdRestraintAtoms> RbtStdRestraintAtomsList;
typedef RbtStdRestraintAtomsList::iterator RbtStdRestraintAtomsListIter;
typedef RbtStdRestraintAtomsList::const_iterator
    RbtStdRestraintAtomsListConstIter;

#endif //_RBTNOERESTRAINT_H_
