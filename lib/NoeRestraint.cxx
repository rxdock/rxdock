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

#include "NoeRestraint.h"

using namespace rxdock;

// Returns a string representing the NOE restraint type
std::string rxdock::ConvertNoeRestraintTypeToString(const eNoeType type) {
  switch (type) {
  case NOE_OR:
    return "OR";
  case NOE_MEAN:
    return "MEAN";
  case NOE_AND:
    return "AND";
  default:
    return "UNDEFINED";
  }
}

///////////////////////////////////////////////////////////////////////
// Simple struct for holding one end of an NOE restraint
// Has an atom name string list and a restraint type
bool NoeEndNames::isOK() const {
  return (!names.empty() && (type != NOE_UNDEFINED));
}

// Insertion operator for the above
std::ostream &rxdock::operator<<(std::ostream &s, const NoeEndNames &n) {
  if (n.names.size() == 1)
    s << n.names.front();
  else {
    s << ConvertNoeRestraintTypeToString(n.type);
    s << "( ";
    for (std::vector<std::string>::const_iterator iter = n.names.begin();
         iter != n.names.end(); iter++)
      s << *iter << " ";
    s << ")";
  }
  return s;
}

///////////////////////////////////////////////////////////////////////
// As above, but with a real atom list in place of the string list

// Constructor accepting an NoeEndNames and an atom list to match against
NoeEndAtoms::NoeEndAtoms(const NoeEndNames &n, const AtomList &atomList) {
  atoms = GetMatchingAtomList(atomList, n.names);
  type = n.type;
}
bool NoeEndAtoms::isOK() const {
  return (!atoms.empty() && (type != NOE_UNDEFINED));
}

// Insertion operator for the above
std::ostream &rxdock::operator<<(std::ostream &s, const NoeEndAtoms &n) {
  if (n.atoms.size() == 1)
    s << n.atoms.front()->GetFullAtomName();
  else {
    s << ConvertNoeRestraintTypeToString(n.type);
    s << "( ";
    if (n.atoms.size() <= 4) {
      for (AtomListConstIter iter = n.atoms.begin(); iter != n.atoms.end();
           iter++)
        s << (*iter)->GetFullAtomName() << " ";
    } else {
      s << n.atoms.front()->GetFullAtomName() << " ... "
        << n.atoms.back()->GetFullAtomName() << " <" << n.atoms.size()
        << " atoms> ";
    }
    s << ")";
  }
  return s;
}

///////////////////////////////////////////////////////////////////////
// Simple struct for holding a complete NOE restraint definition
// Has two NoeEndNames, and max distance

bool NoeRestraintNames::isOK() const { return from.isOK() && to.isOK(); }

// Insertion operator for the above
std::ostream &rxdock::operator<<(std::ostream &s,
                                 const NoeRestraintNames &noe) {
  s << noe.from << " - " << noe.to << ": dist<" << noe.maxDist;
  return s;
}

///////////////////////////////////////////////////////////////////////
// As above, but with real NoeEndAtoms in place of NoeEndNames

// Constructor accepting an NoeRestraintNames and an atom list to match
// against
NoeRestraintAtoms::NoeRestraintAtoms(const NoeRestraintNames &n,
                                     const AtomList &atomList)
    : from(n.from, atomList), to(n.to, atomList), maxDist(n.maxDist) {}

bool NoeRestraintAtoms::isOK() const { return from.isOK() && to.isOK(); }

// true if exactly one atom at each end of the NOE (unambiguous)
bool NoeRestraintAtoms::isSimple() const {
  return (from.atoms.size() == 1) && (to.atoms.size() == 1);
}

// Insertion operator for the above
std::ostream &rxdock::operator<<(std::ostream &s,
                                 const NoeRestraintAtoms &noe) {
  s << noe.from << " - " << noe.to << ": dist<" << noe.maxDist;
  return s;
}

///////////////////////////////////////////////////////////////////////
// Simple struct for holding a complete STD restraint definition
// Has one NoeEndNames, and max distance from receptor

bool StdRestraintNames::isOK() const { return from.isOK(); }

// Insertion operator for the above
std::ostream &rxdock::operator<<(std::ostream &s,
                                 const StdRestraintNames &std) {
  s << std.from << " - receptor: dist<" << std.maxDist;
  return s;
}

///////////////////////////////////////////////////////////////////////
// As above, but with real NoeEndAtoms in place of NoeEndNames

// Constructor accepting an StdRestraintNames and an atom list to match
// against
StdRestraintAtoms::StdRestraintAtoms(const StdRestraintNames &s,
                                     const AtomList &atomList)
    : from(s.from, atomList), maxDist(s.maxDist) {}

bool StdRestraintAtoms::isOK() const { return from.isOK(); }

// true if exactly one atom at in the STD definition (unambiguous)
bool StdRestraintAtoms::isSimple() const { return (from.atoms.size() == 1); }

// Insertion operator for the above
std::ostream &rxdock::operator<<(std::ostream &s,
                                 const StdRestraintAtoms &std) {
  s << std.from << " - receptor: dist<" << std.maxDist;
  return s;
}
