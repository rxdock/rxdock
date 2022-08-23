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

// Derived class for pseudo-atoms.

#ifndef _RBTPSEUDOATOM_H_
#define _RBTPSEUDOATOM_H_

#include "rxdock/Atom.h"

#include <nlohmann/json.hpp>

using json = nlohmann::json;

namespace rxdock {

class PseudoAtom : public Atom {
public:
  ///////////////////////////////////////////////
  // Constructors / destructors

  // Constructor supplying atom list of constituent atoms
  PseudoAtom(const AtomList &atomList, int nAtomId = -1,
             std::string strAtomName = "PSEUDO");

  friend void to_json(json &j, const PseudoAtom &pseudoAtom);
  friend void from_json(const json &j, PseudoAtom &pseudoAtom);

  // Default destructor
  virtual ~PseudoAtom();

  ///////////////////////////////////////////////
  // Public accessor functions
  // Virtual function for dumping pseudoatom details to an output stream
  virtual std::ostream &Print(std::ostream &s) const;

  ///////////////////////////////////////////////
  // Public accessor functions
  // 2-D attributes

  // Derived class methods for returning the constituent atom list
  unsigned int GetNumAtoms() const { return m_atomList.size(); }
  AtomList GetAtomList() const { return m_atomList; }

  ///////////////////////////////////////////////
  // Public accessor functions
  // 3-D attributes

  // Coords
  // DM 8 Dec 1998 - override base Atom methods
  // virtual Coord GetCoords() const;
  // virtual Double GetX() const;
  // virtual Double GetY() const;
  // virtual Double GetZ() const;
  //
  // DM 11 Jul 2000 - base GetCoords method is now non-virtual again
  // Pseudoatoms must update their mean coords by calling UpdateCoords
  // each time the underlying atoms are moved
  void UpdateCoords();

  ///////////////////////////////////////////////

  // Other public methods

protected:
private:
  // Private methods
  PseudoAtom();                             // Disable default constructor
  PseudoAtom(const PseudoAtom &pseudoAtom); // Disable copy constructor
  PseudoAtom &
  operator=(const PseudoAtom &pseudoAtom); // Disable copy assignment

private:
  // Private data
  AtomList m_atomList; // List of smart pointers to constituent atoms
};

void to_json(json &j, const PseudoAtom &pseudoAtom);
void from_json(const json &j, PseudoAtom &pseudoAtom);

// Useful typedefs
typedef SmartPtr<PseudoAtom> PseudoAtomPtr;        // Smart pointer
typedef std::vector<PseudoAtomPtr> PseudoAtomList; // Vector of smart pointers
typedef PseudoAtomList::iterator PseudoAtomListIter;
typedef PseudoAtomList::const_iterator PseudoAtomListConstIter;

} // namespace rxdock

#endif //_RBTPSEUDOATOM_H_
