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

#include "rxdock/PseudoAtom.h"

using namespace rxdock;

///////////////////////////////////////////////
// Constructors / destructors

// Constructor supplying atom list of constituent atoms
// Note: the base Atom constructor sets up default attributes for the
// pseudoatom only (the constituent atoms are already constructed so we don't
// need to worry about them) Use Atomic No = -1 to signify a pseudoatom
PseudoAtom::PseudoAtom(const AtomList &atomList, int nAtomId,
                       std::string strAtomName)
    : Atom(nAtomId, -1, strAtomName) {
  // Check atom list isn't empty (should be an assert)
  if (atomList.empty())
    throw BadArgument(_WHERE_,
                      "Null atom list passed to PseudoAtom constructor");

  // Store the atom list
  m_atomList = atomList;

  // Base the pseudoatom attributes on the first atom in the list
  // Note: this is a write-once action; if the attributes of the first
  // atom change, the pseudoatom attributes will not.
  // Assumes that all constituent atoms are in the same subunit and segment
  AtomPtr spAtom1(atomList.front()); // 1st atom in list
  AtomPtr spAtomN(atomList.back());  // Last atom in list
  // Pseudo-atom name is i the format [C1..N4]
  strAtomName += "[" + spAtom1->GetName() + ".." + spAtomN->GetName() + "]";
  SetAtomName(strAtomName);
  SetSubunitId(spAtom1->GetSubunitId());
  SetSubunitName(spAtom1->GetSubunitName());
  SetSegmentName(spAtom1->GetSegmentName());
  UpdateCoords();
  _RBTOBJECTCOUNTER_CONSTR_("PseudoAtom");
}

// Default destructor
PseudoAtom::~PseudoAtom() { _RBTOBJECTCOUNTER_DESTR_("PseudoAtom"); }

// Virtual function for dumping atom details to an output stream
// Derived classes (e.g. pseudoatom) can override if required
std::ostream &PseudoAtom::Print(std::ostream &s) const {
  Atom::Print(s); // First call base class Print on the pseudo-atom itself
  s << std::endl
    << "Pseudo-atom with " << m_atomList.size()
    << " constituent atoms:" << std::endl;
  // Now dump the atom details for all the constituent atoms
  for (AtomListConstIter iter = m_atomList.begin(); iter != m_atomList.end();
       iter++)
    s << (**iter) << std::endl;
  return s;
}

///////////////////////////////////////////////
// Public accessor functions
// 3-D attributes

// Coords
// DM 11 Jul 2000 - base GetCoords method is now non-virtual again
// Pseudoatoms must update their mean coords by calling UpdateCoords
// each time the underlying atoms are moved
void PseudoAtom::UpdateCoords() {
  Coord coord;
  unsigned int nSize = m_atomList.size();
  if (nSize > 0) {
    for (AtomListConstIter iter = m_atomList.begin(); iter != m_atomList.end();
         iter++) {
      coord += (*iter)->GetCoords();
    }
    coord /= nSize;
  }
  SetCoords(coord);
}

void rxdock::to_json(json &j, const PseudoAtom &pseudoAtom) {
  json atomList;
  for (const auto &aIter : pseudoAtom.m_atomList) {
    json atom = *aIter;
    atomList.push_back(atom);
  }
  j = json{{"atoms", atomList}};
}

void rxdock::from_json(const json &j, PseudoAtom &pseudoAtom) {
  for (auto &atom : j.at("atoms")) {
    AtomPtr spAtom = AtomPtr(new Atom(atom));
    pseudoAtom.m_atomList.push_back(spAtom);
  }
}