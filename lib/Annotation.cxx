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

#include "rxdock/Annotation.h"
#include "rxdock/PseudoAtom.h"

#include <functional>

using namespace rxdock;

////////////////////////////////////////
// Constructors/destructors
Annotation::Annotation(const Atom *pAtom1, const Atom *pAtom2, double dist,
                       double score)
    : m_pAtom1(pAtom1), m_pAtom2(pAtom2), m_dist(dist), m_score(score) {}

Annotation::~Annotation() {}

////////////////////////////////////////
// Public methods
////////////////
// Get
const Atom *Annotation::GetAtom1Ptr() const { return m_pAtom1; }
const Atom *Annotation::GetAtom2Ptr() const { return m_pAtom2; }
double Annotation::GetDistance() const { return m_dist; }
double Annotation::GetScore() const { return m_score; }

// Get the fully qualified (FQ) residue name for atom 2 (target atom)
std::string Annotation::GetFQResName() const {
  return m_pAtom2->GetSegmentName() + ":" + m_pAtom2->GetSubunitName() + "_" +
         m_pAtom2->GetSubunitId() + ":";
}

// Set
void Annotation::SetAtom1Ptr(const Atom *pAt1) { m_pAtom1 = pAt1; }
void Annotation::SetAtom2Ptr(const Atom *pAt2) { m_pAtom2 = pAt2; }
void Annotation::SetDistance(double d) { m_dist = d; }
void Annotation::SetScore(double s) { m_score = s; }

// Render annotation into string in rDock Viewer format
// Note: string does not contain the first field (annotation name)
// Full string required for SD file is:
// strName + "," + Render();
std::string Annotation::Render() const {
  std::ostringstream ostr;
  // Check if either atom is a pseudoatom
  // If so, arbitrarily write the first non-bridgehead atom ID
  const PseudoAtom *pseudo1 = dynamic_cast<const PseudoAtom *>(m_pAtom1);
  const PseudoAtom *pseudo2 = dynamic_cast<const PseudoAtom *>(m_pAtom2);
  if (pseudo1) {
    AtomList tmpList = pseudo1->GetAtomList();
    AtomList al1 =
        GetAtomListWithPredicate(tmpList, std::not1(isAtomBridgehead()));
    // Check for rare case of multifused rings, where all constituent ring atoms
    // are bridgeheads In this case, just use the first atom in the ring list
    AtomPtr spAt1 = (al1.empty()) ? tmpList.front() : al1.front();
    ostr << spAt1->GetAtomId() << ",";
  } else {
    ostr << m_pAtom1->GetAtomId() << ",";
  }
  if (pseudo2) {
    AtomList tmpList = pseudo2->GetAtomList();
    AtomList al2 =
        GetAtomListWithPredicate(tmpList, std::not1(isAtomBridgehead()));
    // Check for rare case of multifused rings, where all constituent ring atoms
    // are bridgeheads In this case, just use the first atom in the ring list
    AtomPtr spAt2 = (al2.empty()) ? tmpList.front() : al2.front();
    ostr << spAt2->GetAtomId() << ",";
  } else {
    ostr << m_pAtom2->GetAtomId() << ",";
  }
  ostr.setf(std::ios_base::fixed, std::ios_base::floatfield);
  ostr.precision(2);
  ostr << m_dist << "," << m_score;
  std::string retVal(ostr.str());
  return retVal;
}

// Operators
// Special meaning of operator+ for accumulating annotations by residue
void Annotation::operator+=(const Annotation &ann) {
  // Check FQ residue names match, if not then don't accumulate
  if (GetFQResName() != ann.GetFQResName()) {
    return;
  }
  // atom 1 and distance relate to the closest contact between the ligand and
  // the target residue, so update them here if necessary
  if (ann.GetDistance() < m_dist) {
    m_pAtom1 = ann.GetAtom1Ptr();
    m_dist = ann.GetDistance();
  }
  // Score represents the total score between the ligand and the target residue
  m_score += ann.GetScore();
  return;
}
