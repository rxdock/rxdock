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

#include "rxdock/ChromDihedralRefData.h"
#include "rxdock/AtomFuncs.h"
#include "rxdock/Model.h"

#include <loguru.hpp>

#include <functional>

using namespace rxdock;

std::string ChromDihedralRefData::_CT = "ChromDihedralRefData";

ChromDihedralRefData::ChromDihedralRefData(BondPtr spBond,
                                           AtomList tetheredAtoms,
                                           double stepSize,
                                           ChromElement::eMode mode,
                                           double maxDihedral)
    : m_stepSize(stepSize), m_mode(mode), m_maxDihedral(maxDihedral) {
  Setup(spBond, tetheredAtoms);
  m_initialValue = GetModelValue();
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

ChromDihedralRefData::~ChromDihedralRefData() { _RBTOBJECTCOUNTER_DESTR_(_CT); }

double ChromDihedralRefData::GetModelValue() const {
  return BondDihedral(m_atom1, m_atom2, m_atom3, m_atom4);
}

void ChromDihedralRefData::SetModelValue(double dihedralAngle) {
  double delta = dihedralAngle - GetModelValue();
  // Only rotate if delta is non-zero
  if (std::fabs(delta) > 0.001) {
    // Coords of atom 1
    Coord coord1(m_atom2->GetCoords());
    // Vector along the bond between atom 1 and atom 2 (rotation axis)
    Vector bondVector(m_atom3->GetCoords() - coord1);
    Vector toOrigin = -coord1;
    Quat quat(bondVector, delta * M_PI / 180.0);
    for (AtomRListIter iter = m_rotAtoms.begin(); iter != m_rotAtoms.end();
         ++iter) {
      (*iter)->Translate(toOrigin);
      (*iter)->RotateUsingQuat(quat);
      (*iter)->Translate(coord1);
    }
  }
}

void ChromDihedralRefData::Setup(BondPtr spBond,
                                 const AtomList &tetheredAtoms) {
  Atom *pAtom2 = spBond->GetAtom1Ptr();
  Atom *pAtom3 = spBond->GetAtom2Ptr();
  // WARNING - we assume that the bond atoms are registered with a model
  Model *pModel = pAtom2->GetModelPtr();
  AtomList atomList = pModel->GetAtomList();
  BondList bondList = pModel->GetBondList();
  int nAtoms = atomList.size();
  int nTethered = tetheredAtoms.size();
  // The following lines get the bonded atom lists on each end of the rotable
  // bond, taking care not to include the atoms actually in the rotable bond
  AtomList bondedAtoms2 = GetAtomListWithPredicate(
      GetBondedAtomList(pAtom2),
      std::bind(std::not2(isAtomPtr_eq()), std::placeholders::_1, pAtom3));
  AtomList bondedAtoms3 = GetAtomListWithPredicate(
      GetBondedAtomList(pAtom3),
      std::bind(std::not2(isAtomPtr_eq()), std::placeholders::_1, pAtom2));
  // Assertion - check bonded atom lists are not empty
  Assert<Assertion>(!MUT_CHECK ||
                    !(bondedAtoms2.empty() || bondedAtoms3.empty()));

  ToSpin(spBond, atomList, bondList);
  pAtom2->SetSelectionFlag(false);
  pAtom3->SetSelectionFlag(false);
  // If we have selected over half the molecule to rotate then invert the
  // selection to minimise the number of atoms to rotate DM 2 Jul 2002 - in
  // tethered mode we want to ensure that the minimum number of *tethered* atoms
  // are rotated (preferably none). i.e. we rotate the free end of the bond,
  // even
  // if this is over half the molecule
  int nSelected = (nTethered == 0) ? GetNumSelectedAtomsInList(atomList)
                                   : GetNumSelectedAtomsInList(tetheredAtoms);
  int nHalf = (nTethered == 0) ? (nAtoms - 2) / 2 : (nTethered - 2) / 2;
  if (nSelected > nHalf) {
    LOG_F(1, "Over half the molecule selected: {} atoms", nSelected);
    InvertAtomSelectionFlags(atomList);
    pAtom2->SetSelectionFlag(false);
    pAtom3->SetSelectionFlag(false);
    m_atom1 = bondedAtoms3.front();
    m_atom2 = pAtom3;
    m_atom3 = pAtom2;
    m_atom4 = bondedAtoms2.front();
    LOG_F(1, "Inverted: {} atoms now selected",
          GetNumSelectedAtomsInList(atomList));
    LOG_F(1, "Dihedral spec: {}\t{}\t{}\t{}", bondedAtoms3.front()->GetName(),
          pAtom3->GetName(), pAtom2->GetName(),
          bondedAtoms2.front()->GetName());
  } else {
    m_atom1 = bondedAtoms2.front();
    m_atom2 = pAtom2;
    m_atom3 = pAtom3;
    m_atom4 = bondedAtoms3.front();
    LOG_F(1, "Dihedral spec: {}\t{}\t{}\t{}", bondedAtoms2.front()->GetName(),
          pAtom2->GetName(), pAtom3->GetName(),
          bondedAtoms3.front()->GetName());
  }

  // Store the smaller atom list (or free atom list in tethered mode) for this
  // rotable bond
  m_rotAtoms.clear();
  std::copy_if(atomList.begin(), atomList.end(), std::back_inserter(m_rotAtoms),
               isAtomSelected());
}
