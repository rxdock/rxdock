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

#include "ModelMutator.h"
#include "AtomFuncs.h"
#include "Model.h"

#include <functional>

using namespace rxdock;

#ifdef _DEBUG
namespace Tmp {
void PrintAtomList(AtomRListConstIter i1, AtomRListConstIter i2) {
  for (AtomRListConstIter iter = i1; iter != i2; iter++) {
    std::cout << "\t" << (*iter)->GetFullAtomName() << " ("
              << (*iter)->GetAtomId() << ")";
  }
  std::cout << std::endl;
}
} // namespace Tmp
#endif //_DEBUG

////////////////////////////////////////
// Constructors/destructors
ModelMutator::ModelMutator(Model *pModel, const BondList &rotBonds,
                           const AtomList &tetheredAtoms)
    : m_pModel(pModel), m_rotBonds(rotBonds), m_tetheredAtoms(tetheredAtoms) {
  if (m_pModel != nullptr) {
    Setup();
  }
  _RBTOBJECTCOUNTER_CONSTR_("ModelMutator");
}

ModelMutator::~ModelMutator() { _RBTOBJECTCOUNTER_DESTR_("ModelMutator"); }

////////////////////////////////////////
// Public methods
////////////////

const AtomRListList &ModelMutator::GetFlexIntns() const { return m_flexIntns; }

const AtomRListList &ModelMutator::GetFlexAtoms() const { return m_rotAtoms; }

BondList ModelMutator::GetFlexBonds() const { return m_rotBonds; }

////////////////////////////////////////
// Private methods
////////////////
void ModelMutator::Setup() {
  m_dih1Atoms.clear();
  m_dih2Atoms.clear();
  m_dih3Atoms.clear();
  m_dih4Atoms.clear();
  m_rotAtoms.clear();
  m_flexIntns.clear();

  if (m_pModel == nullptr)
    return;

  int nAtoms = m_pModel->GetNumAtoms();
  int nTethered = m_tetheredAtoms.size();
  m_flexIntns = AtomRListList(nAtoms, AtomRList());

  // Loop over each rotable bond
  for (BondListConstIter iter = m_rotBonds.begin(); iter != m_rotBonds.end();
       iter++) {
    Atom *pAtom2 = (*iter)->GetAtom1Ptr();
    Atom *pAtom3 = (*iter)->GetAtom2Ptr();
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

    ToSpin(*iter, m_pModel->m_atomList, m_pModel->m_bondList);
    pAtom2->SetSelectionFlag(false);
    pAtom3->SetSelectionFlag(false);
    // If we have selected over half the molecule to rotate then invert the
    // selection to minimise the number of atoms to rotate DM 2 Jul 2002 - in
    // tethered mode we want to ensure that the minimum number of *tethered*
    // atoms are rotated (preferably none). i.e. we rotate the free end of the
    // bond, even if this is over half the molecule
    int nSelected = (nTethered == 0)
                        ? GetNumSelectedAtomsInList(m_pModel->m_atomList)
                        : GetNumSelectedAtomsInList(m_tetheredAtoms);
    int nHalf = (nTethered == 0) ? (nAtoms - 2) / 2 : (nTethered - 2) / 2;
    if (nSelected > nHalf) {
      // std::cout << "Over half the molecule selected: " << nSelected << "
      // atoms" << std::endl;
      InvertAtomSelectionFlags(m_pModel->m_atomList);
      pAtom2->SetSelectionFlag(false);
      pAtom3->SetSelectionFlag(false);
      m_dih1Atoms.push_back(bondedAtoms3.front());
      m_dih2Atoms.push_back(pAtom3);
      m_dih3Atoms.push_back(pAtom2);
      m_dih4Atoms.push_back(bondedAtoms2.front());
      // std::cout << "Inverted: " <<
      // GetNumSelectedAtoms(m_pModel->m_atomList)
      // << " atoms now selected" << std::endl; std::cout << "Dihedral spec: "
      // << bondedAtoms3.front()->GetName() << "\t" << pAtom3->GetName()
      // <<
      // "\t"
      //	   << pAtom2->GetName() << "\t" <<
      // bondedAtoms2.front()->GetName() << std::endl;
    } else {
      m_dih1Atoms.push_back(bondedAtoms2.front());
      m_dih2Atoms.push_back(pAtom2);
      m_dih3Atoms.push_back(pAtom3);
      m_dih4Atoms.push_back(bondedAtoms3.front());
      // std::cout << "Dihedral spec: " << bondedAtoms2.front()->GetName()
      // <<
      // "\t" << pAtom2->GetName() << "\t"
      //	   << pAtom3->GetName() << "\t" <<
      // bondedAtoms3.front()->GetName() << std::endl;
    }

    // Store the smaller atom list (or free atom list in tethered mode) for this
    // rotable bond
    AtomRList sList;
    std::copy_if(m_pModel->m_atomList.begin(), m_pModel->m_atomList.end(),
                 std::back_inserter(sList), isAtomSelected());
    m_rotAtoms.push_back(sList);

    // DM 25 Apr 2002
    // Work out the flexible interactions across the bond
    // Method is as follows.
    // sList = selected atoms on one side of the bond
    // uList = unselected atoms on other side of the bond
    // By definition, all atom pairs from {sList}-{uList} are flexible
    // interactions but we need to avoid counting interaction pairs more than
    // once foreach rotable bond {
    //  Determine sList, uList {
    //  foreach sAtom in sList {
    //    foreach uAtom in uList {
    //      append uAtom to sAtom's interaction list
    //      append sAtom to uAtom's interaction list
    //    }
    //  }
    //}
    // foreach atom in model {
    //  sort interaction list for atom
    //  remove duplicates from interaction list
    //}
    InvertAtomSelectionFlags(m_pModel->m_atomList);
    pAtom2->SetSelectionFlag(false);
    pAtom3->SetSelectionFlag(false);
    AtomRList uList;
    std::copy_if(m_pModel->m_atomList.begin(), m_pModel->m_atomList.end(),
                 std::back_inserter(uList), isAtomSelected());
    // Sort the lists of selected and unselected atoms by atom ID
    std::sort(sList.begin(), sList.end(), AtomPtrCmp_AtomId());
    std::sort(uList.begin(), uList.end(), AtomPtrCmp_AtomId());

    // Loop over selected atom list
    // and append the unselected atom list to the interaction list for each atom
    for (AtomRListConstIter sIter = sList.begin(); sIter != sList.end();
         sIter++) {
      unsigned int id = (*sIter)->GetAtomId() - 1;
      // Assertion - check id is within range
      Assert<Assertion>(!MUT_CHECK || (id < m_flexIntns.size()));
      AtomRListListIter lIter = m_flexIntns.begin() + id;
      std::copy(uList.begin(), uList.end(), std::back_inserter(*lIter));
    }

    // Ditto for unselected atom list
    for (AtomRListConstIter uIter = uList.begin(); uIter != uList.end();
         uIter++) {
      unsigned int id = (*uIter)->GetAtomId() - 1;
      // Assertion - check id is within range
      Assert<Assertion>(!MUT_CHECK || (id < m_flexIntns.size()));
      AtomRListListIter lIter = m_flexIntns.begin() + id;
      std::copy(sList.begin(), sList.end(), std::back_inserter(*lIter));
    }
  }
  // Now clean up the interaction lists for each atom to remove duplicates
  for (AtomRListListIter lIter = m_flexIntns.begin();
       lIter != m_flexIntns.end(); lIter++) {
    if (!(*lIter).empty()) {
      std::sort((*lIter).begin(), (*lIter).end(), AtomPtrCmp_AtomId());
      // std::unique will shuffle all the duplicate elements to the end of the
      // vector and return an iterator to the end of the unique elements
      AtomRListIter uniqIter = std::unique((*lIter).begin(), (*lIter).end());
#ifdef _DEBUG
      int id = lIter - m_flexIntns.begin();
      Atom *pAtom = m_pModel->m_atomList[id];
      std::cout << "Atom " << pAtom->GetFullAtomName() << " (" << id + 1
                << "): with dups=" << (*lIter).size();
#endif //_DEBUG
      // Remove duplicates
      (*lIter).erase(uniqIter, (*lIter).end());
#ifdef _DEBUG
      std::cout << "; unique=" << (*lIter).size() << std::endl;
      Tmp::PrintAtomList((*lIter).begin(), (*lIter).end());
#endif //_DEBUG
    }
  }
}
