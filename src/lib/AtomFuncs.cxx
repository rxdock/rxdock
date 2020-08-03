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

#include "AtomFuncs.h"

#include <loguru.hpp>

using namespace rxdock;

// DM 31 Oct 2000
// Given a bond, determines if it is in a ring (cutdown version of ToSpin)
bool rxdock::FindCyclic(BondPtr spBond, AtomList &atomList,
                        BondList &bondList) {
  // Max no. of atoms to process before giving up and return bCyclic=false
  // DM 8 Nov 2000 - effectively remove the cap
  const int maxProc = 999999;

  AtomPtr spAtom1 = spBond->GetAtom1Ptr();
  AtomPtr spAtom2 = spBond->GetAtom2Ptr();

  SetAtomSelectionFlagsInList(atomList,
                              false); // Clear the atom selection flags
  SetBondSelectionFlagsInList(bondList,
                              false); // Clear the bond selection flags

  // Temporary atom list containing atoms to be processed
  // Note: this is a true list (not a vector) as we will be making numerous
  // insertions and deletions
  AtomTrueList pendingAtomList;

  // Seed the pending list with atom 2
  spAtom2->SetSelectionFlag(true);
  spBond->SetSelectionFlag(true);
  pendingAtomList.push_back(spAtom2);

  int nProc = 0;

  // While we still have atoms to process
  while (!pendingAtomList.empty() && (nProc < maxProc)) {
    // Take the last atom from the list and remove it
    AtomPtr spA1 = pendingAtomList.back();
    pendingAtomList.pop_back();
    nProc++;

    // Get the bonds this atom is in
    const BondMap &bondMap = spA1->GetBondMap();

    // Loop round each bond
    for (BondMapConstIter bIter = bondMap.begin(); bIter != bondMap.end();
         bIter++) {
      Bond *pBnd = (*bIter).first; // Remember not to use smart pointers!
      // If we haven't considered this bond before
      if (!pBnd->GetSelectionFlag()) {
        // If this atom is first in the bond (bIter->second == true), get atom 2
        // ptr, else get atom 1 ptr
        AtomPtr spA2 =
            ((*bIter).second) ? pBnd->GetAtom2Ptr() : pBnd->GetAtom1Ptr();
        if (spA2 == spAtom1) {
          LOG_F(1, "{} is cyclic", spAtom1->GetFullAtomName());
          return true;
        } else if (!spA2->GetSelectionFlag()) {
          // Mark the atom and bond as selected and add the atom to the pending
          // list
          spA2->SetSelectionFlag(true);
          pBnd->SetSelectionFlag(true);
          pendingAtomList.push_back(spA2);
        }
      }
    }
  }
  return false;
}

// DM 4 Dec 1998
// Given a bond, set the selection flags for all atoms which are connected to
// atom 2 of the bond Returns true if bond is in a ring (i.e. if atom 1's flag
// gets set also) DM 8 Feb 2000 - standalone version (formerly only available as
// Model method)
bool rxdock::ToSpin(BondPtr spBond, AtomList &atomList, BondList &bondList) {
  AtomPtr spAtom1 = spBond->GetAtom1Ptr();
  AtomPtr spAtom2 = spBond->GetAtom2Ptr();

  // Check if both atoms are actually in the model
  // if ( (spAtom1->GetModelPtr() != this) || (spAtom2->GetModelPtr() != this) )
  //  return false;

  SetAtomSelectionFlagsInList(atomList,
                              false); // Clear the atom selection flags
  SetBondSelectionFlagsInList(bondList,
                              false); // Clear the bond selection flags
  bool bIsCyclic(false);

  // Temporary atom list containing atoms to be processed
  // Note: this is a true list (not a vector) as we will be making numerous
  // insertions and deletions
  AtomTrueList pendingAtomList;

  // Seed the pending list with atom 2
  spAtom2->SetSelectionFlag(true);
  spBond->SetSelectionFlag(true);
  pendingAtomList.push_back(spAtom2);

  LOG_F(1, "ToSpin: Bond ID={}", spBond->GetBondId());

  // While we still have atoms to process
  while (!pendingAtomList.empty()) {
    // Take the last atom from the list and remove it
    AtomPtr spA1 = pendingAtomList.back();
    pendingAtomList.pop_back();
    LOG_F(1, "ToSpin: Checking atom {}", spA1->GetAtomId());

    // Get the bonds this atom is in
    const BondMap &bondMap = spA1->GetBondMap();

    // Loop round each bond
    for (BondMapConstIter bIter = bondMap.begin(); bIter != bondMap.end();
         bIter++) {
      Bond *pBnd = (*bIter).first; // Remember not to use smart pointers!
      // If we haven't considered this bond before
      if (!pBnd->GetSelectionFlag()) {
        // If this atom is first in the bond (bIter->second == true), get atom 2
        // ptr, else get atom 1 ptr
        AtomPtr spA2 =
            ((*bIter).second) ? pBnd->GetAtom2Ptr() : pBnd->GetAtom1Ptr();
        LOG_F(1, "ToSpin: Checking bond {} to atom {}", pBnd->GetBondId(),
              spA2->GetAtomId());
        // If we've got back to atom 1 of the bond passed into ToSpin, then the
        // bond must be in a ring so set the cyclic flag but don't consider this
        // atom any further (we don't want to set the selection flag otherwise
        // we'd end up selecting the entire molecule).
        if (spA2 == spAtom1) {
          bIsCyclic = true;
          LOG_F(1, "ToSpin: We've hit {} so compound must be cyclic",
                spAtom1->GetAtomId());
        } else if (!spA2->GetSelectionFlag()) {
          // Mark the atom and bond as selected and add the atom to the pending
          // list
          spA2->SetSelectionFlag(true);
          pBnd->SetSelectionFlag(true);
          pendingAtomList.push_back(spA2);
          LOG_F(1, "ToSpin: Adding atom {} to the pending list",
                spA2->GetAtomId());
        }
      }
    }
  }
  return bIsCyclic;
}

// DM 7 Dec 1998
// Set the atom and bond cyclic flags for all atoms and bonds in the model
// DM 8 Feb 2000 - standalone version (formerly only available as Model
// method)
void rxdock::SetAtomAndBondCyclicFlags(AtomList &atomList, BondList &bondList) {
  // Clear the flags
  SetAtomCyclicFlagsInList(atomList, false);
  SetBondCyclicFlagsInList(bondList, false);

  // A bond is cyclic if FindCyclic returns true
  // An atom is cyclic if it is involved in a cyclic bond
  // So loop around each bond, call FindCyclic and set the flags for the bond,
  // atom 1 and atom 2 accordingly
  for (BondListIter iter = bondList.begin(); iter != bondList.end(); iter++) {
    BondPtr spBond(*iter);
    if (FindCyclic(spBond, atomList, bondList)) {
      spBond->SetCyclicFlag(true);
      spBond->GetAtom1Ptr()->SetCyclicFlag(true);
      spBond->GetAtom2Ptr()->SetCyclicFlag(true);
    }
  }
}

// Find the smallest ring containing the given atom
// Assumes SetAtomAndBondCyclicFlags has already been called
// Note: if the atom is a member of two equally sized rings, only one will be
// returned 30 Oct 2000 - Standalone version
AtomList rxdock::FindRing(AtomPtr spAtom, BondList &bondList) {
  // Check that atom is cyclic
  if (!spAtom->GetCyclicFlag())
    return AtomList();

  SetBondSelectionFlagsInList(bondList,
                              false); // Clear the bond selection flags

  // Get the cyclic bonds this atom is in
  BondMap cyclicBondMap = spAtom->GetCyclicBondMap();
  BondMapIter bIter = cyclicBondMap.begin();
  (*bIter).first->SetSelectionFlag(true);
  AtomPtr spAtom2 = GetBondedAtomPtr(*bIter);
  // Seed the partial ring lists with the first two atoms
  AtomList firstRing;
  firstRing.push_back(spAtom);
  firstRing.push_back(spAtom2);

  // This is a list of atom lists (actually a vector of atom vectors)
  AtomListList partialRings;
  partialRings.push_back(firstRing);

  // A bit risky, but loop forever and trust that we will find the seed atom
  // eventually After all, it is in a ring! Whichever partial ring makes it
  // first must be the smallest ring
  for (;;) {
    // First replicate any partial rings whose head atom is at a ring branch
    // point (defined by having >1 unselected cyclic bond)
    unsigned int nRings = partialRings.size();
    LOG_F(1, "FindRing: {} rings prior to forking", nRings);

    // Use indices rather than iterators, as we may be increasing the number of
    // partial rings in the ring list (Iterators are unstable if the container
    // increases in size inside the loop)
    for (unsigned int ringIter = 0; ringIter < nRings; ringIter++) {
      // Get the head atom
      AtomPtr spHeadAtom = (partialRings[ringIter]).back();
      // Get the cyclic bonds this atom is in
      BondMap cyclicBondMap = spHeadAtom->GetCyclicBondMap();
      // Count how many are unselected
      int nUnSelected(0);
      // Loop round each bond
      for (BondMapIter bIter = cyclicBondMap.begin();
           bIter != cyclicBondMap.end(); bIter++) {
        if (!(*bIter).first->GetSelectionFlag())
          nUnSelected++;
      }
      LOG_F(1, "FindRing: #Unselected cyclic bonds={} for {} at atom {}",
            nUnSelected, ringIter, spHeadAtom->GetAtomId());

      // Add a new copy of the partial ring for each unselected cyclic bond from
      // this atom above the normal one
      for (int forkIter = 0; forkIter < (nUnSelected - 1); forkIter++) {
        LOG_F(1, "FindRing: Forking ring path {} at atom {}", ringIter,
              spHeadAtom->GetAtomId());
        partialRings.push_back(partialRings[ringIter]);
      }
    }

    // Now grow each partial ring by one atom from the head atom
    for (AtomListListIter rIter = partialRings.begin();
         rIter != partialRings.end(); rIter++) {
      AtomPtr spHeadAtom = (*rIter).back();
      // Get the cyclic bonds this atom is in
      BondMap cyclicBondMap = spHeadAtom->GetCyclicBondMap();
      // Loop round each bond
      for (BondMapIter bIter = cyclicBondMap.begin();
           bIter != cyclicBondMap.end(); bIter++) {
        Bond *pBnd = (*bIter).first; // Remember not to use smart pointers!
        // Check we haven't considered this bond already
        if (!pBnd->GetSelectionFlag()) {
          // If this atom is first in the bond (bIter->second == true), get atom
          // 2 ptr, else get atom 1 ptr
          AtomPtr spA2 = GetBondedAtomPtr(*bIter);
          // We've hit atom 1 so this must be the smallest ring
          if (spA2 == spAtom) {
            LOG_F(1, "FindRing: Found seed atom {}", spAtom->GetAtomId());
            return *rIter;
          } else {
            LOG_F(1, "FindRing: Adding atom {} to {}", spA2->GetAtomId(),
                  spHeadAtom->GetAtomId());
            (*rIter).push_back(spA2);
            pBnd->SetSelectionFlag(true);
            break; // Break because we only want to add one atom to each head
                   // atom
          }
        }
      }
    }
  }
}

// 30 Oct 2000 (DM) - Find all rings, standalone version
void rxdock::FindRings(AtomList &atomList, BondList &bondList,
                       AtomListList &ringList) {
  ringList.clear();
  // Determine the cyclic atoms and bonds
  SetAtomAndBondCyclicFlags(atomList, bondList);

  // First set the selection flags to the inverse of the cyclic flag
  // so we initially have all the cyclic atoms unselected
  for (AtomListIter iter = atomList.begin(); iter != atomList.end(); iter++)
    (*iter)->SetSelectionFlag(!(*iter)->GetCyclicFlag());

  // Loop over all atoms and use each unselected atom as the seed for FindRing
  for (AtomListIter iter = atomList.begin(); iter != atomList.end(); iter++) {
    if (!(*iter)->GetSelectionFlag()) {
      AtomList ringAtoms = FindRing(*iter, bondList);
      if (!ringAtoms.empty()) {
        // Now we've found a ring, set the selection flags to true for each ring
        // member
        SetAtomSelectionFlagsInList(ringAtoms, true);
        // Store in our list of rings
        ringList.push_back(ringAtoms);
      }
    }
  }
}
