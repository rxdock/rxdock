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

#include "NmrSF.h"
#include "NmrRestraintFileSource.h"
#include "WorkSpace.h"

#include <functional>

using namespace rxdock;

// Static data members
std::string NmrSF::_CT("NmrSF");
std::string NmrSF::_FILENAME("FILENAME");
std::string NmrSF::_QUADRATIC("QUADRATIC");

// NB - Virtual base class constructor (BaseSF) gets called first,
// implicit constructors for BaseInterSF is called second
NmrSF::NmrSF(const std::string &strName)
    : BaseSF(_CT, strName), m_bQuadratic(true) {
  // Add parameters
  AddParameter(_FILENAME, std::string("default.noe"));
  AddParameter(_QUADRATIC, m_bQuadratic);
  SetRange(6.0);
  SetTrace(1);
#ifdef _DEBUG
  std::cout << _CT << " parameterised constructor" << std::endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

NmrSF::~NmrSF() {
#ifdef _DEBUG
  std::cout << _CT << " destructor" << std::endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

void NmrSF::SetupReceptor() {
  m_spGrid = NonBondedGridPtr();
  if (GetReceptor().Null())
    return;

  int iTrace = GetTrace();

  // Create indexed grid for STD restraints penalty function
  m_spGrid = CreateNonBondedGrid();
  // Fixed distance range used for indexing the relevant receptor atoms
  double range = GetRange() + GetMaxError();
  DockingSitePtr spDS = GetWorkSpace()->GetDockingSite();

  // Find the non-polar hydrogens (nph) in the receptor
  // If we are using an implicit hydrogen model, need to include extended
  // carbons (xc) also
  AtomList recepAtomList = GetReceptor()->GetAtomList();
  AtomList nphList = GetAtomListWithPredicate(recepAtomList, isAtomicNo_eq(1));
  nphList = GetAtomListWithPredicate(nphList, std::not1(isAtomHBondDonor()));
  AtomList xcList = GetAtomListWithPredicate(recepAtomList, isAtomExtended());
  std::copy(xcList.begin(), xcList.end(), std::back_inserter(nphList));

  // Now limit the list to those atoms in the vicinity of the docking site
  nphList = spDS->GetAtomList(nphList, 0.0, GetCorrectedRange());

  if (iTrace > 0) {
    std::cout << _CT << "::SetupReceptor() - " << nphList.size()
              << " nonpolar Hs and extended carbons near docking site"
              << std::endl;
    for (AtomListConstIter iter = nphList.begin(); iter != nphList.end();
         iter++) {
      std::cout << (**iter) << std::endl;
    }
  }

  for (AtomListConstIter iter = nphList.begin(); iter != nphList.end();
       iter++) {
    m_spGrid->SetAtomLists(*iter, range);
  }
}

void NmrSF::SetupLigand() {
  m_ligAtomList.clear();
  if (GetLigand().Null())
    return;

  m_ligAtomList = GetLigand()->GetAtomList();
}

void NmrSF::SetupScore() {
  m_noeList.clear();
  m_stdList.clear();
  if (GetLigand().Null() || GetReceptor().Null())
    return;

  int iTrace = GetTrace();

  // Create a restraint filesource with the appropriate filename
  std::string strRestrFile =
      GetDataFileName("data/receptors", GetParameter(_FILENAME));
  NmrRestraintFileSourcePtr spRestrSource(
      new NmrRestraintFileSource(strRestrFile));
  std::string strName = GetName();

  // Now convert from restraint names to an actual restraint atoms list by
  // matching against the receptor and ligand atom lists, and setup the noe
  // restraint interactions
  AtomList atomList = GetReceptor()->GetAtomList();
  std::copy(m_ligAtomList.begin(), m_ligAtomList.end(),
            std::back_inserter(atomList));

  // (1) NOE Restraints
  // Get the noe restraint names list from the restraint file source
  NoeRestraintNamesList noeNamesList = spRestrSource->GetNoeRestraintList();
  for (NoeRestraintNamesListIter rIter = noeNamesList.begin();
       rIter != noeNamesList.end(); rIter++) {
    // Create the NOE restraint - constructor matches the atoms in atomList with
    // the names in *rIter
    NoeRestraintAtoms noe(*rIter, atomList);
    if (iTrace > 0) {
      std::cout << std::endl << *rIter << std::endl << noe << std::endl;
    }
    // Only store if NOE is OK - i.e. has found the restraint names in the atom
    // list
    if (noe.isOK()) {
      m_noeList.push_back(noe);
    } else {
      std::cout
          << "** WARNING - unable to match NOE restraint names to atom list"
          << std::endl;
      std::cout << noe << " not added" << std::endl;
    }
  }

  // (2) STD Restraints
  // Get the std restraint names list from the restraint file source
  StdRestraintNamesList stdNamesList = spRestrSource->GetStdRestraintList();
  for (StdRestraintNamesListIter rIter = stdNamesList.begin();
       rIter != stdNamesList.end(); rIter++) {
    // Create the STD restraint - constructor matches the atoms in m_ligAtomList
    // with the names in *rIter Note for STD restraints it only makes senses to
    // check the names against the ligand atoms
    StdRestraintAtoms std(*rIter, m_ligAtomList);
    if (iTrace > 0) {
      std::cout << std::endl << *rIter << std::endl << std << std::endl;
    }
    // Only store if STD is OK - i.e. has found the restraint names in the atom
    // list and maxDist is within the range of the scoring function
    if (std.maxDist > GetRange()) {
      std::cout << "** WARNING - maxDist is beyond the range of the scoring "
                   "function index"
                << std::endl;
      std::cout << "** Increase the RANGE parameter in " << GetFullName()
                << std::endl;
      std::cout << std << " not added" << std::endl;
    } else if (std.isOK()) {
      m_stdList.push_back(std);
    } else {
      std::cout
          << "** WARNING - unable to match STD restraint names to atom list"
          << std::endl;
      std::cout << std << " not added" << std::endl;
    }
  }
}

double NmrSF::RawScore() const {
  double score(0.0);
  // QUADRATIC POTENTIAL
  if (m_bQuadratic) {
    for (NoeRestraintAtomsListConstIter iter = m_noeList.begin();
         iter != m_noeList.end(); iter++) {
      double r = NoeDistance(*iter);
      double dr = r - (*iter).maxDist; // delta(R)
      double s = (dr > 0.0) ? dr * dr : 0.0;
      // std::cout << "(NOE) R,RMAX,DR,S=" << r << "," << (*iter).maxDist << ","
      // << dr << "," << s << std::endl;
      score += s;
    }
    for (StdRestraintAtomsListConstIter iter = m_stdList.begin();
         iter != m_stdList.end(); iter++) {
      double r = StdDistance(*iter);
      double dr = r - (*iter).maxDist; // delta(R)
      double s = (dr > 0.0) ? dr * dr : 0.0;
      // std::cout << "(STD) R,RMAX,DR,S=" << r << "," << (*iter).maxDist << ","
      // << dr << "," << s << std::endl;
      score += s;
    }
  }
  // LINEAR POTENTIAL
  else {
    for (NoeRestraintAtomsListConstIter iter = m_noeList.begin();
         iter != m_noeList.end(); iter++) {
      double r = NoeDistance(*iter);
      double dr = r - (*iter).maxDist; // delta(R)
      double s = (dr > 0.0) ? dr : 0.0;
      // std::cout << "(NOE) R,RMAX,DR,S=" << r << "," << (*iter).maxDist << ","
      // << dr << "," << s << std::endl;
      score += s;
    }
    for (StdRestraintAtomsListConstIter iter = m_stdList.begin();
         iter != m_stdList.end(); iter++) {
      double r = StdDistance(*iter);
      double dr = r - (*iter).maxDist; // delta(R)
      double s = (dr > 0.0) ? dr : 0.0;
      // std::cout << "(STD) R,RMAX,DR,S=" << r << "," << (*iter).maxDist << ","
      // << dr << "," << s << std::endl;
      score += s;
    }
  }
  return score;
}

double NmrSF::NoeDistance(const NoeRestraintAtoms &noe) const {
  // Simple, unambiguous restraint - just return the interatomic distance
  if (noe.isSimple()) {
    return BondLength(noe.from.atoms.front(), noe.to.atoms.front());
  }
  // More complex cases
  // Compile the list of coords at each end
  // If the restraint end is defined of type MEAN, just store the center of mass
  // of the atoms in the list If the restraint end is defined of type OR or AND,
  // store all coords
  else {
    CoordList fromCoords;
    if (noe.from.type == NOE_MEAN)
      fromCoords.push_back(GetCenterOfAtomicMass(noe.from.atoms));
    else
      fromCoords = GetCoordList(noe.from.atoms);

    CoordList toCoords;
    if (noe.to.type == NOE_MEAN)
      toCoords.push_back(GetCenterOfAtomicMass(noe.to.atoms));
    else
      toCoords = GetCoordList(noe.to.atoms);

    double dist_sq(0.0); // Keep track of minimum distance**2
    // Iterate over coords in each list and return the appropriate distance**2
    // between any of them
    for (CoordListConstIter fIter = fromCoords.begin();
         fIter != fromCoords.end(); fIter++) {
      // dist1_sq is the appropriate distance**2 between the current "from"
      // coord and all the coords in the "to" list. i.e. if to.type==NOE_AND,
      // dist1_sq is the max distance**2 to any atom in the "to" list else
      // dist1_sq is the min distance**2 to any atom in the "to" list
      double dist1_sq(0.0);
      for (CoordListConstIter tIter = toCoords.begin(); tIter != toCoords.end();
           tIter++) {
        double r12_sq = Length2(*fIter, *tIter);
        dist1_sq = (tIter == toCoords.begin())
                       ? r12_sq
                       : (noe.to.type == NOE_AND) ? std::max(dist1_sq, r12_sq)
                                                  : std::min(dist1_sq, r12_sq);
      }
      // dist_sq is the appropriate overall distance**2 for all the calculated
      // dist1_sq's i.e. if from.type==NOE_AND, dist_sq is the max of all the
      // dist1_sq's else dist_sq is the min of all the dist1_sq's
      dist_sq = (fIter == fromCoords.begin())
                    ? dist1_sq
                    : (noe.from.type == NOE_AND) ? std::max(dist_sq, dist1_sq)
                                                 : std::min(dist_sq, dist1_sq);
    }
    return std::sqrt(dist_sq);
  }
}

double NmrSF::StdDistance(const StdRestraintAtoms &std) const {
  // Compile the list of coords at each end
  // If the restraint end is defined of type MEAN, just store the center of mass
  // of the atoms in the list If the restraint end is defined of type OR or AND,
  // store all coords
  CoordList fromCoords;
  if (std.from.type == NOE_MEAN)
    fromCoords.push_back(GetCenterOfAtomicMass(std.from.atoms));
  else
    fromCoords = GetCoordList(std.from.atoms);

  double dist_sq(999.9); // Keep track of minimum distance**2
  // Iterate over coords in each list and return the appropriate distance**2
  // between any of them
  for (CoordListConstIter fIter = fromCoords.begin(); fIter != fromCoords.end();
       fIter++) {
    // dist1_sq is the minimum distance**2 between the current "from" coord
    // and all the coords in the "to" list.
    double dist1_sq(999.9);
    // For STD restraints, the list of "to" coords comes from the indexing grid
    const AtomRList &toAtoms = m_spGrid->GetAtomList(*fIter);
    for (AtomRListConstIter tIter = toAtoms.begin(); tIter != toAtoms.end();
         tIter++) {
      double r12_sq = Length2(*fIter, (*tIter)->GetCoords());
      dist1_sq =
          (tIter == toAtoms.begin()) ? r12_sq : std::min(dist1_sq, r12_sq);
    }
    // dist_sq is the appropriate overall distance**2 for all the calculated
    // dist1_sq's i.e. if from.type==NOE_AND, dist_sq is the max of all the
    // dist1_sq's else dist_sq is the min of all the dist1_sq's
    dist_sq = (fIter == fromCoords.begin())
                  ? dist1_sq
                  : (std.from.type == NOE_AND) ? std::max(dist_sq, dist1_sq)
                                               : std::min(dist_sq, dist1_sq);
  }
  return std::sqrt(dist_sq);
}

void NmrSF::ParameterUpdated(const std::string &strName) {
  if (strName == _QUADRATIC) {
    m_bQuadratic = GetParameter(_QUADRATIC);
  } else {
    BaseSF::ParameterUpdated(strName);
  }
}
