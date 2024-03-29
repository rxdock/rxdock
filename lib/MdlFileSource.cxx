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

#include <functional>
#include <iomanip>

#include "rxdock/AtomFuncs.h"
#include "rxdock/FileError.h"
#include "rxdock/MdlFileSource.h"
#include "rxdock/ModelError.h"
#include "rxdock/Plane.h"

#include <fmt/ostream.h>
#include <loguru.hpp>

using namespace rxdock;

MdlFileSource::MdlFileSource(const std::string &fileName, bool bPosIonisable,
                             bool bNegIonisable, bool bImplHydrogens)
    : BaseMolecularFileSource(fileName, IDS_MDL_RECDELIM,
                              "MDL_FILE_SOURCE"), // Call base class constructor
      m_bPosIonisable(bPosIonisable), m_bNegIonisable(bNegIonisable),
      m_bImplHydrogens(bImplHydrogens) {
  // Open an Element data source
  m_spElementData = ElementFileSourcePtr(
      new ElementFileSource(GetDataFileName("data", "elements.json")));
  _RBTOBJECTCOUNTER_CONSTR_("MdlFileSource");
}

// Default destructor
MdlFileSource::~MdlFileSource() { _RBTOBJECTCOUNTER_DESTR_("MdlFileSource"); }

void MdlFileSource::Parse() {
  // Only parse if we haven't already done so
  if (!m_bParsedOK) {
    ClearMolCache(); // Clear current cache
    Read();          // Read the current record

    try {
      FileRecListIter fileIter = m_lineRecs.begin();
      FileRecListIter fileEnd = m_lineRecs.end();

      if ((*fileIter).size() && (*fileIter)[(*fileIter).size() - 1] == '\r') {
        throw FileParseError(_WHERE_, "Windows CRLF line endings detected in " +
                                          GetFileName());
      }

      //////////////////////////////////////////////////////////
      // 1a. Store title lines (first 3)..
      unsigned int nTitleRec = 3;
      m_titleList.reserve(nTitleRec); // Allocate enough memory for the vector
      while ((m_titleList.size() < nTitleRec) && (fileIter != fileEnd)) {
        m_titleList.push_back(*fileIter++);
      }

      // 1b ..and check we read them all before reaching the end of the file
      if (m_titleList.size() != nTitleRec)
        throw FileParseError(_WHERE_,
                             "Incomplete title records in " + GetFileName());

      // 1c DM 20 Dec 1999 - throw model error if record is 2D instead of 3D
      // DM 08 Aug 2000 - check for "3" at position 21 as well as "3D" at
      // position 20
      if ((m_titleList[1].find("3", 21) != 21) &&
          (m_titleList[1].find("3D", 20) != 20)) {
        throw ModelError(_WHERE_, "2D record in " + GetFileName());
      }

      //////////////////////////////////////////////////////////
      // 2. Read number of atoms and bonds
      unsigned int nAtomRec;
      unsigned int nBondRec;
      if (fileIter != fileEnd) {
        // The SD file format only uses a field width of 3 to store nAtoms,
        // nBonds so for values over 99 the two fields coalesce. Workaround is
        // to insert a space between the two fields (or use sscanf)
        if ((*fileIter).size() > 3)
          (*fileIter).insert(3, " ");
        std::istringstream istr(*fileIter++);
        istr >> nAtomRec >> nBondRec;
        LOG_F(1, "MdlFileSource::Parse: {} atoms, {} bonds", nAtomRec,
              nBondRec);
      } else
        throw FileParseError(_WHERE_, "Missing atom and bond information in " +
                                          GetFileName());

      //////////////////////////////////////////////////////////
      // 3a. Create and store atoms
      m_atomList.reserve(nAtomRec); // Allocate enough memory for the vector

      int nAtomId(0);
      Coord coord;                // X,Y,Z coords
      std::string strElementName; // element name
      int nMassDiff;              // mass difference
      int nFormalCharge; // formal charge (MDL stores in a funny way, see below)
      std::string strSegmentName(
          "H");                      // constant (i.e. one segment, one residue)
      std::string strSubunitId("1"); // constant
      std::string strSubunitName("MOL"); // constant

      while ((m_atomList.size() < nAtomRec) && (fileIter != fileEnd)) {
        std::istringstream istr(*fileIter++);
        istr >> coord.xyz(0) >> coord.xyz(1) >> coord.xyz(2) >>
            strElementName >> nMassDiff >> nFormalCharge;

        // Look up the element data
        ElementData elData = m_spElementData->GetElementData(strElementName);
        int nAtomicNo = elData.atomicNo;

        // Correct the formal charge (for non-zero charges, MDL stores as
        // 4-charge) i.e. -3 is stored as 7, -2 as 6, -1 as 5, 0 as 0 or 4 +1 as
        // 3, +2 as 2, +3 as 1
        if (nFormalCharge > 0)
          nFormalCharge = 4 - nFormalCharge;

        // Compose the atom name from element+atomID (i.e. C1, N2, C3 etc)
        nAtomId++;
        std::ostringstream ostr;
        ostr << strElementName << nAtomId;
        std::string strAtomName(ostr.str());

        // Construct a new atom (constructor only accepts the 2D params)
        AtomPtr spAtom(new Atom(nAtomId, nAtomicNo, strAtomName, strSubunitId,
                                strSubunitName, strSegmentName,
                                Atom::UNDEFINED, // Hybridisation state
                                0,               // Num implicit hydrogens
                                nFormalCharge)   // Formal charge
        );

        // Now set the remaining 2-D and 3-D params
        spAtom->SetCoords(coord);
        spAtom->SetPartialCharge(0.0);
        spAtom->SetGroupCharge(nFormalCharge); // Default is group charge=formal
                                               // charge (gets corrected later)
        spAtom->SetAtomicMass(elData.mass + nMassDiff);

        m_atomList.push_back(spAtom);
        m_segmentMap[strSegmentName]++; // increment atom count in segment map
      }

      // 3b ..and check we read them all before reaching the end of the file
      if (m_atomList.size() != nAtomRec)
        throw FileParseError(_WHERE_,
                             "Incomplete atom records in " + GetFileName());

      //////////////////////////////////////////////////////////
      // 4a. Create and store bonds
      m_bondList.reserve(nBondRec); // Allocate enough memory for the vector

      int nBondId(0);
      unsigned int idxAtom1;
      unsigned int idxAtom2;
      int nBondOrder;

      while ((m_bondList.size() < nBondRec) && (fileIter != fileEnd)) {
        // The SD file format only uses a field width of 3 to store atom1,atom2
        // so for values over 99 the two fields coalesce.
        // Workaround is to insert a space between the two fields (or use
        // sscanf)
        if ((*fileIter).size() > 3)
          (*fileIter).insert(3, " ");
        std::istringstream istr(*fileIter++);
        istr >> idxAtom1 >> idxAtom2 >> nBondOrder;
        if ((idxAtom1 > nAtomRec) ||
            (idxAtom2 > nAtomRec)) { // Check for indices in range
          throw FileParseError(_WHERE_,
                               "Atom index out of range in bond records in " +
                                   GetFileName());
        }
        AtomPtr spAtom1(
            m_atomList[idxAtom1 -
                       1]); // Decrement the atom index as the atoms are
                            // numbered from zero in our atom vector
        AtomPtr spAtom2(
            m_atomList[idxAtom2 -
                       1]); // Decrement the atom index as the atoms are
                            // numbered from zero in our atom vector
        BondPtr spBond(
            new Bond(++nBondId, // Store a nominal bond ID starting from 1
                     spAtom1, spAtom2, nBondOrder));
        m_bondList.push_back(spBond);
      }

      // 4b ..and check we read them all before reaching the end of the file
      if (m_bondList.size() != nBondRec)
        throw FileParseError(_WHERE_,
                             "Incomplete bond records in " + GetFileName());

      // DM 12 May 1999 - read data records (if any)
      for (; fileIter != fileEnd; fileIter++) {
        if ((*fileIter).find(">") == 0) { // Found a data record
          std::string::size_type ob =
              (*fileIter).find("<"); // First open bracket
          std::string::size_type cb =
              (*fileIter).rfind(">"); // Last closed bracket
          if ((ob != std::string::npos) && (cb != std::string::npos)) {
            std::string fieldName =
                (*fileIter).substr(ob + 1, cb - ob - 1); // Data field name
            std::vector<std::string> sl; // String list for storing data value
            while ((++fileIter != fileEnd) && !(*fileIter).empty()) {
              sl.push_back(*fileIter);
            }
            m_dataMap[fieldName] = Variant(sl);
          }
        }
      }

      // DM 18 May 1999 - check for <Name> record
      // If not present, define <Name> as the first title record
      if ((m_dataMap.find("Name") == m_dataMap.end()) && !m_titleList.empty()) {
        m_dataMap["Name"] = Variant(m_titleList.front());
      }

      // Setup the atomic params not stored in the file (e.g. hybridisation
      // state etc)
      SetupAtomParams();

      // DM 4 June 1999 - set up distinct segment names for each molecular
      // fragment in the record
      SetupSegmentNames();

      //////////////////////////////////////////////////////////
      // If we get this far everything is OK
      m_bParsedOK = true;
    }

    catch (Error &error) {
      ClearMolCache(); // Clear the molecular cache so we don't return
                       // incomplete atom and bond lists
      throw;           // Rethrow the Error
    }
  }
}

// Sets up all the atomic attributes that are not explicitly stored in the PSF
// file
void MdlFileSource::SetupAtomParams() {
  // In the SD/MOL file we have:
  //  Element names
  //  X,Y,Z coords
  //  Formal bond orders
  //  Formal charges
  //
  // We need to determine:
  //  Hybridisation state - sp3,sp2,sp,trigonal planar
  //  Valency of each atom - is it in range? - use to assign implicit hydrogens
  //  if necessary Corrected vdW radii on H-bonding hydrogens and extended atoms
  //  Simple "forcefield" type - e.g. Osp3-
  //  Interaction group charges

  // DM 2 Aug 1999 - remove all non-polar hydrogens if required
  if (m_bImplHydrogens) {
    RemoveNonPolarHydrogens();
    RenumberAtomsAndBonds(); // Ensure atom and bond IDs are consecutive prior
                             // to autoprotonation
  }

  SetupHybridState();
  SetupTheRest();
  SetCyclicFlags(); // DM 8 Feb 2000 - needed for SetupIonicGroups()
  if (m_bPosIonisable)
    SetupPosIonisableGroups();
  if (m_bNegIonisable)
    SetupNegIonisableGroups();
  SetupIonicGroups();
  RemoveZwitterions(m_atomList); // DM 6 Apr 1999 - might as well
                                 // neutralise any zwitterions here
  RenumberAtomsAndBonds();       // DM 30 Oct 2000 - clean up atom and bond list
                                 // numbering
}

// Helper functions for SetupAtomParams()

// Sets hybridisation state and checks for Nsp3,Osp3 or Ssp3 next to sp2 special
// cases
void MdlFileSource::SetupHybridState() {
  //////////////////
  // FIRST PASS
  //////////////////
  for (AtomListIter iter = m_atomList.begin(); iter != m_atomList.end();
       iter++) {
    // Determine hybridisation state from the maximum bond order to this atom
    // i.e. if it makes a triple bond it is SP, double bond=SP2, single
    // bond=SP3, no bonds=UNDEFINED
    switch ((*iter)->GetMaxFormalBondOrder()) {
    case 3:
      (*iter)->SetHybridState(Atom::SP);
      break;
    case 2:
      (*iter)->SetHybridState(Atom::SP2);
      break;
    case 1:
      (*iter)->SetHybridState(Atom::SP3);
      break;
    case 0:
      (*iter)->SetHybridState(Atom::UNDEFINED);
      break;
    default:
      // throw ModelError(_WHERE_, (*iter)->GetName()+" makes a bond with
      // an unknown bond order");
      LOG_F(ERROR, "{} makes a bond with an unknown bond order",
            (*iter)->GetName());
      (*iter)->SetHybridState(Atom::UNDEFINED);
      break;
    }
  }
  //////////////////
  // SECOND PASS - check for Nsp3, Osp3 and Ssp3 next to a pi atom. Convert to
  // trigonal planar For Nsp3, also check the actual planarity of the atom,
  // based on the distance from the nitrogen to the plane of the substituent. If
  // non-planar, then leave as Nsp3. Needed for amides, pyrroles, furans and
  // thiophenes. Trigonal planar atoms are sp2 hybridised with two electrons in
  // a pi-orbital, and only make formal single bonds. (as opposed to SP2 atoms
  // which make a formal double bond with only one electron in the pi-orbital)
  //////////////////
  //
  // Get a list of all the sp3 atoms which are not positively charged
  AtomList sp3AtomList =
      GetAtomListWithPredicate(m_atomList, isHybridState_eq(Atom::SP3));
  sp3AtomList =
      GetAtomListWithPredicate(sp3AtomList, std::not1(isAtomPosCharged()));

  for (AtomListIter iter = sp3AtomList.begin(); iter != sp3AtomList.end();
       iter++) {
    AtomList bondedAtomList;
    switch ((*iter)->GetAtomicNo()) {
    case 7:  // N
    case 8:  // O
    case 16: // S
      bondedAtomList = GetBondedAtomList(*iter);
      // Is is bonded to a pi-atom ?
      if (FindAtomInList(bondedAtomList, isPiAtom()) != bondedAtomList.end()) {
        // Determine distance from atom to plane of bonded atoms (only nitrogen
        // will have 3 substituents)
        if (bondedAtomList.size() == 3) {
          const Coord &c0 = (*iter)->GetCoords();
          Vector v1 = bondedAtomList[0]->GetCoords() - c0; // Vector to subst 1
          Vector v2 = bondedAtomList[1]->GetCoords() - c0; // Vector to subst 2
          Vector v3 = bondedAtomList[2]->GetCoords() - c0; // Vector to subst 3
          // Take unit vectors along each bond vector to eliminate the effect of
          // different bond lengths
          Plane p = Plane(c0 + v1.Unit(), c0 + v2.Unit(),
                          c0 + v3.Unit()); // Plane of three substs
          // This is the signed distance from the atom to the plane
          double dist = std::fabs(DistanceFromPointToPlane(c0, p));
          LOG_F(1, "Distance from {} to plane of substituents={}",
                (*iter)->GetName(), dist);
          if (dist > 0.1) { // Not planar
            LOG_F(1, "Not planar");
            break;
          }
        }
        (*iter)->SetHybridState(Atom::TRI);
        LOG_F(1, "Changing {} from SP3 to TRI", (*iter)->GetName());
      }
      break;
    default:
      break;
    }
  }
}

// Checks for valency within range
// Adds implicit hydrogens to C,N,S
// Defines vdW radius, correcting for extended atoms and H-bond donor hydrogens
// Defines formal "force-field" type string
void MdlFileSource::SetupTheRest() {
  // Radius increment for atoms with implicit hydrogens
  // DM 22 Jul 1999 - only increase the radius for sp3 atoms with implicit
  // hydrogens For sp2 and aromatic, leave as is
  isHybridState_eq bIsSP3(Atom::SP3);
  double dImplRadIncr = m_spElementData->GetImplicitRadiusIncr();
  // Radius increment and predicate for H-bonding hydrogens
  double dHBondRadIncr = m_spElementData->GetHBondRadiusIncr();
  isAtomHBondDonor bIsHBondDonor;
  // Element data for hydrogen
  ElementData elHData = m_spElementData->GetElementData(1);

  for (AtomListIter iter = m_atomList.begin(); iter != m_atomList.end();
       iter++) {
    // Get the element data for this atom
    int nAtomicNo = (*iter)->GetAtomicNo();
    ElementData elData = m_spElementData->GetElementData(nAtomicNo);
    double vdwRadius = elData.vdwRadius;
    // Determine valency
    int nFormalCharge = (*iter)->GetFormalCharge();
    int nTotalFormalBondOrder = (*iter)->GetTotalFormalBondOrder();
    int nValency = nTotalFormalBondOrder - nFormalCharge;
    // Is it within range ?
    int nMinVal = elData.minVal;
    int nMaxVal = elData.maxVal;
    // If it is too low, then add implicit hydrogens but only for C,N,or S
    if (nValency < nMinVal) {
      switch (nAtomicNo) {
      case 6:
        // case 7:
      case 16:
        (*iter)->SetNumImplicitHydrogens(nMinVal - nValency);
        (*iter)->SetAtomicMass(elData.mass +
                               (nMinVal - nValency) *
                                   elHData.mass); // Adjust atomic mass
        if (bIsSP3(*iter)) {
          vdwRadius +=
              dImplRadIncr; // adjust vdw radius (for sp3 implicit atoms only)
          LOG_F(1, "Adding {} implicit hydrogens to {}", nMinVal - nValency,
                (*iter)->GetName());
        }
        break;
      default:
        // DM 26 Jul 1999 - don't complain about valency errors
        // throw ModelError(_WHERE_, (*iter)->GetName()+" makes too few
        // bonds");
        LOG_F(ERROR, "{} makes too few bonds", (*iter)->GetName());
        break;
      }
    }
    // DM 26 Jul 1999 - don't complain about valency errors
    else if (nValency > nMaxVal) {
      // throw ModelError(_WHERE_, (*iter)->GetName()+" makes too many
      // bonds");
      LOG_F(ERROR, "{} makes too many bonds", (*iter)->GetName());
    }

    // Set the vdw radius, first correcting for H-bonding hydrogens
    if (bIsHBondDonor(*iter))
      vdwRadius += dHBondRadIncr;
    (*iter)->SetVdwRadius(vdwRadius);

    // Now set a nominal "force-field" type of the style OSP3- etc
    std::string strFFType = elData.element;
    strFFType += "_";
    strFFType += ConvertHybridStateToString((*iter)->GetHybridState());
    strFFType += ConvertFormalChargeToString((*iter)->GetFormalCharge());
    (*iter)->SetFFType(strFFType);
  }
}

// Protonates all amines, guanidines and imidazoles
void MdlFileSource::SetupPosIonisableGroups() {
  // Useful predicates
  isFFType_eq bIsN_SP3("N_SP3");
  isFFType_eq bIsN_SP2("N_SP2");
  isFFType_eq bIsN_TRI("N_TRI");
  isFFType_eq bIsC_SP2("C_SP2");

  // 1. Amines
  // Protonate all Nsp3 atoms
  // DM 14 May 1999 - unless adjacent to a pi-atom
  // Most will be N_TRI so won't be protonated anyway, but some will be
  // non-planar so are still N_SP3. However we don't want to protonate these
  // guys as the pKa for e.g. Ph-NH3+ is much lower than for R(aliphatic)-NH3+.
  // This will also prevent e.g. SO2NHR getting protonated.
  AtomList nsp3AtomList = GetAtomListWithPredicate(m_atomList, bIsN_SP3);
  // DM 04 Aug 1999 - only consider Nsp3 atoms which are bonded to at least one
  // Csp3 to avoid protonating hydrazines for e.g. This is the definition of an
  // amine in any case
  nsp3AtomList = GetAtomListWithPredicate(
      nsp3AtomList, std::not1(isCoordinationNumber_eq(0, "C_SP3")));

  for (AtomListIter iter = nsp3AtomList.begin(); iter != nsp3AtomList.end();
       iter++) {
    AtomList bondedAtomList = GetBondedAtomList(*iter);
    // Only protonate if not bonded to pi-atom ?
    if (FindAtomInList(bondedAtomList, isPiAtom()) == bondedAtomList.end()) {
      LOG_F(1, "Protonating neutral amine {}", (*iter)->GetName());
      AddHydrogen(*iter);
    } else {
      LOG_F(1, "Amine adjacent to pi-atom, not protonated: {}",
            (*iter)->GetName());
    }
  }

  // 2. Guanidine, imidazole, amidine
  // Guanidine
  // N_TRI    N_TRI
  //     `    '
  //      C_SP2
  //       ||
  //       N_SP2
  // - protonate the N_SP2
  //
  // Imidazole
  //    C_SP2-N_SP2
  //    //     `
  //   C_SP2   C_SP2
  //     `     '
  //      N_TRI
  //        |
  //        H
  // - protonate the N_SP2
  //

  // Loop over all C_SP2
  AtomList csp2AtomList = GetAtomListWithPredicate(m_atomList, bIsC_SP2);
  for (AtomListIter iter = csp2AtomList.begin(); iter != csp2AtomList.end();
       iter++) {
    AtomList bondedAtomList = GetBondedAtomList(*iter);

    AtomListIter nsp2Iter;
    AtomListIter ntriIter;
    AtomList nsp2BondedAtomList;
    AtomList ntriBondedAtomList;

    //////////////////
    // switch on #N_SP2
    //////////////////
    switch (GetNumAtomsWithPredicate(bondedAtomList, bIsN_SP2)) {

    case 1: // 1 x N_SP2
      nsp2Iter = FindAtomInList(bondedAtomList, bIsN_SP2);

      //////////////////
      // switch on #N_TRI
      //////////////////
      switch (GetNumAtomsWithPredicate(bondedAtomList, bIsN_TRI)) {

      //////////////////////////////////
      // 1 x N_SP2, 2 x N_TRI = guanidine
      //////////////////////////////////
      case 2:
        if (!(*iter)->GetCyclicFlag()) {
          LOG_F(1, "Protonating neutral guanidine {}", (*nsp2Iter)->GetName());
          AddHydrogen(*nsp2Iter);
        } else {
          LOG_F(1, "Cyclic guanidine-like fragment - ignored {}",
                (*nsp2Iter)->GetName());
        }
        break;

      //////////////////////////////////
      // 1 x N_SP2, 1 x N_TRI = possible imidazole, also need to check for 5
      // membered ring
      //////////////////////////////////
      case 1:
        ntriIter = FindAtomInList(bondedAtomList, bIsN_TRI);
        nsp2BondedAtomList = GetBondedAtomList(*nsp2Iter);
        ntriBondedAtomList = GetBondedAtomList(*ntriIter);
        LOG_F(1,
              "Possible imidazole/amidine, found 1xN_SP2, 1xN_TRI bonded to {}",
              (*iter)->GetName());
        // IMIDAZOLE - Check that the N_SP2 is bonded to 2 x C_SP2 and that
        // the N_TRI is bonded to 2 x C_SP2 and 1 hydrogen
        if ((GetNumAtomsWithPredicate(nsp2BondedAtomList, bIsC_SP2) == 2) &&
            (GetNumAtomsWithPredicate(ntriBondedAtomList, bIsC_SP2) == 2) &&
            (GetNumAtomsWithPredicate(ntriBondedAtomList, isAtomicNo_eq(1)) ==
             1)) {
          LOG_F(1, "Possible imidazole, bonding requirements for N_SP2 and "
                   "N_TRI met");
          // Now check if any of the atoms bonded to the N_TRI are 1-2 connected
          // to any of the atoms bonded to the N_SP2. If so, it is a 5-membered
          // ring DM 25 Jul 2002 - also check whether the atoms are bridgeheads
          // or not We don't want to protonate imidazoles which are part of
          // larger fused ring systems (can have very different pKa's)
          AtomList atoms12Conn;
          for (AtomListConstIter iter2 = ntriBondedAtomList.begin();
               iter2 != ntriBondedAtomList.end(); iter2++) {
            AtomList tmpList = GetAtomListWithPredicate(
                nsp2BondedAtomList, isAtom_12Connected(*iter2));
            std::copy(tmpList.begin(), tmpList.end(),
                      std::back_inserter(atoms12Conn));
          }
          if (atoms12Conn.size() == 1) {
            LOG_F(1, "5-membered imidazole ring found");
            // Remove bridgehead atoms (fused rings)
            atoms12Conn = GetAtomListWithPredicate(
                atoms12Conn, std::not1(isAtomBridgehead()));
            if (atoms12Conn.size() == 1) { // Imidazole!!
              LOG_F(1, "Protonating neutral imidazole {}",
                    (*nsp2Iter)->GetName());
              AddHydrogen(*nsp2Iter);
            } else {
              LOG_F(1, "Fused imidazole - ignored {}", (*nsp2Iter)->GetName());
            }
          } else {
            LOG_F(1, "Not a 5-membered ring");
          }
        }
        // AMIDINE Check that the N_SP2 is bonded to 1 hydrogen and the N_TRI is
        // bonded to 2 hydrogens
        else if ((GetNumAtomsWithPredicate(nsp2BondedAtomList,
                                           isAtomicNo_eq(1)) == 1) &&
                 (GetNumAtomsWithPredicate(ntriBondedAtomList,
                                           isAtomicNo_eq(1)) == 2)) {
          LOG_F(1, "Protonating neutral amidine {}", (*nsp2Iter)->GetName());
          AddHydrogen(*nsp2Iter);
        } else {
          LOG_F(1, "Bonding requirements for N_SP2 and N_TRI not met");
        }
        break;
        //////////////////////
        // End of case #N_TRI=1
        //////////////////////

      default:
        break;
      }
      //////////////////////
      // End of switch #N_TRI
      //////////////////////
      break;
      //////////////////////
      // End of case #N_SP2=1
      //////////////////////

    default:
      break;
    }
    //////////////////////
    // End of switch #N_SP2
    //////////////////////
  }
  // End of loop over all C_SP2
}

// Deprotonate all carboxylic acids, phosphates, phosphonates, sulphates,
// sulphonates
void MdlFileSource::SetupNegIonisableGroups() {
  // Useful predicates
  isFFType_eq bIsO_TRI("O_TRI");
  isFFType_eq bIsO_TRIM("O_TRI-");
  isFFType_eq bIsO_SP2("O_SP2");
  isFFType_eq bIsC_SP2("C_SP2");
  isFFType_eq bIsS_SP2("S_SP2");
  isFFType_eq bIsP_SP2("P_SP2");

  // Compile list of all C_SP2, S_SP2 and P_SP2
  AtomList csp2AtomList = GetAtomListWithPredicate(m_atomList, bIsC_SP2);
  AtomList ssp2AtomList = GetAtomListWithPredicate(m_atomList, bIsS_SP2);
  AtomList psp2AtomList = GetAtomListWithPredicate(m_atomList, bIsP_SP2);
  AtomList acidAtomList = csp2AtomList;
  acidAtomList.insert(acidAtomList.end(), ssp2AtomList.begin(),
                      ssp2AtomList.end());
  acidAtomList.insert(acidAtomList.end(), psp2AtomList.begin(),
                      psp2AtomList.end());

  // Loop over all potential acid centers
  for (AtomListIter iter = acidAtomList.begin(); iter != acidAtomList.end();
       iter++) {
    AtomList bondedAtomList = GetBondedAtomList(*iter);
    // Check if there is at least one O_SP2 bonded (i.e.C=O, S=O, P=O)
    // DM 27 Apr 1999 - check there is already an O_TRI- present. We only want
    // to deprotonate one oxygen.
    if ((GetNumAtomsWithPredicate(bondedAtomList, bIsO_SP2) > 0) &&
        (GetNumAtomsWithPredicate(bondedAtomList, bIsO_TRIM) == 0)) {
      // Now get the list of all O_TRI which are bonded (i.e. C-O-, S-O-, P-O-)
      AtomList otriAtomList =
          GetAtomListWithPredicate(bondedAtomList, bIsO_TRI);
      for (AtomListIter otriIter = otriAtomList.begin();
           otriIter != otriAtomList.end(); otriIter++) {
        // Is there a hydrogen bonded to the O_TRI ?
        AtomList otriBondedAtomList = GetBondedAtomList(*otriIter);
        AtomListIter hIter =
            FindAtomInList(otriBondedAtomList, isAtomicNo_eq(1));
        if (hIter != otriBondedAtomList.end()) {
          LOG_F(1, "Removing {} from {}, acid center={}", (*hIter)->GetName(),
                (*otriIter)->GetName(), (*iter)->GetName());
          RemoveAtom(*hIter);
          // Adjust the attributes on the O_TRI
          (*otriIter)->SetFormalCharge(-1);
          (*otriIter)->SetGroupCharge(-1.0);
          (*otriIter)->SetFFType("O_TRI-");
          break; // DM 19 Apr 1999 - only deprotonate the first O-H on any acid
                 // group
        }
      }
    }
  }
}

// Adds a proton to N_SP3 and N_SP2 atoms
void MdlFileSource::AddHydrogen(AtomPtr spAtom) {
  isFFType_eq bIsNSP3("N_SP3");
  isFFType_eq bIsNSP2("N_SP2");

  // Element data for hydrogen
  ElementData elHData = m_spElementData->GetElementData(1);
  // Radius increment for H-bonding hydrogens
  double dHBondRadIncr = m_spElementData->GetHBondRadiusIncr();

  // Attributes for new atom (to be defined)
  Coord hCoord;   // New coord
  double hRadius; // Vdw radius

  // Attributes for heavy atom
  const Coord &c0 = spAtom->GetCoords();
  AtomList bondedAtomList = GetBondedAtomList(spAtom);

  ////////////////
  // Nsp3
  // We determine the plane formed by the three substituents and add the
  // hydrogen along a vector normal to the plane (and in the other direction)
  //    H
  //    |
  //    N
  //   '|`
  ////////////////
  if (bIsNSP3(spAtom)) {
    if (bondedAtomList.size() != 3)
      throw ModelError(_WHERE_,
                       spAtom->GetName() + " has inconsistent atom type");
    double dNHLength = 1.00;                         // N-H bond length
    Vector v1 = bondedAtomList[0]->GetCoords() - c0; // Vector to subst 1
    Vector v2 = bondedAtomList[1]->GetCoords() - c0; // Vector to subst 2
    Vector v3 = bondedAtomList[2]->GetCoords() - c0; // Vector to subst 3
    // Take unit vectors along each bond vector to eliminate the effect of
    // different bond lengths
    Plane p = Plane(c0 + v1.Unit(), c0 + v2.Unit(),
                    c0 + v3.Unit()); // Plane of three substs
    // This is the signed distance from the nitrogen to the plane
    double dist = DistanceFromPointToPlane(c0, p);
    // Take sign so we go in the opposite direction to the other substituents
    int iSign = (dist > 0) ? 1 : -1;
    hCoord = c0 + (iSign * dNHLength * p.VNorm()); // New H coordinate
    hRadius = elHData.vdwRadius + dHBondRadIncr;   // New H vdw radius
    // Adjust the attributes on the Nitrogen
    spAtom->SetFormalCharge(1);
    spAtom->SetGroupCharge(1.0);
    spAtom->SetFFType("N_SP3+");
  }

  ////////////////
  // Nsp2
  // We determine the average of the bond vectors to the two substituents and
  // add the hydrogen along this vector in the other direction
  //   H
  //   |
  //   N
  //  ' `
  ////////////////
  else if (bIsNSP2(spAtom)) {
    if (bondedAtomList.size() != 2)
      throw ModelError(_WHERE_,
                       spAtom->GetName() + " has inconsistent atom type");
    double dNHLength = 1.00;                         // N-H bond length
    Vector v1 = bondedAtomList[0]->GetCoords() - c0; // Vector to subst 1
    Vector v2 = bondedAtomList[1]->GetCoords() - c0; // Vector to subst 2
    // Take unit vectors along each bond vector to eliminate the effect of
    // different bond lengths
    Vector vMean = (v1.Unit() + v2.Unit()) / 2.0;
    hCoord = c0 - (dNHLength * vMean.Unit());    // New H coordinate
    hRadius = elHData.vdwRadius + dHBondRadIncr; // New H vdw radius
    // Adjust the attributes on the Nitrogen
    spAtom->SetFormalCharge(1);
    spAtom->SetGroupCharge(1.0);
    spAtom->SetFFType("N_SP2+");
  }

  ////////////////
  // Don't consider other element types for now
  ////////////////
  else
    return;

  // Construct the new hydrogen atom (constructor only accepts the 2D params)
  int nAtomId = m_atomList.size() + 1;
  std::ostringstream ostr;
  ostr << "H" << nAtomId;
  std::string strAtomName(ostr.str());
  AtomPtr spHAtom(new Atom(nAtomId,
                           1, // nAtomicNo,
                           strAtomName,
                           "1",       // strSubunitId,
                           "MOL",     // strSubunitName,
                           "H",       // strSegmentName,
                           Atom::SP3, // Hybridisation state
                           0,         // Num implicit hydrogens
                           0)         // Formal charge
  );

  // Now set the remaining 2-D and 3-D params
  spHAtom->SetCoords(hCoord);
  spHAtom->SetPartialCharge(0.0);
  spHAtom->SetGroupCharge(0.0);
  spHAtom->SetAtomicMass(elHData.mass);
  spHAtom->SetVdwRadius(hRadius);
  m_atomList.push_back(spHAtom);

  // Construct a new bond to the hydrogen
  int nBondId = m_bondList.size() + 1;
  BondPtr spBond(new Bond(nBondId, spAtom, spHAtom, 1)); // Bond order = 1
  m_bondList.push_back(spBond);
}

// SetupIonicGroups sets the "group charge" attribute of the atoms. The idea is
// that the "group charge" is file-format independent, and eliminates the need
// for the scoring function setup to have any knowledge of the input file
// format. For example, PSF files represent COO- as OC-C-OC (both oxygens
// charged) whereas SD files represent it as the formal resonance structure
// O=C-O-. For rDock we actually want the negative charge on the central carbon!
//
// It is envisaged each molecular file source will have a version of
// SetupIonicGroups to convert from the native representation to the rDock
// representation. The FormalCharge and PartialCharge attributes should be left
// untouched to allow the model to be rewritten back in the same format.
void MdlFileSource::SetupIonicGroups() {
  // Group charges have been initialised to the formal charge by Parse() so we
  // only need to check for special cases here.
  SetupNSP3Plus();
  SetupNSP2Plus();
  SetupOTRIMinus();
}

// Helper function to set up protonated amines
void MdlFileSource::SetupNSP3Plus() {
  isFFType_eq bIsN_SP3p("N_SP3+");
  isAtomicNo_eq bIsH(1);

  // Divide up the charge equally amongst the hydrogens bonded to the nitrogen
  AtomList nsp3pAtomList = GetAtomListWithPredicate(m_atomList, bIsN_SP3p);
  for (AtomListIter iter = nsp3pAtomList.begin(); iter != nsp3pAtomList.end();
       iter++) {
    // Get the list of hydrogens
    AtomList hbdAtomList =
        GetAtomListWithPredicate(GetBondedAtomList(*iter), bIsH);
    double nCharge = (*iter)->GetGroupCharge();      // Total charge
    double pCharge = nCharge / (hbdAtomList.size()); // Partial charge
    (*iter)->SetGroupCharge(0.0);                    // Neutralise the nitrogen
    for (AtomListIter hIter = hbdAtomList.begin(); hIter != hbdAtomList.end();
         hIter++) {
      (*hIter)->SetGroupCharge(pCharge); // Charge up the hydrogens
      LOG_F(1, "Transferring charge of {} from {} to {}", pCharge,
            (*iter)->GetName(), (*hIter)->GetName());
    }
  }
}

// Helper function to set up guanidines, imidazoles, amidines
// DM 23 Nov 2000 - also handles protonated pyridines correctly (transfers +1
// charge from N to H)
void MdlFileSource::SetupNSP2Plus() {
  // Useful predicates
  isFFType_eq bIsN_SP2p("N_SP2+");
  isFFType_eq bIsN_TRI("N_TRI");
  isFFType_eq bIsC_SP2("C_SP2");
  isAtomicNo_eq bIsH(1);

  // Loop over all N_SP2+
  AtomList nsp2PlusAtomList = GetAtomListWithPredicate(m_atomList, bIsN_SP2p);
  for (AtomListIter iter = nsp2PlusAtomList.begin();
       iter != nsp2PlusAtomList.end(); iter++) {
    // Now get the list of all bonded C_SP2
    AtomList csp2AtomList =
        GetAtomListWithPredicate(GetBondedAtomList(*iter), bIsC_SP2);
    LOG_F(1, "{} CSP2 atoms found bonded to NSP2+", csp2AtomList.size());
    // Eliminate bridgehead carbons, so we can guarantee we handle imidazoles
    // correctly
    csp2AtomList =
        GetAtomListWithPredicate(csp2AtomList, std::not1(isAtomBridgehead()));
    LOG_F(1, "{} are non-bridgehead", csp2AtomList.size());
    // There may be cases where the guanidinium and imidazole bonding patterns
    // are both present (on different C_SP2 atoms) so we want to check for both.
    // Let the guanidinium moiety take precedence
    AtomListIter guanIter =
        FindAtomInList(csp2AtomList, isCoordinationNumber_eq(2, "N_TRI"));
    AtomListIter imidIter =
        FindAtomInList(csp2AtomList, isCoordinationNumber_eq(1, "N_TRI"));
    AtomListIter csp2Iter =
        csp2AtomList.end(); // This will point to the C_SP2 atom to transfer the
                            // charge to
    if (guanIter != csp2AtomList.end())
      csp2Iter = guanIter;
    else if (imidIter != csp2AtomList.end())
      csp2Iter = imidIter;

    LOG_F(1, "Checking N_SP2+ atom {}", (*iter)->GetName());
    if (guanIter != csp2AtomList.end()) {
      LOG_F(1, "{} is a guanidinium-like carbon", (*guanIter)->GetName());
    }
    if (imidIter != csp2AtomList.end()) {
      LOG_F(1, "{} is a imidazole-like carbon", (*imidIter)->GetName());
    }

    // DM 23 Nov 2000 - move this code out of the following if {} block, to
    // handle protonated pyridines Get the list of hydrogens on the attached
    // nitrogens Start with the hydrogen on the N_SP2+
    AtomList hbdAtomList =
        GetAtomListWithPredicate(GetBondedAtomList(*iter), bIsH);
    LOG_F(1, "Found {} hydrogens bonded to {}", hbdAtomList.size(),
          (*iter)->GetName());

    // DM 23 Nov 2000 - this bit handles guanidines, imidazoles, amidines
    if (csp2Iter != csp2AtomList.end()) {
      // Now add the hydrogens on the (1 or 2) N_TRI atoms
      AtomList ntriAtomList =
          GetAtomListWithPredicate(GetBondedAtomList(*csp2Iter), bIsN_TRI);
      LOG_F(1, "Found {} nitrogens bonded to {}", ntriAtomList.size(),
            (*csp2Iter)->GetName());
      for (AtomListConstIter nIter = ntriAtomList.begin();
           nIter != ntriAtomList.end(); nIter++) {
        AtomList hAtomList =
            GetAtomListWithPredicate(GetBondedAtomList(*nIter), bIsH);
        LOG_F(1, "Found {} hydrogens bonded to {}", hAtomList.size(),
              (*nIter)->GetName());
        hbdAtomList.insert(hbdAtomList.end(), hAtomList.begin(),
                           hAtomList.end());
      }
      hbdAtomList.push_back(
          *csp2Iter); // Add the central carbon to the list of charged atoms
    }

    // Divide up the charge
    int nH = hbdAtomList.size();
    // DM 23 Nov 2000 - only transfer the charge if there is at least one
    // attached H
    if (nH > 0) {
      double nCharge = (*iter)->GetGroupCharge();
      double pCharge = nCharge / nH;
      (*iter)->SetGroupCharge(0.0); // N_SP2+
      for (AtomListIter hIter = hbdAtomList.begin(); hIter != hbdAtomList.end();
           hIter++) {
        (*hIter)->SetGroupCharge(pCharge); // H
        LOG_F(1, "Transferring charge of {} from {} to {}", pCharge,
              (*iter)->GetName(), (*hIter)->GetName());
      }
    }
  }
}

// Helper function to set up deprotonated acids
void MdlFileSource::SetupOTRIMinus() {
  isFFType_eq bIsO_TRIm("O_TRI-");
  isFFType_eq bIsO_SP2("O_SP2");
  isFFType_eq bIsC_SP2("C_SP2");
  isFFType_eq bIsS_SP2("S_SP2");
  isFFType_eq bIsP_SP2("P_SP2");

  // Compile list of all C_SP2, S_SP2 and P_SP2
  AtomList csp2AtomList = GetAtomListWithPredicate(m_atomList, bIsC_SP2);
  AtomList ssp2AtomList = GetAtomListWithPredicate(m_atomList, bIsS_SP2);
  AtomList psp2AtomList = GetAtomListWithPredicate(m_atomList, bIsP_SP2);
  AtomList acidAtomList = csp2AtomList;
  acidAtomList.insert(acidAtomList.end(), ssp2AtomList.begin(),
                      ssp2AtomList.end());
  acidAtomList.insert(acidAtomList.end(), psp2AtomList.begin(),
                      psp2AtomList.end());

  // Look for acid centers which have at least one bonded O_SP2 and at least one
  // bonded O_TRI-
  // e.g. C=O
  //     |
  //     O-
  // Divide total charge equally over the oxygens

  // Loop over all potential acid centers
  for (AtomListConstIter iter = acidAtomList.begin();
       iter != acidAtomList.end(); iter++) {
    // Get the list of bonded O_SP2 and O_TRI-
    AtomList osp2AtomList =
        GetAtomListWithPredicate(GetBondedAtomList(*iter), bIsO_SP2);
    AtomList otrimAtomList =
        GetAtomListWithPredicate(GetBondedAtomList(*iter), bIsO_TRIm);
    if ((osp2AtomList.size() > 0) && (otrimAtomList.size() > 0)) {
      // Combine O_SP2 and O_TRI- lists into one
      AtomList oAtomList = osp2AtomList;
      oAtomList.insert(oAtomList.end(), otrimAtomList.begin(),
                       otrimAtomList.end());
      // Total charge is the number of O_TRI-
      double nCharge =
          static_cast<double>(otrimAtomList.size()); // Total charge
      // Partial charge is total / divided by number of oxygens (O_SP2 + O_TRI-)
      double pCharge = -nCharge / (oAtomList.size()); // Partial charge
      for (AtomListIter oIter = oAtomList.begin(); oIter != oAtomList.end();
           oIter++) {
        (*oIter)->SetGroupCharge(pCharge); // Charge up the oxygens
        LOG_F(1, "Transferring charge of {} from {} to {}", pCharge,
              (*iter)->GetName(), (*oIter)->GetName());
      }
    }
  }
}

// DM 4 June 1999 - search for the separate fragments in the record
// and set distinct segment names for each one
// Segment names are initially H1,H2,H3...
// except that the largest segment is renamed back to H
void MdlFileSource::SetupSegmentNames() {
  std::string strLargestSegName; // Name of the largest segment
  int nMaxSize(0);               // Size of the largest segment
  int nSeg;
  AtomListIter seed;
  for (nSeg = 1, seed = m_atomList.begin(); seed != m_atomList.end();
       nSeg++, seed = FindAtomInList(m_atomList, isSegmentName_eq("H"))) {
    // New segment name (H1, H2 etc)
    std::ostringstream ostr;
    ostr << "H" << nSeg;
    std::string strSegName(ostr.str());
    // Temporary atom list containing atoms to be processed
    // Note: this is a true list (not a vector) as we will be making numerous
    // insertions and deletions
    (*seed)->SetSegmentName(strSegName);
    AtomTrueList pendingAtomList;
    pendingAtomList.push_back(*seed);
    // While we still have atoms to process
    // DM 12 May 2000 - correct nSize for ANSI scoping
    int nSize(0);
    for (nSize = 0; !pendingAtomList.empty(); nSize++) {
      // Take the last atom from the list and remove it
      AtomPtr spAtom = pendingAtomList.back();
      pendingAtomList.pop_back();
      // Get the list of bonded atoms which haven't yet been processed (i.e.
      // those whose segName is still H)
      AtomList atomsToAdd = GetAtomListWithPredicate(GetBondedAtomList(spAtom),
                                                     isSegmentName_eq("H"));
      // Set the segment name for each atom, and add to the queue
      for (AtomListIter iter = atomsToAdd.begin(); iter != atomsToAdd.end();
           iter++) {
        (*iter)->SetSegmentName(strSegName);
        pendingAtomList.push_back(*iter);
      }
    }
    LOG_F(1, "Segment {} contains {} atoms", strSegName, nSize);

    // Keep track of the largest segment
    if (nSize > nMaxSize) {
      nMaxSize = nSize;
      strLargestSegName = strSegName;
    }
  }

  // Rename the largest segment back to H
  LOG_F(1, "Largest segment is {}", strLargestSegName);
  AtomList largestAtomList =
      GetAtomListWithPredicate(m_atomList, isSegmentName_eq(strLargestSegName));
  for (AtomListIter iter = largestAtomList.begin();
       iter != largestAtomList.end(); iter++) {
    (*iter)->SetSegmentName("H");
  }

  // Now update the segment map
  m_segmentMap.clear();
  for (AtomListIter iter = m_atomList.begin(); iter != m_atomList.end();
       iter++) {
    m_segmentMap[(*iter)->GetSegmentName()]++;
  }
}

// DM 2 Aug 1999 - remove all non-polar hydrogens
void MdlFileSource::RemoveNonPolarHydrogens() {
  // Get list of all hydrogens
  AtomList hList = GetAtomListWithPredicate(m_atomList, isAtomicNo_eq(1));
  // Get list of all hydrogens bonded to carbon
  AtomList implHList =
      GetAtomListWithPredicate(hList, isCoordinationNumber_eq(1, 6));

  // Remove them one by one (not very efficient)
  for (AtomListIter iter = implHList.begin(); iter != implHList.end(); iter++) {
    RemoveAtom(*iter);
  }
  for (BondListConstIter bIter = m_bondList.begin(); bIter != m_bondList.end();
       bIter++) {
    LOG_F(1, "Bond ID={} ({}-{})", (*bIter)->GetBondId(),
          (*bIter)->GetAtom1Ptr()->GetName(),
          (*bIter)->GetAtom2Ptr()->GetName());
  }
}

// DM 8 Feb 2000 - setup atom and bond cyclic flags (previously in Model)
// Set the atom and bond cyclic flags for all atoms and bonds in the source
void MdlFileSource::SetCyclicFlags() {
  SetAtomAndBondCyclicFlags(m_atomList, m_bondList);
}

void rxdock::to_json(json &j, const MdlFileSource &mdlFileSrc) {
  j = json{
      {"sp-element-data", *mdlFileSrc.m_spElementData},
      {"pos-is-ionisable", mdlFileSrc.m_bPosIonisable},
      {"neg-ionisable", mdlFileSrc.m_bNegIonisable},
      {"imp-hydrogens", mdlFileSrc.m_bImplHydrogens},
  };
}

void rxdock::from_json(const json &j, MdlFileSource &mdlFileSrc) {
  j.at("sp-element-data").get_to(*mdlFileSrc.m_spElementData);
  j.at("pos-is-ionisable").get_to(mdlFileSrc.m_bPosIonisable);
  j.at("neg-ionisable").get_to(mdlFileSrc.m_bNegIonisable);
  j.at("imp-hydrogens").get_to(mdlFileSrc.m_bImplHydrogens);
}