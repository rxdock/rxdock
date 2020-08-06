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

// XB cctype for check atom in MOL2 is number
//#include <cctype>
#include <functional>
#include <sstream>

#include "rxdock/AtomFuncs.h"
#include "rxdock/FileError.h"
#include "rxdock/MOL2FileSource.h"
#include "rxdock/ModelError.h"

#include <fmt/ostream.h>
#include <loguru.hpp>

using namespace rxdock;

std::string MOL2FileSource::_CT("MOL2FileSource");
// record delimiter strings
std::string MOL2FileSource::_TRIPOS_DELIM("@<TRIPOS>");
unsigned int
    MOL2FileSource::_TRIPOS_DELIM_SIZE(MOL2FileSource::_TRIPOS_DELIM.size());
std::string MOL2FileSource::_IDS_MOL2_RECDELIM("@<TRIPOS>MOLECULE");

// straightforward ctor
MOL2FileSource::MOL2FileSource(const std::string &fileName, bool bImplHydrogens)
    : BaseMolecularFileSource(fileName, _IDS_MOL2_RECDELIM, "MOL2_FILE_SOURCE"),
      m_bImplHydrogens(bImplHydrogens) {
  m_spElementData = ElementFileSourcePtr(
      new ElementFileSource(GetDataFileName("data", "Elements.dat")));
  m_spParamSource = ParameterFileSourcePtr(
      new ParameterFileSource(GetDataFileName("data/sf", "IonicAtoms.prm")));
  m_NL = 0; // make line counter zero
  // initialize MOLECULE descriptors
  nAtoms = nBonds = nSubstructures = nFeatures = nSets = 0;
  _RBTOBJECTCOUNTER_CONSTR_("MOL2FileSource");
}

MOL2FileSource::~MOL2FileSource() {
  _RBTOBJECTCOUNTER_DESTR_("MOL2FileSource");
}

// parsing file and setting up atom descriptors
// at the beginning nothing is set, at the end everything should be set
// correctly
void MOL2FileSource::Parse(void) {
  // start parsing if we haven't already done so
  if (!m_bParsedOK) {
    // clear current cache
    ClearMolCache();
    // read the record, setting the delimiter to the beginning
    Read(false);
    // create a switchboard map - each member is a pointer to a function
    // (a sort of functor, but I prefer to call using a meaningful string
    // instead a number)
    std::map<std::string, fcPtr> theSwitchBoard;
    theSwitchBoard["MOLECULE"] = &MOL2FileSource::ParseRecordMOLECULE;
    theSwitchBoard["ATOM"] = &MOL2FileSource::ParseRecordATOM;
    theSwitchBoard["BOND"] = &MOL2FileSource::ParseRecordBOND;
    theSwitchBoard["SUBSTRUCTURE"] = &MOL2FileSource::ParseRecordSUBSTRUCTURE;
    theSwitchBoard["UNSUPPORTED"] = &MOL2FileSource::ParseRecordUNSUPPORTED;
    // the real parser
    try {
      // it would be nice to have a regexp stuff here
      // but the patterns are rather simple so it is not that bad
      // --------------------------------------------------------
      // declare the callback switched parser
      void (MOL2FileSource::*theCurrentParser)(const std::string &);
      // first assume the parser function is that that handles the MOLECULE
      // record
      theCurrentParser = &MOL2FileSource::ParseRecordMOLECULE;
      // than get lines one-by-one
      for (FileRecListIter lineIter = m_lineRecs.begin();
           lineIter < m_lineRecs.end(); ++lineIter) {
        // Ignore blank lines and comment lines
        if ((*lineIter).empty() || ((*lineIter).at(0) == '#')) {
          LOG_F(1, "Skipping blank/comment record");
          continue;
        }
        // check wether is it a record field delimiter (tag like:
        // @<TRIPOS>BLAHBLAH )line
        std::string theNewTag = GetMOL2Tag((*lineIter));
        if (!theNewTag.empty()) {
          // check key existence
          std::map<std::string, fcPtr>::const_iterator i =
              theSwitchBoard.find(theNewTag);
          if (i == theSwitchBoard.end()) {
            LOG_F(1, "Skipping unsupported record: {}", theNewTag);
            theCurrentParser = theSwitchBoard["UNSUPPORTED"];
          } else {
            LOG_F(INFO, "Reading {}", theNewTag);
            theCurrentParser = theSwitchBoard[theNewTag];
          }
        } else { // call the parser
          (*this.*theCurrentParser)((*lineIter));
        }
      }

      // we are ready with parsing
      // let's assign chain ids to atoms
      // and sort atomlist by substructure ID
      m_atomList.clear();
      for (std::map<int, AtomList>::const_iterator iter = m_ssAtoms.begin();
           iter != m_ssAtoms.end(); iter++) {
        int subst_id = iter->first;
        AtomList ssAtomList = iter->second;
        LOG_F(1, "Found {} atoms in SUBSTRUCTURE #{}", ssAtomList.size(),
              subst_id);
        std::copy(ssAtomList.begin(), ssAtomList.end(),
                  std::back_inserter(m_atomList));
        MOL2SubstructureMapConstIter ssIter = m_ssInfo.find(subst_id);
        if (ssIter != m_ssInfo.end()) {
          MOL2Substructure ss = ssIter->second;
          std::string subst_name = ss.GetName();
          std::string chain = ss.GetChain();
          std::string sub_type = ss.GetType();
          std::string sID, sName;
          GetSSIDandName(subst_name, subst_id, sID, sName);
          for (AtomListIter aIter = ssAtomList.begin();
               aIter != ssAtomList.end(); aIter++) {
            std::string atom_sID = (*aIter)->GetSubunitId();
            std::string atom_sName = (*aIter)->GetSubunitName();
            if ((sID != atom_sID) || (sName != atom_sName)) {
              LOG_F(WARNING,
                    "MOL2FileSource::Parse: Mismatch between SUBSTRUCTURE ({} "
                    "{}) and ATOM ({} {}) info",
                    sName, sID, atom_sName, atom_sID);
            }
            if ((atom_sName == "UNK") && (sub_type != "UNK")) {
              LOG_F(INFO,
                    "MOL2FileSource::Parse: Could use sub_type ({}) to update "
                    "SubunitName for {}",
                    sub_type, (*aIter)->GetFullAtomName());
            }
            (*aIter)->SetSegmentName(chain);
            LOG_F(1, "Updating segment name to {} for {}", chain,
                  (*aIter)->GetFullAtomName());
          }
        } else {
          LOG_F(1,
                "MOL2FileSource::Parse: No SUBSTRUCTURE record found for "
                "subst_id={}",
                subst_id);
        }
        // Set the group charges on a residue-by-residue basis
        // so we can check for presence or absence of all required atoms
        // e.g. in HIS need both polar Hs (HD1/HE2) to be present
        SetupPartialIonicGroups(ssAtomList, m_spParamSource);
      }

      if (m_atomList.size() != nAtoms) {
        LOG_F(WARNING,
              "MOL2FileSource::Parse: Mismatch in #ATOMS read: Expected = {}; "
              "Actual = {}",
              nAtoms, m_atomList.size());
      }
      if (m_bondList.size() != nBonds) {
        LOG_F(WARNING,
              "MOL2FileSource::Parse: Mismatch in #BONDS read: Expected = {}; "
              "Actual = {}",
              nBonds, m_bondList.size());
      }
      if (m_ssInfo.size() != nSubstructures) {
        LOG_F(WARNING,
              "MOL2FileSource::Parse: Mismatch in #SUBSTRUCTURES read: "
              "Expected = {}; Actual = {}",
              nSubstructures, m_ssInfo.size());
      }

      SetupAtomParams();

      for (AtomListConstIter iter = m_atomList.begin();
           iter != m_atomList.end(); ++iter) {
        m_segmentMap[(*iter)->GetSegmentName()]++; // increment atom count in
                                                   // segment map
        TriposAtomType::eType tt_file = (*iter)->GetTriposType();
        TriposAtomType::eType tt_auto = m_typer(*iter, true);
        if (tt_file == TriposAtomType::UNDEFINED) {
          LOG_F(WARNING, "MOL2FileSource::Parse: Undefined Tripos type for {}",
                (*iter)->GetFullAtomName());
        }
        if (tt_file != tt_auto) {
          LOG_F(INFO,
                "MOL2FileSource::Parse: Disagreement in Tripos atom types "
                "for {} File={}; Auto={}",
                (*iter)->GetFullAtomName(), m_typer.Type2Str(tt_file),
                m_typer.Type2Str(tt_auto));
        }
      }

      m_bParsedOK = true;
    } catch (Error &error) {
      ClearMolCache();
      throw;
    }
  }
}

// @<TRIPOS>MOLECULE
// six data lines (after tag)
// ------------------------------
// 1 name
// ------------------------------
// 2 numer of
// atoms [ bonds substructures features sets ]
// ------------------------------
// 3 type
// [SMALL|PROTEIN|NUCLEIC_ACID|SACCHARIDE]
// ------------------------------
// 4 type of charges
// [NO_CHARGES|DEL_RE|GASTEIGER|GAST_HUCK|HUCKEL|PULLMAN|GAUSS80_CHARGES
// |AMPAC_CHARGES|MULLIKEN_CHARGES|DICT_CHARGES|MMFF94_CHARGES|USER_CHARGES]
// ------------------------------
// 5 status_bits internal SYBYL stuff, should never be set by the user
// ------------------------------
// 6 one comment line
// ------------------------------
//
// it is the m_titleList that is set up here only
//
void MOL2FileSource::ParseRecordMOLECULE(const std::string &aLine) {
  LOG_F(1, "MOLECULE {}\t{}", aLine, m_NL);
  switch (m_NL) {
  case 0:                         // first line
    m_titleList.push_back(aLine); // put it into title
    m_NL++;                       // increment record line counter
    break;
  case 1: // second line
    ParseCountFields(aLine);
    m_NL++;
    break;
  default: // do nothing
    break;
  }
}

void MOL2FileSource::ParseRecordATOM(const std::string &aLine) {
  std::vector<std::string> tokens;
  Tokenize(aLine, tokens);
  // there must be at least 6 fields, rest is optional
  if (tokens.size() < 6)
    throw FileParseError(_WHERE_,
                         "Corrupted MOL2 file: not enough fields in ATOM ");

  // Compulsory fields
  int atom_id = std::atoi(tokens[0].c_str());
  std::string atom_name = tokens[1];
  Coord atom_coord;
  atom_coord.xyz(0) = std::atof(tokens[2].c_str());
  atom_coord.xyz(1) = std::atof(tokens[3].c_str());
  atom_coord.xyz(2) = std::atof(tokens[4].c_str());
  std::string atom_type = tokens[5];

  // Optional fields.
  int subst_id = (tokens.size() > 6) ? std::atoi(tokens[6].c_str()) : 1;
  std::string subst_name = (tokens.size() > 7) ? tokens[7] : "****";
  std::string sID, sName;
  GetSSIDandName(subst_name, subst_id, sID, sName);
  double charge = (tokens.size() > 8) ? std::atof(tokens[8].c_str()) : 0.0;
  // XB reweighting parameters
  //	Double wxb = (tokens.size() > 9) ? std::atof(tokens[9].c_str())
  //: 1.0;
  // XB mod for only doing if number, not characters
  //	Double wxb;
  //  if (tokens.size() > 9 && isdigit(tokens[9][0])){
  //	  wxb = std::atof(tokens[9].c_str());
  //	} else {
  //	  wxb = 1.0;
  //	}
  // end XB reweighting parameters

  // Derived atom params (some may be updated later)
  TriposAtomType::eType tt = m_typer.Str2Type(atom_type);
  int atomic_number = m_typer.Type2AtomicNo(tt);
  // Workaround for metals such as Zn, Mn, Mg
  // We should add Tripos types for these metals, but the most important thing
  // is to pick up the atomic number (can do this from the atom type string)
  if (tt == TriposAtomType::UNDEFINED) {
    std::string el = atom_type.substr(0, 2);
    LOG_F(INFO,
          "MOL2FileSource::ParseRecordATOM: Attempting to identify element by "
          "1st 2 chars of FFType string: {}",
          el);
    if (m_spElementData->isElementNamePresent(el)) {
      atomic_number = m_spElementData->GetElementData(el).atomicNo;
      LOG_F(INFO, "MOL2FileSource::ParseRecordATOM: Atomic number = {}",
            atomic_number);
    } else {
      LOG_F(INFO, "MOL2FileSource::ParseRecordATOM: Element not found");
    }
  }
  ElementData elementData = m_spElementData->GetElementData(atomic_number);
  Atom::eHybridState hybrid_state = (Atom::eHybridState)m_typer.Type2Hybrid(tt);
  AtomPtr newAtom(new Atom(atom_id, atomic_number, atom_name,
                           sID,     // subunit id undef
                           sName,   // subunit name undef
                           "UNK")); // segment name undefined at this
                                    // point (need SUBSTRUCTURE records)
  newAtom->SetCoords(atom_coord);
  newAtom->SetFFType(atom_type); // Original Tripos type in string format
  newAtom->SetTriposType(
      tt); // Same type, but stored as an enum (may be updated later)
  newAtom->SetPartialCharge(charge);
  newAtom->SetHybridState(hybrid_state);
  newAtom->SetAtomicMass(elementData.mass);
  //  newAtom->SetReweight(wxb); //XB

  m_atomList.push_back(newAtom);
  m_ssAtoms[subst_id].push_back(newAtom);
  LOG_F(1, "MOL2FileSource::ParseRecordATOM: {}", *newAtom);
}

void MOL2FileSource::ParseRecordBOND(const std::string &aLine) {
  std::vector<std::string> tokens;
  Tokenize(aLine, tokens);
  // there must be at least 4 fields, rest is optional
  if (tokens.size() < 4)
    throw FileParseError(_WHERE_,
                         "Corrupted MOL2 file: not enough fields in BOND");

  // Compulsory fields
  unsigned int bond_id = std::atoi(tokens[0].c_str());
  unsigned int origin_bond_id = std::atoi(tokens[1].c_str());
  unsigned int target_bond_id = std::atoi(tokens[2].c_str());
  std::string bond_type = tokens[3];

  // In rDock we only have integer bond order
  // Aromatics should be defined as alternate double/single bonds as in MDL mol
  // files Can't do much here so define everything as single, double or triple
  // We have the hybridisation state of each atoms already (from the Tripos
  // type) so are not so dependent on bond orders to deduce hybridisation
  int bond_order = (bond_type == "3") ? 3 : (bond_type == "2") ? 2 : 1;
  // Check that atom IDs are in range
  // This assumes that the Bond records are read after the Atom records (should
  // be always the case)
  if ((origin_bond_id < 1) || (origin_bond_id > m_atomList.size()) ||
      (target_bond_id < 1) || (target_bond_id > m_atomList.size()))
    throw FileParseError(_WHERE_,
                         "atom_id out of range in BOND record\n" + aLine);

  // Decrement the atom indices as the atoms are numbered from zero in our atom
  // vector
  AtomPtr spAtom1(m_atomList[origin_bond_id - 1]);
  AtomPtr spAtom2(m_atomList[target_bond_id - 1]);
  BondPtr spBond(new Bond(bond_id, spAtom1, spAtom2, bond_order));
  m_bondList.push_back(spBond);
  LOG_F(1, "MOL2FileSource::ParseRecordBOND: {}", *spBond);
}

void MOL2FileSource::ParseRecordSUBSTRUCTURE(const std::string &aLine) {
  std::vector<std::string> tokens;
  Tokenize(aLine, tokens);
  // there must be at least 3 fields, rest is optional
  if (tokens.size() < 3)
    throw FileParseError(
        _WHERE_, "Corrupted MOL2 file: not enough fields in SUBSTRUCTURE");

  // Compulsory fields
  int subst_id = std::atoi(tokens[0].c_str());
  std::string subst_name = tokens[1];
  int root_atom = std::atoi(tokens[2].c_str());

  // Optional fields
  std::string chain =
      ((tokens.size() > 5) && (tokens[5] != "****")) ? tokens[5] : "UNK";
  std::string sub_type =
      ((tokens.size() > 6) && (tokens[6] != "****")) ? tokens[6] : "UNK";

  MOL2Substructure ss(subst_name, root_atom, chain, sub_type);
  m_ssInfo[subst_id] = ss;
  LOG_F(1, "MOL2FileSource::ParseRecordSUBSTRUCTURE: {},{},{},{},{}", subst_id,
        subst_name, root_atom, chain, sub_type);
}

void MOL2FileSource::ParseRecordUNSUPPORTED(const std::string &aLine) {
  //	do nothing...
}

std::string MOL2FileSource::GetMOL2Tag(const std::string &aLine) {
  // template based tokenizer is at boost.org, but we have to wait until
  // they are ready to provide Makefile instead of the crap jam they have
  // so, look for a line starting with "@<TRIPOS>"
  if (aLine.find(_TRIPOS_DELIM) != std::string::npos) {
    // tokenize remainder line in case of garbage at the rest of line
    std::vector<std::string> tokens;
    Tokenize(aLine, tokens);
    // get rid of the "@<TRIPOS>" part
    if (tokens[0].size() <= _TRIPOS_DELIM_SIZE) {
      LOG_F(ERROR, "Corrupted MOL2: no tag identifier after @<TRIPOS>.");
      // why the heck is a segfault here?
      throw FileParseError(
          _WHERE_, "Corrupted MOL2: no tag identifier after @<TRIPOS>.");
    } else {
      tokens[0].erase(0, _TRIPOS_DELIM_SIZE);
      return tokens[0]; // return with the very first string
    }
  } else
    // return empty string if it is not a MOL2 tag
    return "";
}

void MOL2FileSource::ParseCountFields(const std::string &aLine) {
  std::vector<std::string> tokens;
  Tokenize(aLine, tokens);
  if (tokens.size() < 1)
    throw FileParseError(_WHERE_,
                         "Corrupted MOL2: no atom numbers in MOLECULE record.");
  // number of atoms:
  nAtoms = std::atoi(tokens[0].c_str());
  if (tokens.size() > 1) // we have bonds number as well
    nBonds = std::atoi(tokens[1].c_str());
  if (tokens.size() > 2) // substructures
    nSubstructures = std::atoi(tokens[2].c_str());
  if (tokens.size() > 3)
    nFeatures = std::atoi(tokens[3].c_str());
  if (tokens.size() > 4) // substructures
    nSets = std::atoi(tokens[4].c_str());
  LOG_F(1, "MOL2FileSource::ParseCountFields: nAtoms {}", nAtoms);
  LOG_F(1, "MOL2FileSource::ParseCountFields: nBonds {}", nBonds);
  LOG_F(1, "MOL2FileSource::ParseCountFields: nSubstructures {}",
        nSubstructures);
  LOG_F(1, "MOL2FileSource::ParseCountFields: nFeatures {}", nFeatures);
  LOG_F(1, "MOL2FileSource::ParseCountFields: nSets {}", nSets);
}

// Correct all the atom attributes that could not be set explicitly when
// parsing the atom records
void MOL2FileSource::SetupAtomParams() {
  FixImplicitHydrogenCount();
  if (m_bImplHydrogens) {
    RemoveNonPolarHydrogens();
    // Quick fix: remove all "Lone Pair" atoms as well
    AtomList lpList = GetAtomListWithPredicate(m_atomList, isFFType_eq("LP"));
    for (AtomListIter iter = lpList.begin(); iter != lpList.end(); iter++) {
      LOG_F(1, "MOL2FileSource::SetupAtomParams: Removing Lone Pair {}",
            (*iter)->GetFullAtomName());
      RemoveAtom(*iter);
    }
    RenumberAtomsAndBonds();
  }
  FixHybridState();
  FixTriposTypes();
  SetupVdWRadii();
  RenumberAtomsAndBonds();
}

// We don't have the benefit of accurate bond orders (as we haven't bothered to
// separate the aromatic bonds into alternate double/single formal bonds) but we
// do have accurate hybridisation states from the Tripos type lookup. Can detect
// implicit hydrogens in the input file by a simple check of the number of
// connections to each carbon
void MOL2FileSource::FixImplicitHydrogenCount() {
  // Get list of all carbons
  AtomList cList = GetAtomListWithPredicate(m_atomList, isAtomicNo_eq(6));
  for (AtomListIter iter = cList.begin(); iter != cList.end(); iter++) {
    Atom::eHybridState hyb = (*iter)->GetHybridState();
    // Compare actual and expected bond count based on hybridisation state
    // Asssume SP3 (4 bonds) if hyb state is undefined
    int actual = (*iter)->GetNumBonds();
    int expected =
        (hyb == Atom::SP3)
            ? 4
            : (hyb == Atom::SP2)
                  ? 3
                  : (hyb == Atom::AROM) ? 3 : (hyb == Atom::SP) ? 2 : 4;
    int nImplH = expected - actual;
    // Pick up fragmented aromatic rings here
    // Cannot have more than one explicit or implicit hydrogen to an aromatic
    // carbon
    if ((hyb == Atom::AROM) && (nImplH > 0)) {
      int nExplH = (*iter)->GetCoordinationNumber(1);
      if ((nImplH + nExplH) > 1) {
        LOG_F(INFO,
              "MOL2FileSource::FixImplicitHydrogenCount: Too few bonds "
              "detected to aromatic carbon {}; #hydrogens capped at 1",
              (*iter)->GetFullAtomName());
        nImplH = 1 - nExplH;
      }
    }
    if (nImplH > 0) {
      (*iter)->SetNumImplicitHydrogens(nImplH);
    } else if (nImplH < 0) {
      LOG_F(INFO,
            "MOL2FileSource::FixImplicitHydrogenCount: Too many bonds detected "
            "to {}",
            (*iter)->GetFullAtomName());
    }
  }
}

// Convert {O,S}_SP3 to TRI, if and only if the atom is bonded to pi-atom
void MOL2FileSource::FixHybridState() {
  AtomList sp3AtomList =
      GetAtomListWithPredicate(m_atomList, isHybridState_eq(Atom::SP3));
  sp3AtomList =
      GetAtomListWithPredicate(sp3AtomList, std::not1(isAtomPosCharged()));

  for (AtomListIter iter = sp3AtomList.begin(); iter != sp3AtomList.end();
       iter++) {
    AtomList bondedAtomList;
    switch ((*iter)->GetAtomicNo()) {
    case 8:  // O
    case 16: // S
      bondedAtomList = GetBondedAtomList(*iter);
      if (FindAtomInList(bondedAtomList, isPiAtom()) != bondedAtomList.end()) {
        (*iter)->SetHybridState(Atom::TRI);
        LOG_F(1, "Changing {} from SP3 to TRI", (*iter)->GetFullAtomName());
      }
      break;
    default:
      break;
    }
  }
}

// Correct the Tripos types read from file, for extended carbons and polar
// hydrogens
void MOL2FileSource::FixTriposTypes() {
  AtomList cList = GetAtomListWithPredicate(m_atomList, isAtomicNo_eq(6));
  for (AtomListIter iter = cList.begin(); iter != cList.end(); iter++) {
    TriposAtomType::eType t = (*iter)->GetTriposType();
    TriposAtomType::eType xt = t;
    int nImplH = (*iter)->GetNumImplicitHydrogens();
    if (nImplH > 0) {
      switch (t) {
      case TriposAtomType::C_3:
        xt = (nImplH == 3)
                 ? TriposAtomType::C_3_H3
                 : (nImplH == 2) ? TriposAtomType::C_3_H2
                                 : (nImplH == 1) ? TriposAtomType::C_3_H1
                                                 : TriposAtomType::C_3;
        break;
      case TriposAtomType::C_2:
        xt = (nImplH == 2)
                 ? TriposAtomType::C_2_H2
                 : (nImplH == 1) ? TriposAtomType::C_2_H1 : TriposAtomType::C_2;
        break;
      case TriposAtomType::C_ar:
        xt = (nImplH == 1) ? TriposAtomType::C_ar_H1 : TriposAtomType::C_ar;
        break;
      case TriposAtomType::C_1:
        xt = (nImplH == 1) ? TriposAtomType::C_1_H1 : TriposAtomType::C_1;
        break;
      default:
        break;
      }
      if (xt != t) {
        LOG_F(1, "Correcting Tripos type for {} from {} to {}",
              (*iter)->GetFullAtomName(), m_typer.Type2Str(t),
              m_typer.Type2Str(xt));
        (*iter)->SetTriposType(xt);
      }
    }
  }
  AtomList hpList = GetAtomListWithPredicate(m_atomList, isAtomHBondDonor());
  for (AtomListIter iter = hpList.begin(); iter != hpList.end(); iter++) {
    TriposAtomType::eType t = (*iter)->GetTriposType();
    TriposAtomType::eType xt = TriposAtomType::H_P;
    if (xt != t) {
      LOG_F(1, "Correcting Tripos type for {} from {} to {}",
            (*iter)->GetFullAtomName(), m_typer.Type2Str(t),
            m_typer.Type2Str(xt));
      (*iter)->SetTriposType(xt);
    }
  }
}

void MOL2FileSource::RemoveNonPolarHydrogens() {
  AtomList cList = GetAtomListWithPredicate(m_atomList, isAtomicNo_eq(6));

  AtomList removeList;
  for (auto &cIter : cList) {
    // Get list of all bonded hydrogens
    AtomList hList =
        GetAtomListWithPredicate(GetBondedAtomList(cIter), isAtomicNo_eq(1));
    int nH = hList.size();
    if (nH > 0) {
      for (auto &hIter : hList) {
        RemoveAtom(hIter);
      }
      // Adjust number of implicit hydrogens
      cIter->SetNumImplicitHydrogens(cIter->GetNumImplicitHydrogens() + nH);
      LOG_F(1,
            "MOL2FileSource::RemoveNonPolarHydrogens: Removing {} hydrogens "
            "from {}",
            nH, cIter->GetFullAtomName());
    }
  }

  // Remove orphan hydrogens (from fragmented residues for e.g.)
  AtomList hList = GetAtomListWithPredicate(m_atomList, isAtomicNo_eq(1));
  hList = GetAtomListWithPredicate(hList, isCoordinationNumber_eq(0));
  for (const auto &hIter : hList) {
    LOG_F(
        1,
        "MOL2FileSource::RemoveNonPolarHydrogens: Removing orphan hydrogen {}",
        hIter->GetFullAtomName());
    RemoveAtom(hIter);
  }
}

// Defines nominal vdW radius (used for cavity mapping)
// correcting for extended atoms and H-bond donor hydrogens
void MOL2FileSource::SetupVdWRadii() {
  // Radius increment for atoms with implicit hydrogens
  // DM 22 Jul 1999 - only increase the radius for sp3 atoms with implicit
  // hydrogens For sp2 and aromatic, leave as is
  isHybridState_eq bIsSP3(Atom::SP3);
  double dImplRadIncr = m_spElementData->GetImplicitRadiusIncr();
  // Element data for hydrogen
  ElementData elHData = m_spElementData->GetElementData(1);
  // Radius increment and predicate for H-bonding hydrogens
  isAtomHBondDonor bIsHBondDonor;
  double dHBondRadius =
      elHData.vdwRadius + m_spElementData->GetHBondRadiusIncr();

  for (auto &iter : m_atomList) {
    // Get the element data for this atom
    int nAtomicNo = iter->GetAtomicNo();
    ElementData elData = m_spElementData->GetElementData(nAtomicNo);
    double vdwRadius = elData.vdwRadius;
    int nImplH = iter->GetNumImplicitHydrogens();
    // Adjust atomic mass and vdw radius for atoms with implicit hydrogens
    if (nImplH > 0) {
      iter->SetAtomicMass(elData.mass +
                          (nImplH * elHData.mass)); // Adjust atomic mass
      if (bIsSP3(iter)) {
        vdwRadius +=
            dImplRadIncr; // adjust vdw radius (for sp3 implicit atoms only)
      }
    }
    // Adjust vdw radius for H-bonding hydrogens
    else if (bIsHBondDonor(iter)) {
      vdwRadius = dHBondRadius;
    }
    // Finally we can set the radius
    iter->SetVdwRadius(vdwRadius);
    LOG_F(1, "{}: #H={}; vdwR={}; mass={}", iter->GetFullAtomName(), nImplH,
          iter->GetVdwRadius(), iter->GetAtomicMass());
  }
}

// Splits the Tripos MOL2 subst_name field into appropriate residue name (sName)
// and residue ID (sID) components If subst_name is **** or in the style <148>,
// then substitutes default name (UNK) and default ID (subst_id)
void MOL2FileSource::GetSSIDandName(const std::string &subst_name, int subst_id,
                                    std::string &sID, std::string &sName) {
  // subst_names beginning with * or < do not follow the normal pattern of ALA99
  // etc
  const std::string nonstandard("*<");
  const std::string numeric("0123456789");
  const std::string defaultName("UNK");
  const std::string defaultID =
      Variant(subst_id); // Quick way to convert from int to string

  if (subst_name.empty() || (subst_name.find_first_of(nonstandard) == 0)) {
    sID = defaultID;
    sName = defaultName;
    return;
  }
  // Find first numeric char
  std::string::size_type nFirst = subst_name.find_first_of(numeric);
  // Split into name and ID:
  // Name is up to the beginning of 1st numeric char
  // ID is from the first numeric char onwards (may have non-numeric suffix as
  // well e.g 100B)
  sID = (nFirst != std::string::npos) ? subst_name.substr(nFirst) : defaultID;
  sName = (nFirst > 0) ? subst_name.substr(0, nFirst) : defaultName;
  return;
}

// simple tokenizer that splits at whitspace
void MOL2FileSource::Tokenize(const std::string &aString,
                              std::vector<std::string> &aTokensBuf) {
  std::string buf;
  std::stringstream ss(aString);
  while (ss >> buf)
    aTokensBuf.push_back(buf);
}
