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

#include "BaseMolecularFileSource.h"
#include "FileError.h"

#include <functional>

#include <loguru.hpp>

using namespace rxdock;

// Constructors
// BaseMolecularFileSource::BaseMolecularFileSource(const char* fileName,
// const char* sourceName) :
//  BaseFileSource(fileName), BaseMolecularDataSource(sourceName)
//{
//  _RBTOBJECTCOUNTER_CONSTR_("BaseMolecularFileSource");
//}

// Single-record version
BaseMolecularFileSource::BaseMolecularFileSource(const std::string &fileName,
                                                 const std::string &sourceName)
    : BaseFileSource(fileName), BaseMolecularDataSource(sourceName) {
  _RBTOBJECTCOUNTER_CONSTR_("BaseMolecularFileSource");
}

// Multi-record version, with record delimiter passed as an argument
BaseMolecularFileSource::BaseMolecularFileSource(const std::string &fileName,
                                                 const std::string &strRecDelim,
                                                 const std::string &sourceName)
    : BaseFileSource(fileName, strRecDelim),
      BaseMolecularDataSource(sourceName) {
  _RBTOBJECTCOUNTER_CONSTR_("BaseMolecularFileSource");
}

// Destructors
BaseMolecularFileSource::~BaseMolecularFileSource() {
  ClearMolCache(); // Tidy up the atom, bond, title and segment lists
  _RBTOBJECTCOUNTER_DESTR_("BaseMolecularFileSource");
}

// Public methods from BaseMolecularDataSource

// Reset source so that next time it is used, a new set of atom and bond objects
// are created
void BaseMolecularFileSource::Reset() { ClearMolCache(); }

int BaseMolecularFileSource::GetNumTitles() {
  if (isTitleListSupported()) {
    Parse();
    return m_titleList.size();
  } else
    throw InvalidRequest(_WHERE_, "Title list not supported by " + GetName());
}

int BaseMolecularFileSource::GetNumAtoms() {
  if (isAtomListSupported()) {
    Parse();
    if (isSegmentFilterMapDefined())
      return GetNumAtomsWithFilter();
    else
      return m_atomList.size();
  } else
    throw InvalidRequest(_WHERE_, "Atom list not supported by " + GetName());
}

int BaseMolecularFileSource::GetNumBonds() {
  if (isBondListSupported()) {
    Parse();
    if (isSegmentFilterMapDefined())
      return GetNumBondsWithFilter();
    else
      return m_bondList.size();
  } else
    throw InvalidRequest(_WHERE_, "Bond list not supported by " + GetName());
}

int BaseMolecularFileSource::GetNumSegments() {
  if (isAtomListSupported()) {
    Parse();
    return m_segmentMap.size();
  } else
    throw InvalidRequest(_WHERE_, "Atom list not supported by " + GetName());
}

std::vector<std::string> BaseMolecularFileSource::GetTitleList() {
  if (isTitleListSupported()) {
    Parse();
    return m_titleList;
  } else
    throw InvalidRequest(_WHERE_, "Title list not supported by " + GetName());
}

AtomList BaseMolecularFileSource::GetAtomList() {
  if (isAtomListSupported()) {
    Parse();
    if (isSegmentFilterMapDefined())
      return GetAtomListWithFilter();
    else
      return m_atomList;
  } else
    throw InvalidRequest(_WHERE_, "Atom list not supported by " + GetName());
}

BondList BaseMolecularFileSource::GetBondList() {
  if (isBondListSupported()) {
    Parse();
    if (isSegmentFilterMapDefined())
      return GetBondListWithFilter();
    else
      return m_bondList;
  } else
    throw InvalidRequest(_WHERE_, "Bond list not supported by " + GetName());
}

SegmentMap BaseMolecularFileSource::GetSegmentMap() {
  if (isAtomListSupported()) {
    Parse();
    return m_segmentMap;
  } else
    throw InvalidRequest(_WHERE_, "Atom list not supported by " + GetName());
}

// DM 12 May 1999 - support for data records (e.g. SD file)
// Get number of data fields
int BaseMolecularFileSource::GetNumData() {
  if (isDataSupported()) {
    Parse();
    return m_dataMap.size();
  } else
    throw InvalidRequest(_WHERE_, "Data records not supported by " + GetName());
}

// Get list of field names as string list
std::vector<std::string> BaseMolecularFileSource::GetDataFieldList() {
  if (isDataSupported()) {
    Parse();
    std::vector<std::string> dataFieldList;
    dataFieldList.reserve(m_dataMap.size());
    for (StringVariantMapConstIter iter = m_dataMap.begin();
         iter != m_dataMap.end(); iter++)
      dataFieldList.push_back((*iter).first);
    return dataFieldList;
  } else
    throw InvalidRequest(_WHERE_, "Data records not supported by " + GetName());
}

// Get all data as map of key=field name, value=variant (double,string or string
// list)
StringVariantMap BaseMolecularFileSource::GetDataMap() {
  if (isDataSupported()) {
    Parse();
    return m_dataMap;
  } else
    throw InvalidRequest(_WHERE_, "Data records not supported by " + GetName());
}

// Query as to whether a particular data field name is present
bool BaseMolecularFileSource::isDataFieldPresent(
    const std::string &strDataField) {
  if (isDataSupported()) {
    Parse();
    return m_dataMap.find(strDataField) != m_dataMap.end();
  } else
    throw InvalidRequest(_WHERE_, "Data records not supported by " + GetName());
}

// Get a particular data value
Variant BaseMolecularFileSource::GetDataValue(const std::string &strDataField) {
  if (isDataSupported()) {
    Parse();
    StringVariantMapConstIter iter = m_dataMap.find(strDataField);
    if (iter != m_dataMap.end())
      return (*iter).second;
    else
      throw FileMissingParameter(
          _WHERE_, strDataField + " data field not found in " + GetFileName());
  } else
    throw InvalidRequest(_WHERE_, "Data records not supported by " + GetName());
}

// These methods allow filtering of the data source by segment name
// So for example, we could just read the segment named TAR from the source
// The filter modifies the behaviour of getNumAtoms(), GetAtomList,
// GetNumBonds(), GetBondList()
SegmentMap BaseMolecularFileSource::GetSegmentFilterMap() {
  return m_segmentFilterMap;
}

void BaseMolecularFileSource::SetSegmentFilterMap(
    const SegmentMap &segmentFilterMap) {
  m_segmentFilterMap = segmentFilterMap;
}

void BaseMolecularFileSource::ClearSegmentFilterMap() {
  m_segmentFilterMap.clear();
}

bool BaseMolecularFileSource::isSegmentFilterMapDefined() {
  return !m_segmentFilterMap.empty();
}

////////////////////////////////////////////////////////////////////
// Protected methods

// Clear the atom, bond, title and segment lists
void BaseMolecularFileSource::ClearMolCache() {
  m_titleList.clear();
  m_atomList.clear();
  m_bondList.clear();
  m_segmentMap.clear();
  m_dataMap.clear();
  m_bParsedOK = false;
}

// DM 16 Feb 2000 - formerly in MdlFileSource
// Removes an atom and all bonds to the atom, from the atom and bond lists
// Does not adjust the attributes of the atoms whose bonds have been broken
// DM 30 Oct 2000 - no longer renumbers the atoms and bonds
// RenumberAtomsAndBonds should be called after all atoms have been removed
void BaseMolecularFileSource::RemoveAtom(AtomPtr spAtom) {
  // DM 31 Oct 2000 - return by const ref
  const BondMap bondMap = spAtom->GetBondMap();

  // First remove all bonds from the atom
  for (const auto &mapIter : bondMap) {
    auto bIter = FindBond(m_bondList, isBond_eq((mapIter).first));
    if (bIter != m_bondList.end()) {
      LOG_F(1, "Removing bond #{} ({}-{})", (*bIter)->GetBondId(),
            (*bIter)->GetAtom1Ptr()->GetName(),
            (*bIter)->GetAtom2Ptr()->GetName());
      m_bondList.erase(bIter); // Erase the bond
    }
  }

  // Now we have an isolated atom we can remove it
  // Find the atom in the FileSource atom list
  // DM 2 Aug 1999 - search by atom  attributes, not by memory location (v
  // risky) auto aIter =
  // FindAtomInList(m_atomList, std::bind(isAtomPtr_eq(),
  // std::placeholders::_1, spAtom));
  auto aIter = FindAtomInList(
      m_atomList, std::bind(isAtom_eq(), std::placeholders::_1, spAtom));
  if (aIter != m_atomList.end()) {
    LOG_F(1, "Removing atom #{}, {}", (*aIter)->GetAtomId(),
          (*aIter)->GetName());
    m_atomList.erase(aIter); // Erase the atom
  }
}

// DM 30 Oct 2000 - now independent from RemoveAtom
void BaseMolecularFileSource::RenumberAtomsAndBonds() {
  // Renumber the bond and atom IDs from 1
  int nAtomId(1);
  for (AtomListIter aIter = m_atomList.begin(); aIter != m_atomList.end();
       aIter++, nAtomId++) {
    (*aIter)->SetAtomId(nAtomId);
  }
  int nBondId(1);
  for (BondListIter bIter = m_bondList.begin(); bIter != m_bondList.end();
       bIter++, nBondId++) {
    (*bIter)->SetBondId(nBondId);
  }
}

////////////////////////////////////////////////////////////////////
// Private methods

// This version of GetNumAtoms applies the segment filter
// Called by the public GetNumAtoms
// Assume that Parse has already been called
int BaseMolecularFileSource::GetNumAtomsWithFilter() {
  int nAtoms(0);
  // The segment map data member (m_segmentMap) already contains the atom counts
  // for each segment present but we need to iterate over the segment filter map
  // to see how many of those requested are actually present
  for (SegmentMapIter iter = m_segmentFilterMap.begin();
       iter != m_segmentFilterMap.end(); iter++) {
    SegmentMapIter match =
        m_segmentMap.find((*iter).first); // try to find a match
    if (match != m_segmentMap.end())
      nAtoms += (*match).second; // The value (.second) of the map element is
                                 // the number of atoms in the segment
  }
  return nAtoms;
}

// This version of GetNumBonds applies the segment filter
// Called by the public GetNumBonds
// Assume that Parse has already been called
int BaseMolecularFileSource::GetNumBondsWithFilter() {
  int nBonds(0);
  // More tricky this one, there's no option but to iterate over all the bonds
  // to find if the atoms are in one of the segments in the filter map
  for (BondListIter iter = m_bondList.begin(); iter != m_bondList.end();
       iter++) {
    if ((m_segmentFilterMap.find((*iter)->GetAtom1Ptr()->GetSegmentName()) !=
         m_segmentFilterMap.end()) &&
        (m_segmentFilterMap.find((*iter)->GetAtom2Ptr()->GetSegmentName()) !=
         m_segmentFilterMap.end()))
      nBonds++;
  }
  return nBonds;
}

// This version of GetAtomList applies the segment filter
// Called by the public GetAtomList
// Assume that Parse has already been called
AtomList BaseMolecularFileSource::GetAtomListWithFilter() {
  AtomList atomList;
  // Loop over all atoms and check whether they belong to one of the requested
  // segments DM 29 Jul 1999 - set consecutive atom ID's starting from 1,
  // irrespective of ID in file
  int nAtomId(0);
  for (AtomListIter iter = m_atomList.begin(); iter != m_atomList.end();
       iter++) {
    if (m_segmentFilterMap.find((*iter)->GetSegmentName()) !=
        m_segmentFilterMap.end()) {
      (*iter)->SetAtomId(++nAtomId);
      atomList.push_back(*iter); // Got a match, so store the atom pointer
    }
  }
  return atomList;
}

// This version of GetBondList applies the segment filter
// Called by the public GetBondList
// Assume that Parse has already been called
BondList BaseMolecularFileSource::GetBondListWithFilter() {
  BondList bondList;
  // Loop over all bonds and check whether their atoms  belong to one of the
  // requested segments
  for (BondListIter iter = m_bondList.begin(); iter != m_bondList.end();
       iter++) {
    if ((m_segmentFilterMap.find((*iter)->GetAtom1Ptr()->GetSegmentName()) !=
         m_segmentFilterMap.end()) &&
        (m_segmentFilterMap.find((*iter)->GetAtom2Ptr()->GetSegmentName()) !=
         m_segmentFilterMap.end()))
      bondList.push_back(*iter); // Got a match, so store the bond pointer;
  }
  return bondList;
}

// Moved from MOL2FileSource - sets partial charges on a per-residue basis
// Automatically determines protonation states of ASP/GLU/HIS using Mandatory
// and Forbidden atoms in spParamSource
void BaseMolecularFileSource::SetupPartialIonicGroups(
    AtomList &atoms, ParameterFileSourcePtr spParamSource) {
  const std::string _MANDATORY("MANDATORY");
  const std::string _FORBIDDEN("FORBIDDEN");
  if (atoms.empty()) {
    LOG_F(WARNING, "SetupPartialIonicGroups: Empty atom list");
    return;
  }
  AtomPtr leadAtom = atoms.front();
  std::string subunitName = leadAtom->GetSubunitName();
  std::string match = leadAtom->GetSegmentName() + ":" + subunitName + "_" +
                      leadAtom->GetSubunitId() + ":";
  if (GetNumMatchingAtoms(atoms, match) != atoms.size()) {
    LOG_F(WARNING,
          "SetupPartialIonicGroups: Inconsistent subunit names in "
          "atom list headed by {}",
          leadAtom->GetFullAtomName());
    return;
  }
  std::vector<std::string> resList(spParamSource->GetSectionList());
  if (std::find(resList.begin(), resList.end(), subunitName) == resList.end()) {
    LOG_F(INFO, "SetupPartialIonicGroups: No section for residue {}",
          subunitName);
    return;
  }

  spParamSource->SetSection(subunitName);
  std::vector<std::string> atList = spParamSource->GetParameterList();
  if (std::find(atList.begin(), atList.end(), _MANDATORY) != atList.end()) {
    std::string mandatory =
        spParamSource->GetParameterValueAsString(_MANDATORY);
    std::vector<std::string> mandAtoms =
        ConvertDelimitedStringToList(mandatory);
    unsigned int nPresent = GetNumMatchingAtoms(atoms, mandAtoms);
    if (nPresent != mandAtoms.size()) {
      LOG_F(INFO,
            "SetupPartialIonicGroups: Only {} out of {} mandatory atoms "
            "present in atom list headed by {}",
            nPresent, mandAtoms.size(), leadAtom->GetFullAtomName());
      return;
    }
    std::remove(atList.begin(), atList.end(), _MANDATORY);
  }
  if (std::find(atList.begin(), atList.end(), _FORBIDDEN) != atList.end()) {
    std::string forbidden =
        spParamSource->GetParameterValueAsString(_FORBIDDEN);
    std::vector<std::string> forbAtoms =
        ConvertDelimitedStringToList(forbidden);
    int nPresent = GetNumMatchingAtoms(atoms, forbAtoms);
    if (nPresent > 0) {
      LOG_F(INFO,
            "SetupPartialIonicGroups: {} forbidden atoms present in "
            "atom list headed by {}",
            nPresent, leadAtom->GetFullAtomName());
      return;
    }
    std::remove(atList.begin(), atList.end(), _FORBIDDEN);
  }

  for (std::vector<std::string>::const_iterator aIter = atList.begin();
       aIter != atList.end(); aIter++) {
    double partialCharge(spParamSource->GetParameterValue(
        *aIter)); // Get the partial charge value
    LOG_F(1, "SetupPartialIonicGroups: Trying to match {}", *aIter);
    // Find the atoms which match the specifier string
    AtomList selectedAtoms = GetMatchingAtomList(atoms, *aIter);
    // Now we've got the matching atoms, set the group charge on each atom
    for (AtomListIter iter = selectedAtoms.begin(); iter != selectedAtoms.end();
         iter++) {
      (*iter)->SetGroupCharge(partialCharge);
      LOG_F(INFO, "Setting group charge on {} to {}",
            (*iter)->GetFullAtomName(), partialCharge);
    }
  }
}
