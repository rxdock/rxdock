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

#include <sstream>

#include "rxdock/CharmmDataSource.h"
#include "rxdock/CharmmTypesFileSource.h"
#include "rxdock/ElementFileSource.h"
#include "rxdock/FileError.h"

using namespace rxdock;

////////////////////////////////////////
// Constructors/destructors

// DM 30 Apr 1999 - pass in masses.rtf file name as parameter
// (default=GetFileName("data","masses.rtf"))
CharmmDataSource::CharmmDataSource(const std::string &strMassesFile) {
  // Need element data source to convert element names to atomic numbers
  ElementFileSourcePtr spElementData(
      new ElementFileSource(GetDataFileName("data", "elements.json")));

  // Initialize m_atomTypes, m_atomicNumber and m_hybridState from masses.rtf
  CharmmTypesFileSourcePtr spTypesFile(
      new CharmmTypesFileSource(strMassesFile));
  CharmmTypeList typeList = spTypesFile->GetTypeList();
  for (CharmmTypeListIter iter = typeList.begin(); iter != typeList.end();
       iter++) {
    m_atomTypes[iter->nAtomType] = iter->strAtomType;
    ElementData elData = spElementData->GetElementData(iter->element);
    m_atomicNumber[iter->strAtomType] = elData.atomicNo;
    // The hybridisation state has been encoded in the comments field in the
    // rxdock version of masses.rtf
    m_hybridState[iter->strAtomType] =
        ConvertCommentStringToHybridState(iter->comment);
  }

  // Initialize m_implicitHydrogens
  m_implicitHydrogens["CH1E"] = 1;
  m_implicitHydrogens["CH2E"] = 2;
  m_implicitHydrogens["CH3E"] = 3;
  m_implicitHydrogens["C5RE"] = 1;
  m_implicitHydrogens["C6RE"] = 1;
  m_implicitHydrogens["SH1E"] = 1;
  m_implicitHydrogens["NH2E"] =
      2; // Added DM 4 Jan 1999 - support older Charmm atom types in XPLOR

  // Initialize m_formalCharge
  // Note: these formal charges are not intended to sum to the correct overall
  // charge on the molecule It is only the sign of the charge which is correct,
  // and indicates the atom can take part in a charged interaction
  m_formalCharge["OC"] = -1; // e.g. in carboxylate, phosphate
  m_formalCharge["HC"] = +1; // e.g. in ammonium, guanidinium
  m_formalCharge["NC"] = +1; // e.g. in guanidinium

  _RBTOBJECTCOUNTER_CONSTR_("CharmmDataSource");
}
// Default destructor
CharmmDataSource::~CharmmDataSource() {
  _RBTOBJECTCOUNTER_DESTR_("CharmmDataSource");
}

////////////////////////////////////////
// Public methods
////////////////
std::string CharmmDataSource::AtomTypeString(int nAtomType) {
  std::map<int, std::string>::iterator iter = m_atomTypes.find(nAtomType);
  if (iter != m_atomTypes.end())
    return (*iter).second;
  else {
    std::ostringstream ostr;
    ostr << "CHARMm atom type #" << nAtomType
         << " not found in CharmmDataSource" << std::ends;
    std::string strError(ostr.str());
    // delete ostr.str();
    throw FileMissingParameter(_WHERE_, strError);
  }
}

int CharmmDataSource::ImplicitHydrogens(const std::string &strFFType) {
  std::map<std::string, int>::iterator iter =
      m_implicitHydrogens.find(strFFType);
  if (iter != m_implicitHydrogens.end())
    return (*iter).second;
  else
    return 0;
}

int CharmmDataSource::AtomicNumber(const std::string &strFFType) {
  std::map<std::string, int>::iterator iter = m_atomicNumber.find(strFFType);
  if (iter != m_atomicNumber.end())
    return (*iter).second;
  else
    throw FileMissingParameter(_WHERE_, "CHARMm atom type " + strFFType +
                                            " not found in CharmmDataSource");
}

int CharmmDataSource::FormalCharge(const std::string &strFFType) {
  std::map<std::string, int>::iterator iter = m_formalCharge.find(strFFType);
  if (iter != m_formalCharge.end()) {
    return (*iter).second;
  } else
    return 0;
}

// DM 8 Dec 1998
Atom::eHybridState CharmmDataSource::HybridState(const std::string &strFFType) {
  StringHybridStateMapIter iter = m_hybridState.find(strFFType);
  if (iter != m_hybridState.end())
    return (*iter).second;
  else
    throw FileMissingParameter(_WHERE_, "CHARMm atom type " + strFFType +
                                            " not found in CharmmDataSource");
}

////////////////////////////////////////
// Private methods
////////////////////////////////////////

// DM 8 Dec 1998 Searches for Hybridisation state string in the masses.rtf
// comment field and returns the matching Atom::eHybridState enum Valid
// strings are (RBT::SP), (RBT::SP2), (RBT::SP3), (RBT::AROM), (RBT::TRI)
//(brackets are important)
Atom::eHybridState CharmmDataSource::ConvertCommentStringToHybridState(
    const std::string &strComment) {
  if (strComment.find("(RBT::SP)") != std::string::npos)
    return Atom::SP;
  else if (strComment.find("(RBT::SP2)") != std::string::npos)
    return Atom::SP2;
  else if (strComment.find("(RBT::SP3)") != std::string::npos)
    return Atom::SP3;
  else if (strComment.find("(RBT::AROM)") != std::string::npos)
    return Atom::AROM;
  else if (strComment.find("(RBT::TRI)") != std::string::npos)
    return Atom::TRI;
  else
    return Atom::UNDEFINED;
}

void rxdock::to_json(json &j, const CharmmDataSource &charDatSrc) {
  j = json{{"atom-types", charDatSrc.m_atomTypes},
           {"imp-hydrog", charDatSrc.m_implicitHydrogens},
           {"at-num", charDatSrc.m_atomicNumber},
           {"form-charg", charDatSrc.m_formalCharge},
           {"hyb-state", charDatSrc.m_hybridState}};
}

void rxdock::from_json(const json &j, CharmmDataSource &charDatSrc) {
  j.at("atom-types").get_to(charDatSrc.m_atomTypes);
  j.at("imp-hydrog").get_to(charDatSrc.m_implicitHydrogens);
  j.at("at-num").get_to(charDatSrc.m_atomicNumber);
  j.at("form-charg").get_to(charDatSrc.m_formalCharge);
  j.at("hyb-state").get_to(charDatSrc.m_hybridState);
}