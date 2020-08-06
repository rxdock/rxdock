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

#include <algorithm> //For std::find

#include "FileError.h"
#include "ParameterFileSource.h"

#include <loguru.hpp>

using namespace rxdock;

// Constructors
ParameterFileSource::ParameterFileSource(const char *fileName)
    : BaseFileSource(fileName) {
  _RBTOBJECTCOUNTER_CONSTR_("ParameterFileSource");
}

ParameterFileSource::ParameterFileSource(const std::string &fileName)
    : BaseFileSource(fileName) {
  _RBTOBJECTCOUNTER_CONSTR_("ParameterFileSource");
}

// Destructor
ParameterFileSource::~ParameterFileSource() {
  ClearParamsCache();
  _RBTOBJECTCOUNTER_DESTR_("ParameterFileSource");
}

////////////////////////////////////////
// Public methods
std::string ParameterFileSource::GetTitle() {
  Parse();
  return m_strTitle;
}

std::string ParameterFileSource::GetVersion() {
  Parse();
  return m_strVersion;
}

// DM 06 June 2000 - limit #parameters to those in current section
unsigned int ParameterFileSource::GetNumParameters() {
  Parse();
  return GetParameterList().size();
}

// DM 06 June 2000 - limit parameters to those in current section
std::vector<std::string> ParameterFileSource::GetParameterList() {
  Parse();
  std::vector<std::string> paramList;
  std::string prmPrefix(GetSection() + "::");
  int prmPrefixLen = prmPrefix.size();
  for (StringVariantMapConstIter iter = m_paramsMap.begin();
       iter != m_paramsMap.end(); iter++) {
    if (((*iter).first).find(prmPrefix) == 0) {
      paramList.push_back(((*iter).first).substr(prmPrefixLen));
    }
  }
  return paramList;
}

// DM 4 Feb 1999 Get a particular named parameter value as a double
double ParameterFileSource::GetParameterValue(const std::string &strParamName) {
  Parse();
  if (!strParamName.empty()) {
    std::string strFullParamName = GetFullParameterName(strParamName);
    StringVariantMapConstIter iter = m_paramsMap.find(strFullParamName);
    if (iter != m_paramsMap.end())
      return (*iter).second;
    else
      throw FileMissingParameter(_WHERE_, strFullParamName +
                                              " parameter not found in " +
                                              GetFileName());
  }
  return 0.0;
}

// DM 12 Feb 1999 Get a particular named parameter value as a string
std::string ParameterFileSource::GetParameterValueAsString(
    const std::string &strParamName) {
  Parse();
  if (!strParamName.empty()) {
    std::string strFullParamName = GetFullParameterName(strParamName);
    StringVariantMapConstIter iter = m_paramsMap.find(strFullParamName);
    if (iter != m_paramsMap.end())
      return (*iter).second;
    else
      throw FileMissingParameter(_WHERE_, strFullParamName +
                                              " parameter not found in " +
                                              GetFileName());
  }
  return "";
}

// DM 11 Feb 1999 Check if parameter is present
bool ParameterFileSource::isParameterPresent(const std::string &strParamName) {
  Parse();
  std::string strFullParamName = GetFullParameterName(strParamName);
  return m_paramsMap.find(strFullParamName) != m_paramsMap.end();
}

// DM 11 Feb 1999 - section handling
// Parameters can be grouped into named sections
// such that the same parameter name can appear in multiple sections
// Default namespace is the unnamed section.
// Main use is for simulation protocols which may need to define a variable
// number of phases - e.g. high temperature sampling, cooling phase,
// minimisation phase and need the same parameters to appear in each

// Number of named sections
int ParameterFileSource::GetNumSections() {
  Parse();
  return m_sectionNames.size();
}

// List of section names
std::vector<std::string> ParameterFileSource::GetSectionList() {
  Parse();
  return m_sectionNames;
}

// Get current section name
// This is essentially just a prefix for each parameter name
// so we don't actually need to parse the file to get/set the section name
std::string ParameterFileSource::GetSection() const { return m_strSection; }

// Set current section name
void ParameterFileSource::SetSection(const std::string &strSection) {
  m_strSection = strSection;
}

// Private methods

// Pure virtual in BaseFileSource - needs to be defined here
void ParameterFileSource::Parse() {
  const std::string strKey = "RBT_PARAMETER_FILE_V1.00";
  const std::string strTitleKey = "TITLE ";
  const std::string strVersionKey = "VERSION ";
  const std::string strSectionKey = "SECTION ";
  const std::string strEndSectionKey = "END_SECTION";

  // Only parse if we haven't already done so
  if (!m_bParsedOK) {
    ClearParamsCache(); // Clear current cache
    Read();             // Read the file

    try {
      FileRecListIter fileIter = m_lineRecs.begin();
      FileRecListIter fileEnd = m_lineRecs.end();

      //////////////////////////////////////////////////////////
      // 1. Check for RBT string on line 1
      if ((*fileIter).find(strKey) != 0)
        throw FileParseError(_WHERE_, "Missing " + strKey + " string in " +
                                          GetFileName());

      //////////////////////////////////////////////////////////
      // 2. Parse the rest of file
      while (++fileIter != fileEnd) {
        // Ignore blank lines and comment lines
        if (((*fileIter).length() == 0) || ((*fileIter).at(0) == '#')) {
          LOG_F(1, "ParameterFileSource::Parse: Comment record");
          continue;
        }
        // Check for Title record
        else if ((*fileIter).find(strTitleKey) == 0) {
          m_strTitle = *fileIter;
          m_strTitle.erase(0, strTitleKey.length());
          LOG_F(1, "ParameterFileSource::Parse: Title = ", m_strTitle);
        }
        // Check for Version record
        else if ((*fileIter).find(strVersionKey) == 0) {
          m_strVersion = *fileIter;
          m_strVersion.erase(0, strVersionKey.length());
          LOG_F(1, "ParameterFileSource::Parse: Version = ", m_strVersion);
        }
        // Check for Section record
        else if ((*fileIter).find(strSectionKey) == 0) {
          std::string strSection;
          std::string strDummy;
          std::istringstream istr(*fileIter);
          istr >> strDummy >> strSection;
          AddSection(strSection);
        }
        // Check for End Section record
        else if ((*fileIter).find(strEndSectionKey) == 0) {
          // Just have to set an empty section name to return to the
          // default unnamed section
          SetSection();
        }
        // Assume everything else is a Key=Parameter name, Value=parameter value
        // pair
        else {
          std::string strParamName;
          // DM 12 May 1999 - read as string then convert to variant
          std::string strParamValue;
          std::istringstream istr(*fileIter);
          istr >> strParamName >> strParamValue;
          // Prefix the parameter name with the section name and ::
          // Hopefully, this will ensure unique parameter names between sections
          strParamName = GetFullParameterName(strParamName);
          m_paramsMap[strParamName] = Variant(strParamValue);
          LOG_F(1, "ParameterFileSource::Parse: {} = {}", strParamName,
                strParamValue);
        }
      }
      //////////////////////////////////////////////////////////
      // If we get this far everything is OK
      m_bParsedOK = true;
      SetSection(); // Reset to the unnamed section, in case the final
                    // END_SECTION record is missing
    }

    catch (Error &error) {
      ClearParamsCache();
      throw; // Rethrow the Error
    }
  }
}

void ParameterFileSource::ClearParamsCache() {
  m_paramsMap.clear();
  m_sectionNames.clear();
  m_strTitle = "";
  m_strVersion = "";
  m_strSection = "";
}

// Add a new section name
void ParameterFileSource::AddSection(const std::string &strSection) {
  // Only add if the name is not already in the list, and is not the empty
  // string
  if (strSection.empty())
    throw FileParseError(_WHERE_, "Missing SECTION name in " + GetFileName());
  else if (std::find(m_sectionNames.begin(), m_sectionNames.end(),
                     strSection) != m_sectionNames.end())
    throw FileParseError(_WHERE_, "Duplicate " + strSection +
                                      " SECTION name in " + GetFileName());
  else {
    m_sectionNames.push_back(strSection);
    SetSection(strSection);
  }
}

// Returns the fully qualified parameter name (<section>::<parameter name>)
// Checks if name already contains a section name, if so just returns the name
// unchanged If not, prefixes the name with the current section name
std::string
ParameterFileSource::GetFullParameterName(const std::string &strParamName) {
  // Already has section name present, so return unchanged
  if (strParamName.find("::") != std::string::npos)
    return strParamName;
  // Else adorn with current section name
  else
    return GetSection() + "::" + strParamName;
}
