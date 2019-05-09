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

#include "RbtFileError.h"
#include "RbtParameterFileSource.h"

// Constructors
RbtParameterFileSource::RbtParameterFileSource(const char *fileName)
    : RbtBaseFileSource(fileName) {
  _RBTOBJECTCOUNTER_CONSTR_("RbtParameterFileSource");
}

RbtParameterFileSource::RbtParameterFileSource(const std::string &fileName)
    : RbtBaseFileSource(fileName) {
  _RBTOBJECTCOUNTER_CONSTR_("RbtParameterFileSource");
}

// Destructor
RbtParameterFileSource::~RbtParameterFileSource() {
  ClearParamsCache();
  _RBTOBJECTCOUNTER_DESTR_("RbtParameterFileSource");
}

////////////////////////////////////////
// Public methods
std::string RbtParameterFileSource::GetTitle() {
  Parse();
  return m_strTitle;
}

std::string RbtParameterFileSource::GetVersion() {
  Parse();
  return m_strVersion;
}

// DM 06 June 2000 - limit #parameters to those in current section
unsigned int RbtParameterFileSource::GetNumParameters() {
  Parse();
  return GetParameterList().size();
}

// DM 06 June 2000 - limit parameters to those in current section
RbtStringList RbtParameterFileSource::GetParameterList() {
  Parse();
  RbtStringList paramList;
  std::string prmPrefix(GetSection() + "::");
  int prmPrefixLen = prmPrefix.size();
  for (RbtStringVariantMapConstIter iter = m_paramsMap.begin();
       iter != m_paramsMap.end(); iter++) {
    if (((*iter).first).find(prmPrefix) == 0) {
      paramList.push_back(((*iter).first).substr(prmPrefixLen));
    }
  }
  return paramList;
}

// DM 4 Feb 1999 Get a particular named parameter value as a double
double
RbtParameterFileSource::GetParameterValue(const std::string &strParamName) {
  Parse();
  if (!strParamName.empty()) {
    std::string strFullParamName = GetFullParameterName(strParamName);
    RbtStringVariantMapConstIter iter = m_paramsMap.find(strFullParamName);
    if (iter != m_paramsMap.end())
      return (*iter).second;
    else
      throw RbtFileMissingParameter(_WHERE_, strFullParamName +
                                                 " parameter not found in " +
                                                 GetFileName());
  }
  return 0.0;
}

// DM 12 Feb 1999 Get a particular named parameter value as a string
std::string RbtParameterFileSource::GetParameterValueAsString(
    const std::string &strParamName) {
  Parse();
  if (!strParamName.empty()) {
    std::string strFullParamName = GetFullParameterName(strParamName);
    RbtStringVariantMapConstIter iter = m_paramsMap.find(strFullParamName);
    if (iter != m_paramsMap.end())
      return (*iter).second;
    else
      throw RbtFileMissingParameter(_WHERE_, strFullParamName +
                                                 " parameter not found in " +
                                                 GetFileName());
  }
  return "";
}

// DM 11 Feb 1999 Check if parameter is present
bool RbtParameterFileSource::isParameterPresent(
    const std::string &strParamName) {
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
int RbtParameterFileSource::GetNumSections() {
  Parse();
  return m_sectionNames.size();
}

// List of section names
RbtStringList RbtParameterFileSource::GetSectionList() {
  Parse();
  return m_sectionNames;
}

// Get current section name
// This is essentially just a prefix for each parameter name
// so we don't actually need to parse the file to get/set the section name
std::string RbtParameterFileSource::GetSection() const { return m_strSection; }

// Set current section name
void RbtParameterFileSource::SetSection(const std::string &strSection) {
  m_strSection = strSection;
}

// Private methods

// Pure virtual in RbtBaseFileSource - needs to be defined here
void RbtParameterFileSource::Parse() {
  const std::string strRbtKey = "RBT_PARAMETER_FILE_V1.00";
  const std::string strTitleKey = "TITLE ";
  const std::string strVersionKey = "VERSION ";
  const std::string strSectionKey = "SECTION ";
  const std::string strEndSectionKey = "END_SECTION";

  // Only parse if we haven't already done so
  if (!m_bParsedOK) {
    ClearParamsCache(); // Clear current cache
    Read();             // Read the file

    try {
      RbtFileRecListIter fileIter = m_lineRecs.begin();
      RbtFileRecListIter fileEnd = m_lineRecs.end();

      //////////////////////////////////////////////////////////
      // 1. Check for RBT string on line 1
      if ((*fileIter).find(strRbtKey) != 0)
        throw RbtFileParseError(_WHERE_, "Missing " + strRbtKey +
                                             " string in " + GetFileName());

      //////////////////////////////////////////////////////////
      // 2. Parse the rest of file
      while (++fileIter != fileEnd) {
        // Ignore blank lines and comment lines
        if (((*fileIter).length() == 0) || ((*fileIter).at(0) == '#')) {
#ifdef _DEBUG
          // std::cout << "Comment record" << std::endl;
#endif //_DEBUG
          continue;
        }
        // Check for Title record
        else if ((*fileIter).find(strTitleKey) == 0) {
          m_strTitle = *fileIter;
          m_strTitle.erase(0, strTitleKey.length());
#ifdef _DEBUG
          // std::cout << "Title = " << m_strTitle << std::endl;
#endif //_DEBUG
        }
        // Check for Version record
        else if ((*fileIter).find(strVersionKey) == 0) {
          m_strVersion = *fileIter;
          m_strVersion.erase(0, strVersionKey.length());
#ifdef _DEBUG
          // std::cout << "Version = " << m_strVersion << std::endl;
#endif //_DEBUG
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
          m_paramsMap[strParamName] = RbtVariant(strParamValue);
#ifdef _DEBUG
          // std::cout << strParamName<< " = " << dParamValue << std::endl;
#endif //_DEBUG
        }
      }
      //////////////////////////////////////////////////////////
      // If we get this far everything is OK
      m_bParsedOK = true;
      SetSection(); // Reset to the unnamed section, in case the final
                    // END_SECTION record is missing
    }

    catch (RbtError &error) {
      ClearParamsCache();
      throw; // Rethrow the RbtError
    }
  }
}

void RbtParameterFileSource::ClearParamsCache() {
  m_paramsMap.clear();
  m_sectionNames.clear();
  m_strTitle = "";
  m_strVersion = "";
  m_strSection = "";
}

// Add a new section name
void RbtParameterFileSource::AddSection(const std::string &strSection) {
  // Only add if the name is not already in the list, and is not the empty
  // string
  if (strSection.empty())
    throw RbtFileParseError(_WHERE_,
                            "Missing SECTION name in " + GetFileName());
  else if (std::find(m_sectionNames.begin(), m_sectionNames.end(),
                     strSection) != m_sectionNames.end())
    throw RbtFileParseError(_WHERE_, "Duplicate " + strSection +
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
RbtParameterFileSource::GetFullParameterName(const std::string &strParamName) {
  // Already has section name present, so return unchanged
  if (strParamName.find("::") != std::string::npos)
    return strParamName;
  // Else adorn with current section name
  else
    return GetSection() + "::" + strParamName;
}
