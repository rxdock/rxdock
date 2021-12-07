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

#include "rxdock/FileError.h"
#include "rxdock/ParameterFileSource.h"

#include <loguru.hpp>

using namespace rxdock;

// Constructors
ParameterFileSource::ParameterFileSource(const char *fileName)
    : BaseFileSource(fileName) {
  std::string name(fileName);
  m_inputFileName = name;
  std::string ext = name.substr(name.size() - 5, 5);
  if (ext == ".json") {
    m_json = true;
  } else {
    m_json = false;
  }
  _RBTOBJECTCOUNTER_CONSTR_("ParameterFileSource");
}

ParameterFileSource::ParameterFileSource(const std::string &fileName)
    : BaseFileSource(fileName) {
  m_inputFileName = fileName;
  std::string ext = fileName.substr(fileName.size() - 5, 5);
  if (ext == ".json") {
    m_json = true;
  } else {
    m_json = false;
  }
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
  // Only parse if we haven't already done so
  if (!m_bParsedOK) {
    const std::string typeKey = "media-type";
    const std::vector<std::string> expectedTypeValues{
        "application/vnd.rxdock.parameters",
        "application/vnd.rxdock.run-script",
        "application/vnd.rxdock.scoring-function"};

    ClearParamsCache(); // Clear current cache
    std::ifstream inputFile(m_inputFileName);
    json params;
    inputFile >> params;
    inputFile.close();

    try {
      // Check the file type
      std::string actualTypeValue;
      params.at(typeKey).get_to(actualTypeValue);
      params.erase(typeKey);
      if (std::find(expectedTypeValues.begin(), expectedTypeValues.end(),
                    actualTypeValue) == expectedTypeValues.end())
        throw FileParseError(
            _WHERE_,
            "Wrong value for " + typeKey + " in " + GetFileName() +
                ", expected any of " +
                fmt::format("{}", fmt::join(expectedTypeValues, ", ")));

      // Load title
      params.at("title").get_to(m_strTitle);
      params.erase("title");
      LOG_F(1, "ParameterFileSource::Parse: Title = {}", m_strTitle);

      // Load version
      params.at("version").get_to(m_strVersion);
      params.erase("version");
      LOG_F(1, "ParameterFileSource::Parse: Version = {}", m_strVersion);

      // Load array of section names
      if (!params.contains("sections"))
        throw FileParseError(_WHERE_, "Missing array of section names in " + GetFileName());

      std::vector<std::string> sectionKeys;
      params.at("sections").get_to(sectionKeys);
      params.erase("sections");
      LOG_F(1, "ParameterFileSource::Parse: Sections = {}", fmt::join(sectionKeys, ", "));

      // Load elements
      for (std::string sectionKey : sectionKeys) {
        if (!params.contains(sectionKey))
          throw FileParseError(_WHERE_, "Wrongly listed section " + sectionKey + " in " + GetFileName());

        json sectionValue = params.at(sectionKey);
        LOG_F(2,
              "ParameterFileSource::Parse: section {} contains parameters {}",
              sectionKey, sectionValue.dump());
        if (sectionKey.rfind("_comment_", 0) == std::string::npos &&
            sectionKey.rfind("_disabled_", 0) == std::string::npos) {
          AddSection(sectionKey);
          for (auto &param : sectionValue.items()) {
            if (param.key().rfind("_comment_", 0) == std::string::npos &&
                param.key().rfind("_previous_", 0) == std::string::npos) {
              std::string paramName = GetFullParameterName(param.key());
              Variant paramVariant;
              json paramValue = param.value();
              if (paramValue.is_boolean()) {
                paramVariant = paramValue.get<bool>();
              } else if (paramValue.is_number()) {
                paramVariant = paramValue.get<double>();
              } else if (paramValue.is_string()) {
                paramVariant = paramValue.get<std::string>();
              } else if (paramValue.is_array()) {
                paramVariant = paramValue.dump();
              } else {
                throw FileParseError(
                    _WHERE_,
                    "Wrong type for " + param.key() + " in " + GetFileName() +
                        ", expected boolean, number, string, or array");
              }
              m_paramsMap[paramName] = paramVariant;
              LOG_F(1, "ParameterFileSource::Parse: {} = {}", paramName,
                    paramValue.dump());
            } else {
              LOG_F(1, "ParameterFileSource::Parse: not parsing {} = {}",
                    param.key(), param.value().dump());
            }
          }
        } else {
          LOG_F(1, "ParameterFileSource::Parse: not parsing {} = {}",
                sectionKey, sectionValue.dump());
        }
        SetSection();
      }

      // If we get this far everything is OK
      m_bParsedOK = true;
    } catch (Error &error) {
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
    throw FileParseError(_WHERE_, "Missing section name in " + GetFileName());
  else if (std::find(m_sectionNames.begin(), m_sectionNames.end(),
                     strSection) != m_sectionNames.end())
    throw FileParseError(_WHERE_, "Duplicate " + strSection +
                                      " section name in " + GetFileName());
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
