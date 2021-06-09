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

// File reader for rxdock parameter files

#ifndef _RBTPARAMETERFILESOURCE_H_
#define _RBTPARAMETERFILESOURCE_H_

#include "rxdock/BaseFileSource.h"
#include "rxdock/Variant.h"

namespace rxdock {

class ParameterFileSource : public BaseFileSource {
public:
  // Constructors
  ParameterFileSource(const char *fileName);
  RBTDLL_EXPORT ParameterFileSource(const std::string &fileName);

  // Destructor
  virtual ~ParameterFileSource();

  ////////////////////////////////////////
  // Public methods
  RBTDLL_EXPORT std::string GetTitle();
  std::string GetVersion();
  // DM 06 June 2000 - limit #parameters to those in current section
  unsigned int GetNumParameters();
  // DM 4 Feb 1999 Renamed from GetParameters()
  // std::map<std::string, double> GetParameterMap();
  // DM 12 Feb 1999 - only return the list of parameter names, not their values
  // DM 06 Jun 2000 - limits params to those in named section
  RBTDLL_EXPORT std::vector<std::string> GetParameterList();
  // DM 4 Feb 1999 Get a particular named parameter value as a double
  double GetParameterValue(const std::string &strParamName);
  // DM 12 Feb 1999 Get a particular named parameter value as a string
  RBTDLL_EXPORT std::string
  GetParameterValueAsString(const std::string &strParamName);
  // DM 11 Feb 1999 Check if parameter is present
  bool isParameterPresent(const std::string &strParamName);

  // DM 11 Feb 1999 - section handling
  // Parameters can be grouped into named sections
  // such that the same parameter name can appear in multiple sections
  // Default namespace is the unnamed section.
  // Main use is for simulation protocols which may need to define a variable
  // number of phases - e.g. high temperature sampling, cooling phase,
  // minimisation phase and need the same parameters to appear in each
  int GetNumSections();                      // Number of named sections
  std::vector<std::string> GetSectionList(); // List of section names
  std::string GetSection() const;            // Get current section name
  RBTDLL_EXPORT void
  SetSection(const std::string &strSection = ""); // Set current section name

protected:
  // Protected methods

private:
  // Private methods
  ParameterFileSource(); // Disable default constructor
  ParameterFileSource(
      const ParameterFileSource &); // Copy constructor disabled by default
  ParameterFileSource &
  operator=(const ParameterFileSource &); // Copy assignment disabled by default

  // Pure virtual in BaseFileSource - needs to be defined here
  virtual void Parse();
  void ClearParamsCache();
  // Add a new section name
  // Throws an error if section name is empty, or is a duplicate
  void AddSection(const std::string &strSection);
  // Returns the fully qualified parameter name (<section>::<parameter name>)
  // Checks if name already contains a section name, if so just returns the name
  // unchanged If not, prefixes the name with the current section name
  std::string GetFullParameterName(const std::string &strParamName);

protected:
  // Protected data

private:
  // Private data
  std::string m_inputFileName;
  std::string m_strTitle;
  std::string m_strVersion;
  StringVariantMap m_paramsMap;
  std::vector<std::string> m_sectionNames; // List(vector) of section names
  std::string m_strSection;                // Current section
  bool m_json;
};

// useful typedefs
typedef SmartPtr<ParameterFileSource> ParameterFileSourcePtr; // Smart pointer

} // namespace rxdock

#endif //_RBTPARAMETERFILESOURCE_H_
