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

#include "rxdock/SiteMapperFactory.h"

#include "rxdock/FileError.h"
#include "rxdock/LigandSiteMapper.h"
#include "rxdock/SphereSiteMapper.h"

using namespace rxdock;

// Parameter name which identifies a scoring function definition
const std::string SiteMapperFactory::_MAPPER = "site-mapper";

////////////////////////////////////////
// Constructors/destructors
SiteMapperFactory::SiteMapperFactory() {}

SiteMapperFactory::~SiteMapperFactory() {}

////////////////////////////////////////
// Public methods
////////////////
// Creates a single site mapper object of named class
SiteMapper *SiteMapperFactory::Create(const std::string &_strMapperClass,
                                      const std::string &strName) {
  std::string strMapperClass = _strMapperClass;
  // compatibility with old format Rbt*
  if (strMapperClass.substr(0, 3) == "Rbt")
    strMapperClass = _strMapperClass.substr(3, std::string::npos);

  if (strMapperClass == SphereSiteMapper::_CT)
    return new SphereSiteMapper(strName);
  if (strMapperClass == LigandSiteMapper::_CT)
    return new LigandSiteMapper(strName);

  throw BadArgument(_WHERE_, "Unknown site mapper " + strMapperClass);
}

// Creates a single site mapper object whose class is defined
// as the value of the SITE_MAPPER parameter in the strName section of the
// parameter file source argument. Also sets the site mapper parameters from the
// remaining parameter values in the current section
// Note: the current section is restored to its original value upon exit
SiteMapper *
SiteMapperFactory::CreateFromFile(ParameterFileSourcePtr spPrmSource,
                                  const std::string &strName) {
  std::string strOrigSection(spPrmSource->GetSection());
  spPrmSource->SetSection(strName);
  if (spPrmSource->isParameterPresent(_MAPPER)) {
    std::string strMapperClass(spPrmSource->GetParameterValueAsString(_MAPPER));
    // Create new site mapper according to the string value of _MAPPER parameter
    SiteMapper *pSiteMapper = Create(strMapperClass, strName);
    // Set all the mapper parameters from the rest of the parameters listed
    std::vector<std::string> prmList = spPrmSource->GetParameterList();
    for (std::vector<std::string>::const_iterator prmIter = prmList.begin();
         prmIter != prmList.end(); prmIter++) {
      if ((*prmIter) != _MAPPER) { // Skip _SF parameter
        pSiteMapper->SetParameter(
            *prmIter, spPrmSource->GetParameterValueAsString(*prmIter));
      }
    }
    spPrmSource->SetSection(strOrigSection);
    return pSiteMapper;
  } else {
    spPrmSource->SetSection(strOrigSection);
    throw FileMissingParameter(_WHERE_, "Missing " + _MAPPER +
                                            " parameter in section " + strName);
  }
}
