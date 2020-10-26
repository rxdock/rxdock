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

// Factory class for creating docking site mapper objects.

#ifndef _RBTSITEMAPPERFACTORY_H_
#define _RBTSITEMAPPERFACTORY_H_

#include "rxdock/Config.h"
#include "rxdock/ParameterFileSource.h"
#include "rxdock/SiteMapper.h"

namespace rxdock {

class SiteMapperFactory {
  // Parameter name which identifies a site mapper definition
  static std::string _MAPPER;

public:
  ////////////////////////////////////////
  // Constructors/destructors

  RBTDLL_EXPORT SiteMapperFactory(); // Default constructor
  virtual ~SiteMapperFactory();      // Default destructor

  ////////////////////////////////////////
  // Public methods
  ////////////////

  // Creates a single site mapper object of named class
  virtual SiteMapper *Create(const std::string &_strMapperClass,
                             const std::string &strName);

  // Creates a single site mapper object whose class is defined
  // as the value of the SITE_MAPPER parameter in the strName section of the
  // parameter file source argument. Also sets the site mapper parameters from
  // the remaining parameter values in the current section Note: the current
  // section is restored to its original value upon exit
  virtual SiteMapper *CreateFromFile(ParameterFileSourcePtr spPrmSource,
                                     const std::string &strName);

protected:
  ////////////////////////////////////////
  // Protected methods
  ///////////////////

private:
  ////////////////////////////////////////
  // Private methods
  /////////////////

  SiteMapperFactory(
      const SiteMapperFactory &); // Copy constructor disabled by default
  SiteMapperFactory &
  operator=(const SiteMapperFactory &); // Copy assignment disabled by default

protected:
  ////////////////////////////////////////
  // Protected data
  ////////////////

private:
  ////////////////////////////////////////
  // Private data
  //////////////
};

// Useful typedefs
typedef SmartPtr<SiteMapperFactory> SiteMapperFactoryPtr; // Smart pointer

} // namespace rxdock

#endif //_RBTSITEMAPPERFACTORY_H_
