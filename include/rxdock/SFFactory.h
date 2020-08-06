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

// Factory class for creating scoring function objects. Scoring functions are
// specified by string names (as defined in each class's static _CT string)

#ifndef _RBTSFFACTORY_H_
#define _RBTSFFACTORY_H_

#include "rxdock/Config.h"
#include "rxdock/ParameterFileSource.h"
#include "rxdock/SFAgg.h"

namespace rxdock {

class SFFactory {
  // Parameter name which identifies a scoring function definition
  static std::string _SF;

public:
  ////////////////////////////////////////
  // Constructors/destructors

  RBTDLL_EXPORT SFFactory(); // Default constructor

  virtual ~SFFactory(); // Default destructor

  ////////////////////////////////////////
  // Public methods
  ////////////////

  // Creates a single scoring function object of type strSFClass, and name
  // strName e.g. strSFClass = HBondIntnSF
  virtual BaseSF *Create(const std::string &strSFClass,
                         const std::string &strName);

  // Creates an aggregate scoring function from a parameter file source
  // Each component SF is in a named section, which should minimally contain a
  // SCORING_FUNCTION parameter whose value is the class name to instantiate
  // strSFClasses contains a comma-delimited list of SF class names to
  // instantiate If strSFClasses is empty, all named sections in spPrmSource are
  // scanned for valid scoring function definitions SF parameters are set from
  // the list of parameters in each named section
  virtual SFAgg *
  CreateAggFromFile(ParameterFileSourcePtr spPrmSource,
                    const std::string &strName,
                    const std::string &strSFClasses = std::string());

protected:
  ////////////////////////////////////////
  // Protected methods
  ///////////////////

private:
  ////////////////////////////////////////
  // Private methods
  /////////////////

  SFFactory(const SFFactory &); // Copy constructor disabled by default

  SFFactory &
  operator=(const SFFactory &); // Copy assignment disabled by default

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
typedef SmartPtr<SFFactory> SFFactoryPtr; // Smart pointer

} // namespace rxdock

#endif //_RBTSFFACTORY_H_
