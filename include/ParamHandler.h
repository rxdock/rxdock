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

// Base implementation class for managing a collection of named parameters
// Any class requiring parameter handling can derive from ParamHandler
// Parameters are stored as Variants (double, string or stringlist)
// Only derived classes can add and delete parameters from the collection

#ifndef _RBTPARAMHANDLER_H_
#define _RBTPARAMHANDLER_H_

#include "Config.h"
#include "Variant.h"

namespace rxdock {

class ParamHandler {
public:
  ////////////////////////////////////////
  // Constructors/destructors
  virtual ~ParamHandler(); // Default destructor

  ////////////////////////////////////////
  // Public methods
  ////////////////
  // Get number of stored parameters
  unsigned int GetNumParameters() const;
  // Get a named parameter, throws error if name not found
  Variant GetParameter(const std::string &strName) const;
  // Check if named parameter is present
  bool isParameterValid(const std::string &strName) const;
  // Get list of all parameter names
  std::vector<std::string> GetParameterNames() const;
  // Get list of all parameter
  StringVariantMap GetParameters() const;

  // Set named parameter to new value, throws error if name not found
  RBTDLL_EXPORT void SetParameter(const std::string &strName,
                                  const Variant &vValue);

  // Virtual function for dumping parameters to an output stream
  // Called by operator <<
  virtual void Print(std::ostream &s) const;

protected:
  ////////////////////////////////////////
  // Protected methods
  ///////////////////
  ParamHandler(); // Default constructor

  // Only derived classes can mess with the parameter list
  void AddParameter(const std::string &strName, const Variant &vValue);
  void DeleteParameter(const std::string &strName);
  void ClearParameters();
  // DM 25 Oct 2000 - ParameterUpdated is invoked whenever SetParameter is
  // called to allow derived classes to manage a data member which tracks the
  // param value Useful for performance purposes as there is quite an overhead
  // in finding a string in a map, then converting from a Variant to the native
  // datatype Base class version does nothing
  virtual void ParameterUpdated(const std::string &strName) {}

private:
  ////////////////////////////////////////
  // Private methods
  /////////////////

  ParamHandler(const ParamHandler &); // Copy constructor disabled by default
  ParamHandler &
  operator=(const ParamHandler &); // Copy assignment disabled by default

protected:
  ////////////////////////////////////////
  // Protected data
  ////////////////

private:
  ////////////////////////////////////////
  // Private data
  //////////////
  StringVariantMap m_parameters;
};

////////////////////////////////////////
// Non-member functions

// Insertion operator (primarily for debugging)
RBTDLL_EXPORT std::ostream &operator<<(std::ostream &s, const ParamHandler &ph);

} // namespace rxdock

#endif //_RBTPARAMHANDLER_H_
