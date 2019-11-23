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

#include "RbtParamHandler.h"

RbtParamHandler::RbtParamHandler() {
  _RBTOBJECTCOUNTER_CONSTR_("RbtParamHandler");
}

RbtParamHandler::~RbtParamHandler() {
  _RBTOBJECTCOUNTER_DESTR_("RbtParamHandler");
}

////////////////////////////////////////
// Public methods
////////////////
// Get number of stored parameters
unsigned int RbtParamHandler::GetNumParameters() const {
  return m_parameters.size();
}

// Get a named parameter, throws error if name not found
RbtVariant RbtParamHandler::GetParameter(const std::string &strName) const {
  RbtStringVariantMapConstIter iter = m_parameters.find(strName);
  if (iter == m_parameters.end()) {
    throw RbtBadArgument(_WHERE_, "Undefined parameter " + strName);
  } else
    return (*iter).second;
}

// Check if named parameter is present
bool RbtParamHandler::isParameterValid(const std::string &strName) const {
  return m_parameters.find(strName) != m_parameters.end();
}

// Get list of all parameter names
std::vector<std::string> RbtParamHandler::GetParameterNames() const {
  std::vector<std::string> nameList;
  for (RbtStringVariantMapConstIter iter = m_parameters.begin();
       iter != m_parameters.end(); iter++) {
    nameList.push_back((*iter).first);
  }
  return nameList;
}

// Get map of all parameters
RbtStringVariantMap RbtParamHandler::GetParameters() const {
  return m_parameters;
}

// Set named parameter to new value, throws error if name not found
void RbtParamHandler::SetParameter(const std::string &strName,
                                   const RbtVariant &vValue) {
  if (!isParameterValid(strName)) {
    throw RbtBadArgument(_WHERE_, "Undefined parameter " + strName);
  } else {
    m_parameters[strName] = vValue;
    // DM 25 Oct 2000 - notify derived class that parameter has changed
    ParameterUpdated(strName);
  }
}

////////////////////////////////////////
// Protected methods
///////////////////
// Only derived classes can mess with the parameter list
void RbtParamHandler::AddParameter(const std::string &strName,
                                   const RbtVariant &vValue) {
  m_parameters[strName] = vValue;
  // DM 25 Oct 2000 - notify derived class that parameter has changed
  // DM 12 Apr 2002 - no need to call ParameterUpdated here
  // Causes problems during construction of derived classes with virtual base
  // classes Better solution is for derived classes to synchronise the initial
  // values of private data members themselves during construction
  // ParameterUpdated(strName);
}

void RbtParamHandler::DeleteParameter(const std::string &strName) {
  m_parameters.erase(strName);
}

void RbtParamHandler::ClearParameters() { m_parameters.clear(); }

// Virtual function for dumping parameters to an output stream
// Called by operator <<
void RbtParamHandler::Print(std::ostream &s) const {
  for (RbtStringVariantMapConstIter iter = m_parameters.begin();
       iter != m_parameters.end(); iter++) {
    std::cout << (*iter).first << "\t" << (*iter).second << std::endl;
  }
}

////////////////////////////////////////
// Non-member functions

// Insertion operator (primarily for debugging)
std::ostream &operator<<(std::ostream &s, const RbtParamHandler &ph) {
  ph.Print(s);
  return s;
}
