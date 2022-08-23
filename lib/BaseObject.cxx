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

#include "rxdock/BaseObject.h"
#include "rxdock/SFRequest.h"
#include "rxdock/WorkSpace.h"

#include <loguru.hpp>

using namespace rxdock;

// Static data members
const std::string BaseObject::_CT = "BaseObject";
const std::string BaseObject::_CLASS = "CLASS";
const std::string BaseObject::_NAME = "NAME";
const std::string BaseObject::_ENABLED = "ENABLED";

////////////////////////////////////////
// Constructors/destructors
BaseObject::BaseObject(const std::string &strClass, const std::string &strName)
    : m_workspace(nullptr), m_enabled(true) {
  LOG_F(2, "BaseObject parameterised constructor for {}", strClass);
  // Add parameters
  AddParameter(_CLASS, strClass);
  AddParameter(_NAME, strName);
  AddParameter(_ENABLED, m_enabled);
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

// Dummy default constructor for virtual base subclasses
// Should never get called
BaseObject::BaseObject() {
  LOG_F(WARNING, "BaseObject default constructor");
  //_RBTOBJECTCOUNTER_CONSTR_(_CT);
}

BaseObject::~BaseObject() {
  LOG_F(2, "BaseObject destructor");
  Unregister(); // Unregister object from workspace
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

////////////////////////////////////////
// Public methods
////////////////

// Class name (e.g. ConstSF)
std::string BaseObject::GetClass() const { return GetParameter(_CLASS); }

std::string BaseObject::GetName() const { return GetParameter(_NAME); }
void BaseObject::SetName(const std::string &strName) {
  SetParameter(_NAME, strName);
}
// Fully qualified name (should be overridden by subclasses which can be
// aggregated to prefix the name with the parent's name)
std::string BaseObject::GetFullName() const { return GetName(); }

void BaseObject::Enable() { SetParameter(_ENABLED, true); }
void BaseObject::Disable() { SetParameter(_ENABLED, false); }
bool BaseObject::isEnabled() const { return m_enabled; }

// WorkSpace handling methods
// Register scoring function with a workspace
// Base class version just registers itself
void BaseObject::Register(WorkSpace *pWorkSpace) {
  // First unregister with current workspace
  Unregister();
  m_workspace = pWorkSpace;
  LOG_F(1, "BaseObject::Register: Registering {} with workspace", GetName());
  m_workspace->Attach(this);
}

// Unregister with a workspace
// Base class version just unregisters itself
void BaseObject::Unregister() {
  if (m_workspace) {
    LOG_F(1, "BaseObject::Unregister: Unregistering {} from workspace",
          GetName());
    m_workspace->Detach(this);
  }
  m_workspace = nullptr;
}

// Get workspace pointer
WorkSpace *BaseObject::GetWorkSpace() const { return m_workspace; }

// Override Observer method
// Notify observer that subject is about to be deleted
void BaseObject::Deleted(Subject *theDeletedSubject) {
  if (theDeletedSubject == m_workspace) {
    LOG_F(1,
          "BaseObject::Deleted: {} received notice that subject "
          "is to be deleted",
          GetName());
    Unregister();
  }
}

// Request Handling method
void BaseObject::HandleRequest(RequestPtr spRequest) {
  // Base class can handle ENABLE, DISABLE and SETPARAM
  VariantList params = spRequest->GetParameters();
  switch (spRequest->GetID()) {

    // Enable object
    // params[0] is fully-qualified object name
  case ID_REQ_SF_ENABLE:
    Assert<Assertion>(!REQ_CHECK || params.size() == 1);
    if (!params.empty() && (params[0].GetString() == GetFullName())) {
      LOG_F(1, "BaseObject::HandleRequest: Enabling {}", GetFullName());
      Enable();
    }
    break;

    // Disable object
    // params[0] is fully-qualified object name
  case ID_REQ_SF_DISABLE:
    Assert<Assertion>(!REQ_CHECK || params.size() == 1);
    if (!params.empty() && (params[0].GetString() == GetFullName())) {
      LOG_F(1, "BaseObject::HandleRequest: Disabling {}", GetFullName());
      Disable();
    }
    break;

    // Set parameter has two flavours
  case ID_REQ_SF_SETPARAM:
    // Either:
    // params[0] is fully-qualified object name
    // params[1] is parameter name
    // params[2] is parameter value
    // Check if object has this named parameter before calling SetParameter
    if ((params.size() == 3) && (params[0].GetString() == GetFullName()) &&
        isParameterValid(params[1])) {
      LOG_F(1, "BaseObject::HandleRequest: Setting parameter {} to {} for {}",
            params[1].GetString(), params[2].GetString(), GetFullName());
      SetParameter(params[1], params[2]);
    }
    // Or:
    // params[0] is parameter name
    // params[1] is parameter value
    // Check if object has this named parameter before calling SetParameter
    else if ((params.size() == 2) && isParameterValid(params[0])) {
      LOG_F(1, "BaseObject::HandleRequest: Setting parameter {} to {} for {}",
            params[0].GetString(), params[1].GetString(), GetFullName());
      SetParameter(params[0], params[1]);
    }
    break;

    // Ignore all other requests
  default:
    LOG_F(1, "BaseObject::HandleRequest: {} is ignoring request",
          GetFullName());
    break;
  }
}
// Virtual function for dumping parameters to an output stream
// Called by operator <<
void BaseObject::Print(std::ostream &s) const {
  LOG_F(2, "BaseObject::Print: {}", GetFullName());
  ParamHandler::Print(s);
}

// DM 25 Oct 2000 - track changes to parameter values in local data members
// ParameterUpdated is invoked by ParamHandler::SetParameter
void BaseObject::ParameterUpdated(const std::string &strName) {
  if (strName == _ENABLED) {
    m_enabled = GetParameter(_ENABLED);
  } else {
    ParamHandler::ParameterUpdated(strName);
  }
}

void rxdock::to_json(json &j, const BaseObject &baseObject) {
  // workspace pointer skipped
  j = json{{"enabled", baseObject.m_enabled}};
}

void rxdock::from_json(const json &j, BaseObject &baseObject) {
  // model pointer skipped
  j.at("enabled").get_to(baseObject.m_enabled);
}