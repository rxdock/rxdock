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

#include "rxdock/BaseTransform.h"
#include "rxdock/BaseSF.h"
#include "rxdock/WorkSpace.h"

#include <loguru.hpp>

using namespace rxdock;

// Static data members
const std::string BaseTransform::_CT = "BaseTransform";

////////////////////////////////////////
// Constructors/destructors
BaseTransform::BaseTransform(const std::string &strClass,
                             const std::string &strName)
    : BaseObject(strClass, strName), m_parent(nullptr) {
  LOG_F(2, "BaseTransform parameterised constructor");
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

BaseTransform::~BaseTransform() {
  LOG_F(2, "BaseTransform destructor");
  Orphan(); // Remove object from parent aggregate
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

////////////////////////////////////////
// Public methods
////////////////

// Fully qualified name, prefixed by all ancestors
std::string BaseTransform::GetFullName() const {
  // If we have a parent, prefix our short name with our parents full name,
  // else full name is just our short name
  return (m_parent) ? m_parent->GetFullName() + "." + GetName() : GetName();
}

// Actually apply the transform
void BaseTransform::Go() {
  if (isEnabled()) {
    // Send any stored Scoring Function requests (e.g. to change any params)
    SendSFRequests();
    Execute();
  }
}

// Aggregate handling (virtual) methods
// Base class throws an InvalidRequest error

void BaseTransform::Add(BaseTransform *) {
  throw InvalidRequest(_WHERE_, "Add() invalid for non-aggregate transforms");
}
void BaseTransform::Remove(BaseTransform *) {
  throw InvalidRequest(_WHERE_,
                       "Remove() invalid for non-aggregate transforms");
}

BaseTransform *BaseTransform::GetTransform(unsigned int i) const {
  throw InvalidRequest(_WHERE_,
                       "GetTransform() invalid for non-aggregate transforms");
}
bool BaseTransform::isAgg() const { return false; }
unsigned int BaseTransform::GetNumTransforms() const { return 0; }

// Aggregate handling (concrete) methods
BaseTransform *BaseTransform::GetParentTransform() const { return m_parent; }

// Force removal from the parent aggregate
void BaseTransform::Orphan() {
  if (m_parent) {
    LOG_F(1, "BaseTransform::Orphan: Removing {} from {}", GetName(),
          m_parent->GetName());
    m_parent->Remove(this);
  }
}

// Scoring function request handling
// Transforms can store up a list of requests to send to the workspace
// scoring function each time Go() is executed
void BaseTransform::AddSFRequest(RequestPtr spRequest) {
  m_SFRequests.push_back(spRequest);
}
void BaseTransform::ClearSFRequests() { m_SFRequests.clear(); }
void BaseTransform::SendSFRequests() {
  // Get the current scoring function from the workspace
  WorkSpace *pWorkSpace = GetWorkSpace();
  if (pWorkSpace == nullptr) // Return if this transform is not registered
    return;
  BaseSF *pSF = pWorkSpace->GetSF();
  if (pSF == nullptr) // Return if workspace does not have a scoring function
    return;
  // Send each request to the scoring function
  for (RequestListConstIter iter = m_SFRequests.begin();
       iter != m_SFRequests.end(); iter++) {
    pSF->HandleRequest(*iter);
  }
}

void rxdock::to_json(json &j, const BaseTransform &bt) {
  // skipping parent pointer
  json requestList;
  for (const auto &aIter : bt.m_SFRequests) {
    json request = *aIter;
    requestList.push_back(request);
  }

  j = json{{"request-list", requestList}};
}

void rxdock::from_json(const json &j, BaseTransform &bt) {
  // skipping parent pointer
  for (auto &request : j.at("request-list")) {
    RequestPtr spRequest = RequestPtr(new Request(request));
    bt.m_SFRequests.push_back(spRequest);
  }
}