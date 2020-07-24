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

// Scoring function specific requests

#ifndef _RBTSFREQUEST_H_
#define _RBTSFREQUEST_H_

#include "Request.h"

namespace rxdock {

// Request IDs
const RequestID ID_REQ_SF_ENABLE = 1;
const RequestID ID_REQ_SF_DISABLE = 2;
const RequestID ID_REQ_SF_PARTITION = 3;
const RequestID ID_REQ_SF_SETPARAM = 4;

class SFEnableRequest : public Request {
public:
  // Request to enable a particular named scoring function
  SFEnableRequest(const std::string &sfName) : Request(ID_REQ_SF_ENABLE) {
    AddParameter(sfName);
  }
};

class SFDisableRequest : public Request {
public:
  // Request to disable a particular named scoring function
  SFDisableRequest(const std::string &sfName) : Request(ID_REQ_SF_DISABLE) {
    AddParameter(sfName);
  }
};

class SFPartitionRequest : public Request {
public:
  // Request to partition all scoring functions to a particular distance
  SFPartitionRequest(double dist) : Request(ID_REQ_SF_PARTITION) {
    AddParameter(dist);
  }
  // Request to partition a particular named scoring function to a particular
  // distance
  SFPartitionRequest(const std::string &sfName, double dist)
      : Request(ID_REQ_SF_PARTITION) {
    AddParameter(sfName);
    AddParameter(dist);
  }
};

class SFSetParamRequest : public Request {
public:
  // Request to set a named parameter of all scoring functions to a new value
  // (e.g. TRACE)
  SFSetParamRequest(const std::string &paramName, const Variant &paramValue)
      : Request(ID_REQ_SF_SETPARAM) {
    AddParameter(paramName);
    AddParameter(paramValue);
  }
  // Request to set a named parameter of a named scoring function to a new value
  SFSetParamRequest(const std::string &sfName, const std::string &paramName,
                    const Variant &paramValue)
      : Request(ID_REQ_SF_SETPARAM) {
    AddParameter(sfName);
    AddParameter(paramName);
    AddParameter(paramValue);
  }
};

} // namespace rxdock

#endif //_RBTSFREQUEST_H_
