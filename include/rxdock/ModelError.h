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

// rxdock model exceptions (e.g. topological errors)

#ifndef _RBTMODELERROR_H_
#define _RBTMODELERROR_H_

#include "rxdock/Error.h"

namespace rxdock {

const std::string IDS_MODEL_ERROR = "RBT_MODEL_ERROR";

// Unspecified model error
class ModelError : public Error {
public:
  ModelError(const std::string &strFile, int nLine,
             const std::string &strMessage = "")
      : Error(IDS_MODEL_ERROR, strFile, nLine, strMessage) {}
  // Protected constructor to allow derived model error classes to set error
  // name
protected:
  ModelError(const std::string &strName, const std::string &strFile, int nLine,
             const std::string &strMessage = "")
      : Error(strName, strFile, nLine, strMessage) {}
};

} // namespace rxdock

#endif //_RBTMODELERROR_H_
