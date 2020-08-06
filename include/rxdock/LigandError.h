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

// rxdock ligand exceptions, the model is correct, but the ligand is not
// appropiate for this case. (e.g. doesn't comply with all the pharmacophore
// constraints)

#ifndef _RBTLIGANDERROR_H_
#define _RBTLIGANDERROR_H_

#include "rxdock/Error.h"

namespace rxdock {

const std::string IDS_LIGAND_ERROR = "RBT_LIGAND_ERROR";

// Unspecified ligand error
class LigandError : public Error {
public:
  LigandError(const std::string &strFile, int nLine,
              const std::string &strMessage = "")
      : Error(IDS_LIGAND_ERROR, strFile, nLine, strMessage) {}
  // Protected constructor to allow derived ligand error classes to set error
  // name
protected:
  LigandError(const std::string &strName, const std::string &strFile, int nLine,
              const std::string &strMessage = "")
      : Error(strName, strFile, nLine, strMessage) {}
};

} // namespace rxdock

#endif //_RBTLIGANDERROR_H_
