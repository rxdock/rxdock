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

#ifndef _RBTPMFDIRSOURCE_H_
#define _RBTPMFDIRSOURCE_H_

#include "RbtDirectorySource.h"
#include "RbtPMF.h"

namespace rxdock {

class RbtPMFDirSource : public RbtDirectorySource {
public:
  RbtPMFDirSource(const std::string &);
  virtual ~RbtPMFDirSource() {}

  virtual void ReadFiles(std::vector<std::vector<RbtPMFValue>> *aVect,
                         std::vector<std::string> *aNameVect,
                         std::vector<RbtPMFValue> *aSlopeVect);
  virtual void ParseLines(std::vector<std::string> anStrVect,
                          std::vector<RbtPMFValue> *aValueVect);
};

} // namespace rxdock

#endif // _RBTPMFDIRSOURCE_H_
