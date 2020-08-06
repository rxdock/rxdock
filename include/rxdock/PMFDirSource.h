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

#include "DirectorySource.h"
#include "PMF.h"

namespace rxdock {

class PMFDirSource : public DirectorySource {
public:
  PMFDirSource(const std::string &);
  virtual ~PMFDirSource() {}

  virtual void ReadFiles(std::vector<std::vector<PMFValue>> *aVect,
                         std::vector<std::string> *aNameVect,
                         std::vector<PMFValue> *aSlopeVect);
  virtual void ParseLines(std::vector<std::string> anStrVect,
                          std::vector<PMFValue> *aValueVect);
};

} // namespace rxdock

#endif // _RBTPMFDIRSOURCE_H_
