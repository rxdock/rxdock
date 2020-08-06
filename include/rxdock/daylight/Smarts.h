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

// SMARTS class. Defines a function that takes a model and a smarts string as
// input, and returns a list of lists of atoms that match the pattern

#ifndef _RBT_SMARTS_H_
#define _RBT_SMARTS_H_

#include "rxdock/Model.h"

namespace rxdock {

namespace daylight {
AtomListList QueryModel(ModelPtr spModel, const std::string &strSmart,
                        std::string &strSmiles);
}

// Map of (key=integer, value=integer)
typedef std::map<int, int> IntIntMap;
typedef IntIntMap::iterator IntIntMapIter;
typedef IntIntMap::const_iterator IntIntMapConstIter;

} // namespace rxdock

#endif
