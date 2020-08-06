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

#ifdef _RBTOBJECTCOUNTER
#include "rxdock/Debug.h"

using namespace rxdock;

ObjectCounter::ObjectCounter() : nConstr(0), nCopyConstr(0), nDestr(0) {}

ObjectCounter::~ObjectCounter() {}

ostream &operator<<(std::ostream &s, const ObjectCounter &counter) {
  return s << "Constr=" << counter.nConstr
           << ", Copy constr=" << counter.nCopyConstr
           << ", Destr=" << counter.nDestr << ", Net total="
           << counter.nConstr + counter.nCopyConstr - counter.nDestr;
}

// The one and only ObjectCounterMap definition
ObjectCounterMap rxdock::theObjectCounterMap;

#endif //_RBTOBJECTCOUNTER
