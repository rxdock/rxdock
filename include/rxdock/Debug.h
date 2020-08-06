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

// Tools for tracking of object construction and destruction for
// debugging purposes

#ifndef _RBTDEBUG_H_
#define _RBTDEBUG_H_

#ifdef _RBTOBJECTCOUNTER

#include <map>

namespace rxdock {

// Object counter class for debugging no. of object constructions and
// destructions
class ObjectCounter {
public:
  ObjectCounter();
  ~ObjectCounter();
  friend std::ostream &operator<<(std::ostream &s,
                                  const ObjectCounter &counter);
  Int nConstr;
  Int nCopyConstr;
  Int nDestr;
};

// Useful typedefs
typedef std::map<String, ObjectCounter> ObjectCounterMap;
typedef ObjectCounterMap::iterator ObjectCounterMapIter;

// External declaration of the one and only object counter map
extern ObjectCounterMap theObjectCounterMap;

} // namespace rxdock

#define _RBTOBJECTCOUNTER_CONSTR_(class) theObjectCounterMap[(class)].nConstr++
#define _RBTOBJECTCOUNTER_COPYCONSTR_(class)                                   \
  theObjectCounterMap[(class)].nCopyConstr++
#define _RBTOBJECTCOUNTER_DESTR_(class) theObjectCounterMap[(class)].nDestr++

#define _RBTOBJECTCOUNTER_DUMP_(stream)                                        \
  (stream) << std::endl << "RBT OBJECT COUNTERS:" << std::endl << std::endl;   \
  for (ObjectCounterMapIter dumpIter = theObjectCounterMap.begin();            \
       dumpIter != theObjectCounterMap.end(); dumpIter++) {                    \
    (stream) << (*dumpIter).first << ": " << (*dumpIter).second << std::endl;  \
  }

#else //_RBTOBJECTCOUNTER

#define _RBTOBJECTCOUNTER_CONSTR_(class)
#define _RBTOBJECTCOUNTER_COPYCONSTR_(class)
#define _RBTOBJECTCOUNTER_DESTR_(class)
#define _RBTOBJECTCOUNTER_DUMP_(stream)

#endif //_RBTOBJECTCOUNTER

#endif //_RBTDEBUG_H_
