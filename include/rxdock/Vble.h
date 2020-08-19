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

// Vble class. Stores information about the variables used by the filters.

#ifndef _RBT_VBLE_H_
#define _RBT_VBLE_H_

#include "rxdock/geneticprogram/GPTypes.h"

namespace rxdock {

class Vble {
public:
  static std::string _CT;
  enum VbleType { CTE, LIG, SCORE, SITE };
  Vble() : vt(CTE), value(0.0), name("") {
    _RBTOBJECTCOUNTER_CONSTR_(_CT);
  } // default const.
  Vble(std::string s, ReturnType val = 0.0) : value(val), name(s) {
    if (s.find("LIG_") == 0)
      vt = LIG;
    else if (s.find("SITE_") == 0)
      vt = SITE;
    else if (s.find(GetMetaDataPrefix() + "score") == 0)
      vt = SCORE;
    else
      vt = CTE;
    _RBTOBJECTCOUNTER_CONSTR_(_CT);
  }
  Vble(const Vble &v) : vt(v.vt), value(v.value), name(v.name) {
    _RBTOBJECTCOUNTER_COPYCONSTR_(_CT);
  }
  void SetValue(ReturnType val) {
    value = val;
    if (name == "") {
      name = std::to_string(value);
    }
  }
  ReturnType GetValue() const { return value; }
  virtual ~Vble() { _RBTOBJECTCOUNTER_DESTR_(_CT); }

  Vble &operator=(const Vble &v) {
    vt = v.vt;
    value = v.value;
    name = v.name;
    return *this;
  }

  void SetName(std::string nm) { name = nm; }
  std::string GetName() const { return name; }

  bool IsLig() { return (vt == LIG); }
  bool IsScore() { return (vt == SCORE); }
  bool IsSite() { return (vt == SITE); }

private:
  VbleType vt;
  ReturnType value;
  std::string name;
};

// Useful typedefs
typedef SmartPtr<Vble> VblePtr;        // Smart pointer
typedef std::vector<VblePtr> VbleList; // Vector of smart pointers
typedef VbleList::iterator VbleListIter;
typedef VbleList::const_iterator VbleListConstIter;

} // namespace rxdock

#endif //_Vble_H_
