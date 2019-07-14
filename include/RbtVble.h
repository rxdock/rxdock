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

#include "RbtGPTypes.h"

class RbtVble {
public:
  static std::string _CT;
  enum VbleType { CTE, LIG, SCORE, SITE };
  RbtVble() : vt(CTE), value(0.0), name("") {
    _RBTOBJECTCOUNTER_CONSTR_(_CT);
  } // default const.
  RbtVble(std::string s, RbtReturnType val = 0.0) : value(val), name(s) {
    if (s.find("LIG_") == 0)
      vt = LIG;
    else if (s.find("SITE_") == 0)
      vt = SITE;
    else if (s.find("SCORE") == 0)
      vt = SCORE;
    else
      vt = CTE;
    _RBTOBJECTCOUNTER_CONSTR_(_CT);
  }
  RbtVble(const RbtVble &v) : vt(v.vt), value(v.value), name(v.name) {
    _RBTOBJECTCOUNTER_COPYCONSTR_(_CT);
  }
  void SetValue(RbtReturnType val) {
    value = val;
    if (name == "") {
      name = std::to_string(value);
    }
  }
  RbtReturnType GetValue() const { return value; }
  virtual ~RbtVble() { _RBTOBJECTCOUNTER_DESTR_(_CT); }
  void SetName(std::string nm) { name = nm; }
  std::string GetName() const { return name; }

  bool IsLig() { return (vt == LIG); }
  bool IsScore() { return (vt == SCORE); }
  bool IsSite() { return (vt == SITE); }

private:
  VbleType vt;
  RbtReturnType value;
  std::string name;
};

// Useful typedefs
typedef SmartPtr<RbtVble> RbtVblePtr;        // Smart pointer
typedef std::vector<RbtVblePtr> RbtVbleList; // Vector of smart pointers
typedef RbtVbleList::iterator RbtVbleListIter;
typedef RbtVbleList::const_iterator RbtVbleListConstIter;

#endif //_RbtVble_H_
