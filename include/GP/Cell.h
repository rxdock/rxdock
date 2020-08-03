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

// Cell Class. Cell represents a function or input variable of the CGP

#ifndef _RBT_CELL_H_
#define _RBT_CELL_H_

#include "GPTypes.h"

#include <loguru.hpp>

namespace rxdock {

class Cell {
public:
  ///////////////////
  // Constructors
  ///////////////////
  Cell() {
    evaluated = false;
    name = "";
    constant = false;
  }
  Cell(const Cell &c)
      : evaluated(c.evaluated), constant(c.constant), name(c.name),
        result(c.result) {}

  ///////////////////
  // Destructor
  //////////////////
  virtual ~Cell() {}

  Cell &operator=(const Cell &c) {
    evaluated = c.evaluated;
    constant = c.constant;
    name = c.name;
    result = c.result;
    return *this;
  }

  bool Evaluated() const { return evaluated; }
  bool Named() const { return (name != ""); }
  bool Constant() const { return constant; }
  void Clear() {
    if (!constant) {
      evaluated = false;
      name = "";
    }
  }
  ReturnType GetResult() const { return result; }
  std::string GetName() const { return name; }
  void SetName(std::string s) { name = s; }
  void SetResult(ReturnType r) {
    result = r;
    evaluated = true;
  }
  void SetConstant(ReturnType r) {
    result = r;
    evaluated = true;
    constant = true;
    name = std::to_string(r);
  }
  void ResetConstant() {
    evaluated = false;
    name = "";
    if (constant)
      constant = false;
    else
      LOG_F(WARNING, "Cell::ResetConstant: no constant");
  }

private:
  bool evaluated, constant;
  std::string name;
  ReturnType result;
  //    friend std::ostream& operator<<(std::ostream& s, const Cell &p);
};

// Useful typedefs
typedef SmartPtr<Cell> CellPtr;        // Smart pointer
typedef std::vector<CellPtr> CellList; // Vector of smart pointers
typedef CellList::iterator CellListIter;
typedef CellList::const_iterator CellListConstIter;

} // namespace rxdock

#endif //_Cell_H_
