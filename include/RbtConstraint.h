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

// Base class for pharmacophore constraints

#ifndef _RBT_CONSTRAINT_H_
#define _RBT_CONSTRAINT_H_

#include "RbtAtom.h"
#include "RbtCoord.h"
#include "RbtError.h"
#include "RbtModel.h"

class RbtConstraint {
public:
  static std::string _CT;
  ///////////////////
  // Constructors
  ///////////////////
  RbtConstraint(RbtCoord, double);
  RbtConstraint(const RbtConstraint &);
  RbtConstraint();

  ///////////////////
  // Destructor
  //////////////////
  virtual ~RbtConstraint();

  void copy(const RbtConstraint &);
  std::ostream &Print(std::ostream &) const;
  friend std::ostream &operator<<(std::ostream &, const RbtConstraint &);

  RbtCoord GetCoords() const { return coord; }
  double GetTolerance() const { return tolerance; }
  virtual void AddAtomList(RbtModelPtr, bool bCheck = true) = 0;
  double Score() const;

  ////////////////////
  // Private methods
  ////////////////////
protected:
  RbtCoord coord;
  double tolerance;
  RbtAtomList m_atomList;

private:
  /////////////////////
  // Private data
  /////////////////////
};

// Useful typedefs
typedef SmartPtr<RbtConstraint> RbtConstraintPtr;   // Smart pointer
typedef vector<RbtConstraintPtr> RbtConstraintList; // Vector of smart pointers
typedef RbtConstraintList::iterator RbtConstraintListIter;
typedef RbtConstraintList::const_iterator RbtConstraintListConstIter;

namespace Rbt {
RbtConstraintPtr CreateConstraint(RbtCoord &c, double &t, std::string &n,
                                  bool bCount = true);
void ZeroCounters();
void ReadConstraint(std::istream &ifile, RbtConstraintPtr &cnt,
                    bool bCount = true);
void ReadConstraintFromMoe(std::istream &ifile, RbtConstraintPtr &cnt,
                           bool bCount = true);
void ReadStartMoe(std::istream &ifile);
void ReadConstraints(std::istream &ii, RbtConstraintList &cl,
                     bool bCount = true);
} // namespace Rbt

// 7 Feb 2005 (DM, Enspiral Discovery)
// New constraint type - any heavy atom
// Likely to be inefficient for large ligands, but has the advantage
// of being simple to code and run
class RbtHeavyConstraint : public RbtConstraint {
public:
  RbtHeavyConstraint(RbtCoord c, double t) : RbtConstraint(c, t) {}
  void AddAtomList(RbtModelPtr, bool bCheck = true);
  static unsigned int counter;
};

class RbtHBAConstraint : public RbtConstraint {
public:
  RbtHBAConstraint(RbtCoord c, double t) : RbtConstraint(c, t) {}
  void AddAtomList(RbtModelPtr, bool bCheck = true);
  static unsigned int counter;
};

class RbtHBDConstraint : public RbtConstraint {
public:
  RbtHBDConstraint(RbtCoord c, double t) : RbtConstraint(c, t) {}
  void AddAtomList(RbtModelPtr, bool bCheck = true);
  static unsigned int counter;
};

class RbtHydroConstraint : public RbtConstraint {
public:
  RbtHydroConstraint(RbtCoord c, double t) : RbtConstraint(c, t) {}
  void AddAtomList(RbtModelPtr, bool bCheck = true);
  static unsigned int counter;
};

class RbtHydroAliphaticConstraint : public RbtConstraint {
public:
  RbtHydroAliphaticConstraint(RbtCoord c, double t) : RbtConstraint(c, t) {}
  void AddAtomList(RbtModelPtr, bool bCheck = true);
  static unsigned int counter;
};

class RbtHydroAromaticConstraint : public RbtConstraint {
public:
  RbtHydroAromaticConstraint(RbtCoord c, double t) : RbtConstraint(c, t) {}
  void AddAtomList(RbtModelPtr, bool bCheck = true);
  static unsigned int counter;
};

class RbtNegChargeConstraint : public RbtConstraint {
public:
  RbtNegChargeConstraint(RbtCoord c, double t) : RbtConstraint(c, t) {}
  void AddAtomList(RbtModelPtr, bool bCheck = true);
  static unsigned int counter;
};

class RbtPosChargeConstraint : public RbtConstraint {
public:
  RbtPosChargeConstraint(RbtCoord c, double t) : RbtConstraint(c, t) {}
  void AddAtomList(RbtModelPtr, bool bCheck = true);
  static unsigned int counter;
};

class RbtRingAromaticConstraint : public RbtConstraint {
public:
  RbtRingAromaticConstraint(RbtCoord c, double t) : RbtConstraint(c, t) {}
  void AddAtomList(RbtModelPtr, bool bCheck = true);
  static unsigned int counter;
};

#endif //_RbtConstraint_H_
