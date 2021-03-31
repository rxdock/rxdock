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

#include "rxdock/Atom.h"
#include "rxdock/Coord.h"
#include "rxdock/Error.h"
#include "rxdock/Model.h"

namespace rxdock {

class Constraint {
public:
  static const std::string _CT;
  ///////////////////
  // Constructors
  ///////////////////
  Constraint(Coord, double);
  Constraint(const Constraint &);
  Constraint();

  ///////////////////
  // Destructor
  //////////////////
  virtual ~Constraint();

  void copy(const Constraint &);
  std::ostream &Print(std::ostream &) const;
  friend std::ostream &operator<<(std::ostream &, const Constraint &);

  Coord GetCoords() const { return coord; }
  double GetTolerance() const { return tolerance; }
  virtual void AddAtomList(ModelPtr, bool bCheck = true) = 0;
  double Score() const;

  ////////////////////
  // Private methods
  ////////////////////
protected:
  Coord coord;
  double tolerance;
  AtomList m_atomList;

private:
  /////////////////////
  // Private data
  /////////////////////
};

// Useful typedefs
typedef SmartPtr<Constraint> ConstraintPtr;        // Smart pointer
typedef std::vector<ConstraintPtr> ConstraintList; // Vector of smart pointers
typedef ConstraintList::iterator ConstraintListIter;
typedef ConstraintList::const_iterator ConstraintListConstIter;

ConstraintPtr CreateConstraint(Coord &c, double &t, std::string &n,
                               bool bCount = true);
void ZeroCounters();
void ReadConstraint(std::istream &ifile, ConstraintPtr &cnt,
                    bool bCount = true);
void ReadConstraintFromMoe(std::istream &ifile, ConstraintPtr &cnt,
                           bool bCount = true);
void ReadStartMoe(std::istream &ifile);
void ReadConstraints(std::istream &ii, ConstraintList &cl, bool bCount = true);

// 7 Feb 2005 (DM, Enspiral Discovery)
// New constraint type - any heavy atom
// Likely to be inefficient for large ligands, but has the advantage
// of being simple to code and run
class HeavyConstraint : public Constraint {
public:
  HeavyConstraint(Coord c, double t) : Constraint(c, t) {}
  void AddAtomList(ModelPtr, bool bCheck = true);
  static unsigned int counter;
};

class HBAConstraint : public Constraint {
public:
  HBAConstraint(Coord c, double t) : Constraint(c, t) {}
  void AddAtomList(ModelPtr, bool bCheck = true);
  static unsigned int counter;
};

class HBDConstraint : public Constraint {
public:
  HBDConstraint(Coord c, double t) : Constraint(c, t) {}
  void AddAtomList(ModelPtr, bool bCheck = true);
  static unsigned int counter;
};

class HydroConstraint : public Constraint {
public:
  HydroConstraint(Coord c, double t) : Constraint(c, t) {}
  void AddAtomList(ModelPtr, bool bCheck = true);
  static unsigned int counter;
};

class HydroAliphaticConstraint : public Constraint {
public:
  HydroAliphaticConstraint(Coord c, double t) : Constraint(c, t) {}
  void AddAtomList(ModelPtr, bool bCheck = true);
  static unsigned int counter;
};

class HydroAromaticConstraint : public Constraint {
public:
  HydroAromaticConstraint(Coord c, double t) : Constraint(c, t) {}
  void AddAtomList(ModelPtr, bool bCheck = true);
  static unsigned int counter;
};

class NegChargeConstraint : public Constraint {
public:
  NegChargeConstraint(Coord c, double t) : Constraint(c, t) {}
  void AddAtomList(ModelPtr, bool bCheck = true);
  static unsigned int counter;
};

class PosChargeConstraint : public Constraint {
public:
  PosChargeConstraint(Coord c, double t) : Constraint(c, t) {}
  void AddAtomList(ModelPtr, bool bCheck = true);
  static unsigned int counter;
};

class RingAromaticConstraint : public Constraint {
public:
  RingAromaticConstraint(Coord c, double t) : Constraint(c, t) {}
  void AddAtomList(ModelPtr, bool bCheck = true);
  static unsigned int counter;
};

} // namespace rxdock

#endif //_Constraint_H_
