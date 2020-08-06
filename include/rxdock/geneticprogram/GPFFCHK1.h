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

// GPFFCHK1 Class. GPFFCHK1 represents a fitness function
// valid for the CGP. It is used to train with the CHK1 set.

#ifndef _RBT_GPFFCHK1_H_
#define _RBT_GPFFCHK1_H_

#include "rxdock/Context.h"
#include "rxdock/Rand.h"
#include "rxdock/geneticprogram/GPFitnessFunction.h"
#include "rxdock/geneticprogram/GPGenome.h"
#include "rxdock/geneticprogram/GPTypes.h"
// #include <cmath>

namespace rxdock {

namespace geneticprogram {

class GPFFCHK1 : public GPFitnessFunction {
public:
  static std::string _CT;
  GPFFCHK1(ContextPtr c) : contextp(c) { _RBTOBJECTCOUNTER_CONSTR_(_CT); }
  virtual ~GPFFCHK1() { _RBTOBJECTCOUNTER_DESTR_(_CT); }

  void ReadTables(std::istream &, ReturnTypeArray &, ReturnTypeArray &);
  double CalculateFitness(GPGenomePtr, ReturnTypeArray &, ReturnTypeArray &,
                          double, bool);
  double CalculateFitness(GPGenomePtr, ReturnTypeArray &, ReturnTypeArray &,
                          bool);

private:
  void CreateRandomCtes(int);
  std::vector<double> ctes;
  CellContextPtr contextp;
  GPFFCHK1(); // default constructor disabled
  ReturnTypeArray inputTestTable, inputTrainingTable;
  ReturnTypeArray sfTestTable, sfTrainingTable;
};

// Useful typedefs
typedef SmartPtr<GPFFCHK1> GPFFCHK1Ptr;        // Smart pointer
typedef std::vector<GPFFCHK1Ptr> GPFFCHK1List; // Vector of smart pointers
typedef GPFFCHK1List::iterator GPFFCHK1ListIter;
typedef GPFFCHK1List::const_iterator GPFFCHK1ListConstIter;

} // namespace geneticprogram

} // namespace rxdock

#endif //_GPFFCHK1
