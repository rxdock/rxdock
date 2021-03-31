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

// GPFFHSP90 Class. GPFFHSP90 represents a fitness function
// valid for the CGP. It is used to train with the HSP90 set.

#ifndef _RBT_GPFFHSP90_H_
#define _RBT_GPFFHSP90_H_

#include "rxdock/Context.h"
#include "rxdock/Rand.h"
#include "rxdock/geneticprogram/GPFitnessFunction.h"
#include "rxdock/geneticprogram/GPGenome.h"
#include "rxdock/geneticprogram/GPTypes.h"

namespace rxdock {

namespace geneticprogram {

class GPFFHSP90 : public GPFitnessFunction {
public:
  static const std::string _CT;
  GPFFHSP90(ContextPtr c) : contextp(c) { _RBTOBJECTCOUNTER_CONSTR_(_CT); }
  virtual ~GPFFHSP90() { _RBTOBJECTCOUNTER_DESTR_(_CT); }

  void ReadTables(std::istream &, ReturnTypeArray &, ReturnTypeArray &);
  double CalculateFitness(GPGenomePtr, ReturnTypeArray &, ReturnTypeArray &,
                          double, bool);
  double CalculateFitness(GPGenomePtr, ReturnTypeArray &, ReturnTypeArray &,
                          bool);

private:
  void CreateRandomCtes(int);
  std::vector<double> ctes;
  CellContextPtr contextp;
  GPFFHSP90(); // default constructor disabled
  ReturnTypeArray inputTestTable, inputTrainingTable;
  ReturnTypeArray sfTestTable, sfTrainingTable;
};

// Useful typedefs
typedef SmartPtr<GPFFHSP90> GPFFHSP90Ptr;        // Smart pointer
typedef std::vector<GPFFHSP90Ptr> GPFFHSP90List; // Vector of smart pointers
typedef GPFFHSP90List::iterator GPFFHSP90ListIter;
typedef GPFFHSP90List::const_iterator GPFFHSP90ListConstIter;

} // namespace geneticprogram

} // namespace rxdock

#endif //_GPFFHSP90
