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

// GPFFSpike Class. GPFFSpike represents a fitness function
// valid for the CGP. It is used to train with the Spike set.

#ifndef _RBT_GPFFSPIKE_H_
#define _RBT_GPFFSPIKE_H_

#include "Context.h"
#include "GPFitnessFunction.h"
#include "GPGenome.h"
#include "GPTypes.h"
#include "Rand.h"
// #include <cmath>

namespace rxdock {

namespace geneticprogram {

class GPFFSpike : public GPFitnessFunction {
public:
  static std::string _CT;
  GPFFSpike(ContextPtr c) : contextp(c) { _RBTOBJECTCOUNTER_CONSTR_(_CT); }
  virtual ~GPFFSpike() { _RBTOBJECTCOUNTER_DESTR_(_CT); }

  void ReadTables(std::istream &, ReturnTypeArray &, ReturnTypeArray &);
  double CalculateFitness(GPGenomePtr, ReturnTypeArray &, ReturnTypeArray &,
                          double, bool);
  double CalculateFitness(GPGenomePtr, ReturnTypeArray &, ReturnTypeArray &,
                          bool);

private:
  void CreateRandomCtes(int);
  std::vector<double> ctes;
  CellContextPtr contextp;
  GPFFSpike(); // default constructor disabled
};

// Useful typedefs
typedef SmartPtr<GPFFSpike> GPFFSpikePtr;        // Smart pointer
typedef std::vector<GPFFSpikePtr> GPFFSpikeList; // Vector of smart pointers
typedef GPFFSpikeList::iterator GPFFSpikeListIter;
typedef GPFFSpikeList::const_iterator GPFFSpikeListConstIter;

} // namespace geneticprogram

} // namespace rxdock

#endif //_GPFFSpike
