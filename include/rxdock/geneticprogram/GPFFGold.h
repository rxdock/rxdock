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

// GPFFGold Class. GPFFGold represents a fitness function
// valid for the CGP. It is used to train with the Gold set.

#ifndef _RBT_GPFFGOLD_H_
#define _RBT_GPFFGOLD_H_

#include "GPFitnessFunction.h"
#include "GPGenome.h"
#include "GPTypes.h"
#include "Rand.h"
// #include <cmath>

namespace rxdock {

namespace geneticprogram {

class GPFFGold : public GPFitnessFunction {
public:
  static std::string _CT;

  void ReadTables(std::istream &, ReturnTypeArray &, ReturnTypeArray &);
  double CalculateFitness(GPGenomePtr, ReturnTypeArray &, ReturnTypeArray &,
                          double, bool);
  double CalculateFitness(GPGenomePtr, ReturnTypeArray &, ReturnTypeArray &,
                          bool);
};

// Useful typedefs
typedef SmartPtr<GPFFGold> GPFFGoldPtr;        // Smart pointer
typedef std::vector<GPFFGoldPtr> GPFFGoldList; // Vector of smart pointers
typedef GPFFGoldList::iterator GPFFGoldListIter;
typedef GPFFGoldList::const_iterator GPFFGoldListConstIter;

} // namespace geneticprogram

} // namespace rxdock

#endif //_GPFFGold
