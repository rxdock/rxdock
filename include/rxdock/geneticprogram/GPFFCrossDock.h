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

// GPFFCrossDock Class. GPFFCrossDock represents a fitness function
// valid for the CGP. It is used to train with the CrossDock set.

#ifndef _RBT_GPFFCROSSDOCK_H_
#define _RBT_GPFFCROSSDOCK_H_

#include "GPFitnessFunction.h"
#include "GPGenome.h"
#include "GPTypes.h"
#include "Rand.h"
// #include <cmath>

namespace rxdock {

namespace geneticprogram {

class GPFFCrossDock : public GPFitnessFunction {
public:
  static std::string _CT;

  void ReadTables(std::istream &, ReturnTypeArray &, ReturnTypeArray &);
  double CalculateFitness(GPGenomePtr, ReturnTypeArray &, ReturnTypeArray &,
                          double, bool);
  double CalculateFitness(GPGenomePtr, ReturnTypeArray &, ReturnTypeArray &,
                          bool);

private:
  void CreateRandomCtes(int);
  std::vector<double> ctes;
};

// Useful typedefs
typedef SmartPtr<GPFFCrossDock> GPFFCrossDockPtr; // Smart pointer
typedef std::vector<GPFFCrossDockPtr>
    GPFFCrossDockList; // Vector of smart pointers
typedef GPFFCrossDockList::iterator GPFFCrossDockListIter;
typedef GPFFCrossDockList::const_iterator GPFFCrossDockListConstIter;

} // namespace geneticprogram

} // namespace rxdock

#endif //_GPFFCrossDock
