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

// GPFitnessFunction Class. GPFitnessFunction represents a fitness function
// valid for the CGP

#ifndef _RBT_GPFITNESSFUNCTION_H_
#define _RBT_GPFITNESSFUNCTION_H_

#include "rxdock/Context.h"
#include "rxdock/Rand.h"
#include "rxdock/geneticprogram/GPGenome.h"
#include "rxdock/geneticprogram/GPTypes.h"
// #include <cmath>

namespace rxdock {

namespace geneticprogram {

class GPFitnessFunction {
public:
  static const std::string _CT;
  ///////////////////
  // Constructors
  ///////////////////
  GPFitnessFunction();
  GPFitnessFunction(const GPFitnessFunction &);

  GPFitnessFunction &operator=(const GPFitnessFunction &);

  ///////////////////
  // Destructor
  //////////////////
  virtual ~GPFitnessFunction();
  virtual void ReadTables(std::istream &, ReturnTypeArray &,
                          ReturnTypeArray &) = 0;
  double GetObjective() const;
  double GetFitness() const;
  void SetFitness(double);
  virtual double CalculateFitness(GPGenomePtr, ReturnTypeArray &,
                                  ReturnTypeArray &, double, bool) = 0;
  virtual double CalculateFitness(GPGenomePtr, ReturnTypeArray &,
                                  ReturnTypeArray &, bool) = 0;

protected:
  Rand &m_rand;
  double objective; // raw score
  double fitness;   // scaled score
  ReturnTypeArray inputTable, SFTable;

private:
};

// Useful typedefs
typedef SmartPtr<GPFitnessFunction> GPFitnessFunctionPtr; // Smart pointer
typedef std::vector<GPFitnessFunctionPtr>
    GPFitnessFunctionList; // Vector of smart pointers
typedef GPFitnessFunctionList::iterator GPFitnessFunctionListIter;
typedef GPFitnessFunctionList::const_iterator GPFitnessFunctionListConstIter;

} // namespace geneticprogram

} // namespace rxdock

#endif //_GPFitnessFunction
