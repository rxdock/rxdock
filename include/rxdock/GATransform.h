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

// Evolves an existing population using a GA
#ifndef _RBTGATRANSFORM_H_
#define _RBTGATRANSFORM_H_

#include "BaseBiMolTransform.h"
#include "Population.h"
#include "Rand.h"

namespace rxdock {

class GATransform : public BaseBiMolTransform {
public:
  static std::string _CT;
  // New individuals to create each cycle, as fraction of population size
  static std::string _NEW_FRACTION;
  // Probability of crossover (1-probability of mutation)
  static std::string _PCROSSOVER;
  // If true, perform Cauchy mutation after each crossover
  static std::string _XOVERMUT;
  // If true, mutations are from Cauchy distribution; if false, from rect.
  // distribution
  static std::string _CMUTATE;
  // Relative step size for mutations (relative to absolute step sizes defined
  // for each chromosome element)
  static std::string _STEP_SIZE;
  // Two genomes are considered equal if the maximum relative difference
  // between chromosome elements is less than _EQUALITY_THRESHOLD
  static std::string _EQUALITY_THRESHOLD;
  // Maximum number of cycles
  static std::string _NCYCLES;
  // Terminate if the best score does not improve over _NCONVERGENCE
  // consecutive cycles
  static std::string _NCONVERGENCE;
  // Output the best pose every _HISTORY_FREQ cycles.
  static std::string _HISTORY_FREQ;

  ////////////////////////////////////////
  // Constructors/destructors
  ////////////////////////////////////////
  RBTDLL_EXPORT GATransform(const std::string &strName = "GAGENRW");
  virtual ~GATransform();

protected:
  ////////////////////////////////////////
  // Protected methods
  ///////////////////
  virtual void SetupReceptor(); // Called by Update when receptor is changed
  virtual void SetupLigand();   // Called by Update when ligand is changed
  virtual void
  SetupTransform(); // Called by Update when either model has changed
  virtual void Execute();

private:
  ////////////////////////////////////////
  // Private methods
  /////////////////
  GATransform(const GATransform &); // Copy constructor disabled by
                                    // default
  GATransform &
  operator=(const GATransform &); // Copy assignment disabled by default

private:
  Rand &m_rand;
};

} // namespace rxdock

#endif //_RBTGATRANSFORM_H_
