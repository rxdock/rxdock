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

// Wrapper around Randint class
// Function provided to return reference to single instance (singleton) of
// RbtRand

#ifndef _RBTRAND_H_
#define _RBTRAND_H_

// Use standard C++11 RNG on Solaris and Windows/MSVC due to PCG build failure
// Solaris issue: https://github.com/imneme/pcg-cpp/issues/42
// Windows/MSVC issue: https://github.com/imneme/pcg-cpp/issues/11
#if !defined(__sun) && !(defined(_WIN32) && defined(_MSC_VER))
#include <pcg_random.hpp>
#endif
#include <random>

#include "RbtCoord.h"

class RbtRand {
  /////////////
  // Constructor
public:
  RbtRand();
  /////////////
  // Destructor
  ~RbtRand();

  ////////////////
  // Public methods

  // Seed the random number generator
  RBTDLL_EXPORT void Seed(int seed = 0);
  // Seed the random number generator from the random device
  void SeedFromRandomDevice();
  // Returns current seed
  RBTDLL_EXPORT int GetSeed();
  // Get a random double between 0 and 1
  RBTDLL_EXPORT double GetRandom01();
  // Get a random integer between 0 and nMax-1
  int GetRandomInt(int nMax);
  // Get a random unit vector distributed evenly over the surface of a sphere
  RbtVector GetRandomUnitVector();
  double GetGaussianRandom(double, double);
  double GetCauchyRandom(double, double);

private:
#if defined(__sun) || (defined(_WIN32) && defined(_MSC_VER))
  std::default_random_engine m_rng;
#else
  pcg32 m_rng; // Random number generator
#endif
};

///////////////////////////////////////
// Non-member functions in Rbt namespace

namespace Rbt {
// Returns reference to single instance of RbtRand class (singleton)
RBTDLL_EXPORT RbtRand &GetRbtRand();
} // namespace Rbt
#endif //_RBTRAND_H_
