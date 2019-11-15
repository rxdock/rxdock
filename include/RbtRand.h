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

#include "RandInt.h"
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
  // Seed the random number generator from the system clock
  void SeedFromClock();
  // Returns current seed
  RBTDLL_EXPORT int GetSeed();
  // Get a random double between 0 and 1 (inlined)
  double GetRandom01() { return m_rand.fdraw(); }
  // Get a random integer between 0 and nMax-1
  int GetRandomInt(int nMax);
  // Get a random unit vector distributed evenly over the surface of a sphere
  RbtVector GetRandomUnitVector();
  double GetGaussianRandom(double, double);
  double GetCauchyRandom(double, double);

private:
  Randint m_rand; // Random number generator
};

///////////////////////////////////////
// Non-member functions in Rbt namespace

namespace Rbt {
// Returns reference to single instance of RbtRand class (singleton)
RBTDLL_EXPORT RbtRand &GetRbtRand();
} // namespace Rbt
#endif //_RBTRAND_H_
