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

#include "rxdock/Rand.h"
#include "rxdock/Singleton.h"

using namespace rxdock;

/////////////
// Constructor
Rand::Rand() {
  // Seed the random number generator
  // Fixed seed in debug mode
  // Seed from the random device in release mode
#ifdef _DEBUG
  Seed();
#else  //_DEBUG
  SeedFromRandomDevice();
#endif //_DEBUG
  _RBTOBJECTCOUNTER_CONSTR_("Rand");
}

/////////////
// Destructor
Rand::~Rand() { _RBTOBJECTCOUNTER_DESTR_("Rand"); }

////////////////
// Public methods

// Seed the random number generator
void Rand::Seed(int seed) { m_rng.seed(seed); }

// Seed the random number generator from the random device
void Rand::SeedFromRandomDevice() {
#if defined(__sun) || (defined(_WIN32) && defined(_MSC_VER))
  std::random_device rd;
  m_rng.seed(rd());
#else
  pcg_extras::seed_seq_from<std::random_device> seedSource;
  m_rng.seed(seedSource);
#endif
}

// Get a random double between 0 and 1
double Rand::GetRandom01() {
  std::uniform_real_distribution<> uniformDist(0.0, 1.0);
  return uniformDist(m_rng);
}

// Get a random integer between 0 and nMax-1
int Rand::GetRandomInt(int nMax) {
  if (nMax == 0) {
    return nMax;
  }
  std::uniform_int_distribution<> uniformDist(0, nMax - 1);
  return uniformDist(m_rng);
}

// Get a random unit vector distributed evenly over the surface of a sphere
Vector Rand::GetRandomUnitVector() {
  double z = 2.0 * GetRandom01() - 1.0;
  double t = 2.0 * M_PI * GetRandom01();
  double w = std::sqrt(1.0 - z * z);
  double x = w * std::cos(t);
  double y = w * std::sin(t);
  return Vector(x, y, z);
}

// Get a random number from the Normal distribution (mean, variance)
double Rand::GetGaussianRandom(double mean, double variance) {
  std::normal_distribution<> gaussianDist{mean, variance};
  return gaussianDist(m_rng);
}

// Get a random number from the Cauchy distribution (mean, variance)
double Rand::GetCauchyRandom(double mean, double variance) {
  std::cauchy_distribution<> cauchyDist{mean, variance};
  return cauchyDist(m_rng);
}

///////////////////////////////////////
// Non-member functions in rxdock namespace

// Returns reference to single instance of Rand class (singleton)
Rand &rxdock::GetRandInstance() { return Singleton<Rand>::instance(); }
