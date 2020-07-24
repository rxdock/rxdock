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

#include "GPFitnessFunction.h"
#include "Debug.h"
#include "GPGenome.h"
#include "GPParser.h"
#include <fstream>
#include <sstream>

using namespace rxdock;
using namespace rxdock::geneticprogram;

std::string GPFitnessFunction::_CT("GPFitnessFunction");

// Constructors
GPFitnessFunction::GPFitnessFunction()
    : m_rand(GetRandInstance()), objective(0.0), fitness(0.0) {
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

GPFitnessFunction::GPFitnessFunction(const GPFitnessFunction &g)
    : m_rand(GetRandInstance()), objective(g.objective), fitness(g.fitness) {
  _RBTOBJECTCOUNTER_COPYCONSTR_(_CT);
}

GPFitnessFunction &GPFitnessFunction::operator=(const GPFitnessFunction &g) {
  if (this != &g) {
    objective = g.objective;
    fitness = g.fitness;
  }
  return *this;
}

///////////////////
// Destructor
//////////////////
GPFitnessFunction::~GPFitnessFunction() { _RBTOBJECTCOUNTER_DESTR_(_CT); }

double GPFitnessFunction::GetObjective() const { return objective; }

double GPFitnessFunction::GetFitness() const { return fitness; }

void GPFitnessFunction::SetFitness(double f) { fitness = f; }
