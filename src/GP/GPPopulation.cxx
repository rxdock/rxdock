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

#include "GPPopulation.h"
#include "Debug.h"

#include <loguru.hpp>

#include <algorithm>
#include <cfloat>

using namespace rxdock;
using namespace rxdock::geneticprogram;

// Static data members for GPPopulation
std::string GPPopulation::_CT("GPPopulation");

// Constructor

GPPopulation::GPPopulation(int s, int nr, GPFitnessFunctionPtr f,
                           ReturnTypeArray &it, ReturnTypeArray &sft)
    : m_rand(GetRandInstance()) {
  popsize = s;
  nrepl = nr;
  ff = f;
  ittrain = it;
  sfttrain = sft;
  LOG_F(1, "const: {}", ittrain.size());
  c = 2.0; // default value for the sigma truncation multiplier
  pop = GPGenomeList(popsize);
  for (GPGenomeListIter iter = pop.begin(); iter != pop.end(); iter++)
    *iter = GPGenomePtr(new GPGenome());
  newpop = GPGenomeList(nrepl);
  for (GPGenomeListIter iter = newpop.begin(); iter != newpop.end(); iter++)
    *iter = GPGenomePtr(new GPGenome());
  psum = new double[popsize];
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

// Copy constructor
// This is copying the pointers. Do I need to change it to
// copying the actual genomes?
GPPopulation::GPPopulation(const GPPopulation &p) : m_rand(GetRandInstance()) {
  popsize = p.popsize;
  nrepl = p.nrepl;
  ff = p.ff;
  ittrain = p.ittrain;
  sfttrain = p.sfttrain;
  c = 2.0; // default value for the sigma truncation multiplier
  pop = GPGenomeList(popsize);
  for (int i = 0; i < p.popsize; i++)
    pop[i] = p.pop[i];
  newpop = GPGenomeList(p.newpop.size());
  for (unsigned int i = 0; i < p.newpop.size(); i++)
    newpop[i] = p.newpop[i];
  psum = new double[popsize];
  _RBTOBJECTCOUNTER_COPYCONSTR_(_CT);
}

// Destructor

GPPopulation::~GPPopulation() {
  delete[] psum;
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

int GPPopulation::GetSize() { return popsize; }

int GPPopulation::GetNrepl() { return nrepl; }

// Initialises the population by calling the initializer for each genome
void GPPopulation::Initialise(double hitlimit, bool function) {
  double bestFitness = -FLT_MAX, fit;
  for (int i = 0; i < popsize; i++) {
    pop[i]->Initialise();
    fit = ff->CalculateFitness(pop[i], ittrain, sfttrain, function);
    LOG_F(1, "init: {}", ittrain.size());
    if (fit > bestFitness) {
      bestFitness = fit;
      bestInd = i;
    }
  }
}

// calculates the scoring function. This is normally calculated as soon as a new
// indidivual is created. But we need to make a call to this function when the
// slope is changed, to make sure all the individuals get the new value for
// their scoring function before continuing with the GA
void GPPopulation::Eval(double hitlimit, bool function) {
  for (int i = 0; i < popsize; i++)
    ff->CalculateFitness(pop[i], ittrain, sfttrain, function);
}

// Size 2 probabilistic tournament selection, where the winner of the
// tournament is accepted with a probability of tp. Otherwise the
// loser is accepted.
GPGenomePtr GPPopulation::TSelect(double tp) const {
  int ind1, ind2;
  ind1 = m_rand.GetRandomInt(popsize);
  do {
    ind2 = m_rand.GetRandomInt(popsize);
  } while (ind1 == ind2);
  double r = m_rand.GetRandom01();
  if (pop[ind1]->GetFitness() > pop[ind2]->GetFitness())
    // winner is ind1
    if (r < tp)
      return (pop[ind1]);
    else
      return (pop[ind2]);
  else
      // winner is ind2
      if (r < tp)
    return (pop[ind2]);
  else
    return (pop[ind1]);
}

GPGenomePtr GPPopulation::Select(std::string selector) const {
  return TSelect(0.7);
}

void GPPopulation::GAstep(std::string selector, double pcross, double pmut,
                          double mean, double variance, double hitlimit,
                          bool function) {
  GPGenomePtr mother, father;
  for (int i = 0; i < nrepl - 1; i++) {
    mother = Select(selector);
    father = Select(selector);
    int j = 0;
    while (father == mother) {
      father = Select(selector);
      if (j > 100)
        // only if all or almost all individuals of the population are
        // copies of the same one
        throw Error(_WHERE_, "Population failure - not enough diversity");
      j++;
    }
    if (m_rand.GetRandom01() < pcross)
      newpop[i]->UniformCrossover(*mother, *father);
    newpop[i]->Mutate(pmut);
  }
  *(newpop[nrepl - 1]) = *(pop[bestInd]);

  // calculate the objective values and
  // sort newpop  newpop = GPGenomeList (p.newpop.size());
  double bestFitness = -FLT_MAX, fit;
  for (unsigned int i = 0; i < newpop.size(); i++) {
    *(pop[i]) = *(newpop[i]);
    fit = ff->CalculateFitness(pop[i], ittrain, sfttrain, function);
    if (fit > bestFitness) {
      bestFitness = fit;
      bestInd = i;
    }
  }
}

void GPPopulation::EPstep(std::string selector, double pcross, double pmut,
                          double mean, double variance, double hitlimit,
                          bool function) {
  *(newpop[0]) = *(pop[bestInd]);
  for (int i = 1; i < popsize; i++) {
    *(newpop[i]) = *(pop[bestInd]);
    LOG_F(1, "new ind");
    newpop[i]->Mutate(pmut);
  }
  // calculate the objective values and
  // sort newpop  newpop = GPGenomeList (p.newpop.size());
  double bestFitness = -FLT_MAX, fit;
  for (unsigned int i = 0; i < newpop.size(); i++) {
    *(pop[i]) = *(newpop[i]);
    fit = ff->CalculateFitness(pop[i], ittrain, sfttrain, function);
    if (fit >= bestFitness) {
      bestFitness = fit;
      bestInd = i;
    }
  }
}

// the worst individual from the population are replaced by the best from p
GPGenomePtr GPPopulation::Best() const { return (pop[bestInd]); }

std::ostream &GPPopulation::Print(std::ostream &s) const {
  for (unsigned int i = 0; i < pop.size(); i++) {
    pop[i]->Print(s);
    s << pop[i]->GetFitness() << std::endl;
  }
  return s;
}

std::ostream &operator<<(std::ostream &s, const GPPopulation &p) {
  return p.Print(s);
}

// calculates the scaled objectives
// of the current population. To calculate the scaled objectives, we
// use sigma truncating scaling as explained in Golberg's book,
// page 124
void GPPopulation::ScaleFitness() {
  total = 0;
  // calculate the new Fitnesses and the average
  for (int i = 0; i < popsize; i++) {
    total += pop[i]->GetFitness();
  }
  ave = total / popsize;
  // calculate variance
  double ss = 0.0;
  for (int i = 0; i < popsize; i++) {
    ss += (pop[i]->GetFitness() - ave) * (pop[i]->GetFitness() - ave);
  }
  // standard deviation
  stdev = std::sqrt(ss / (popsize - 1));
  // calculate scaled values using sigma truncation
  // Goldberg page 124
  for (int i = 0; i < popsize; i++) {
    double f = pop[i]->GetFitness() - (ave - c * stdev);
    if (f < 0.0)
      f = 0.0;
    pop[i]->SetFitness(f);
  }
}

GPGenomePtr GPPopulation::RwSelect() const {
  double cutoff;
  int i, upper, lower;

  cutoff = m_rand.GetRandom01();
  lower = 0;
  upper = popsize - 1;
  while (upper >= lower) {
    i = lower + (upper - lower) / 2;
    if (psum[i] > cutoff)
      upper = i - 1;
    else
      lower = i + 1;
  }
  // make sure lower is a number between 0 and size - 1
  lower = (lower > (popsize - 1)) ? (popsize - 1) : lower;
  lower = (lower < 0) ? 0 : lower;
  return (pop[lower]);
}

void GPPopulation::Swap(GPPopulationPtr p, int nmigr) {
  for (int i = 0; i < nmigr; i++) {
    GPGenomeListIter where = find(pop.begin(), pop.end(), p->pop[i]);
    // if (!where)
    if (where == pop.end())
      // copy it only if it is not already there
      this->pop[popsize - nmigr + i] = p->pop[i];
  }
}

GPGenomePtr GPPopulation::RkSelect() const {
  double bias = 1.9;
  double index;
  int max = popsize; // index will be between 0 and popsize - 1
  index = max *
          (bias -
           std::sqrt(bias * bias - 4.0 * (bias - 1) * m_rand.GetRandom01())) /
          2.0 / (bias - 1);
  if (index < 0) {
    LOG_F(1, "index={} max={}", index, max);
    index = 0.0;
  }
  // DM 25/04/05 - avoid compiler warning by explicitly converting double to int
  int intIndex = (int)std::floor(index);
  return (pop[intIndex]);
}

/*
Bool GPPopulation::Gen_eq(GPGenome* x, GPGenome* y)
{
    return true; //x->GenEqual(*x, *y);
}

void GPPopulation::QSort(GPGenomeList &p)
{
     sort(p.begin(), p.end(), GenCmp());
}

void GPPopulation::QSort()
{
    sort(pop.begin(), pop.end(), GenCmp());
}

void GPPopulation::MergePops()
{
    GPGenomeList mpop(popsize + newpop.size(), pop[0]);
    LOG_F(1, "pop");
    for (unsigned int k1 = 0 ; k1 < popsize ; k1++)
        LOG_F(1, "{} {} {}", *(pop[k1]), pop[k1]->GetFitness(),
              pop[k1]->GetFitness());
    LOG_F(1, "newpop");
    for (GPGenomeListIter p = newpop.begin() ; p != newpop.end() ; p++)
        LOG_F(1, "{} {} {}", (*p), (*p)->GetFitness(), (*p)->GetFitness());
    merge(pop.begin(), pop.end(), newpop.begin(), newpop.end(),
          mpop.begin(), GenCmp());
    GPGenomeListIter end = unique(mpop.begin(), mpop.end(), Gen_eq);
            // uses operator== , Gen_eq);
    mpop.erase(end, mpop.end());
    LOG_F(1, "mpop");
    for (unsigned int k1 = 0 ; k1 < mpop.size() ; k1++)
        LOG_F(1, "{} {} {}", mpop[k1], mpop[k1]->GetFitness(),
              mpop[k1]->GetFitness());
    pop.clear();
    end = (mpop.size() > popsize) ? (mpop.begin() + popsize) : mpop.end();
    copy(mpop.begin(), end, back_inserter(pop));
    LOG_F(1, "copied pop");
    for (Int k1 = 0 ; k1 < popsize ; k1++)
        LOG_F(1, "{} {} {}", *(pop[k1]), pop[k1]->GetFitness(),
              pop[k1]->GetFitness());
}*/
