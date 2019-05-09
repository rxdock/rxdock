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

// Population class. Defines a population of CGP genomes that represent
// functions.

#ifndef _RBTGPPOPULATION_H_
#define _RBTGPPOPULATION_H_

#include "RbtError.h"
#include "RbtGPFFGold.h"
#include "RbtGPFitnessFunction.h"
#include "RbtGPGenome.h"

class RbtGPPopulation {
public:
  static std::string _CT;
  ///////////////////
  // Constructors
  ///////////////////
  RbtGPPopulation(int, int, RbtGPFitnessFunctionPtr, RbtReturnTypeArray &,
                  RbtReturnTypeArray &);
  void SetTrainingSets(RbtReturnTypeArray &it, RbtReturnTypeArray sf) {
    ittrain = it;
    sfttrain = sf;
  }

  ///////////////////
  // Destructors
  ///////////////////
  RbtGPPopulation(const RbtGPPopulation &);
  virtual ~RbtGPPopulation();
  void Initialise(double, bool);
  void ScaleFitness();
  void Eval(double, bool);
  RbtGPGenomePtr Select(std::string) const;

  void SelectionUpdate(std::string);

  void GAstep(std::string, double, double, double, double, double,
              bool) throw(RbtError);
  void EPstep(std::string, double, double, double, double, double,
              bool) throw(RbtError);
  RbtGPGenomePtr Best() const;
  std::ostream &Print(std::ostream &) const;
  friend std::ostream &operator<<(std::ostream &, const RbtGPPopulation &);
  void Swap(SmartPtr<RbtGPPopulation>, int);
  void QSort(RbtGPGenomeList &);
  void QSort();
  void MergePops();
  static bool Gen_eq(RbtGPGenome *, RbtGPGenome *);
  int GetSize();
  int GetNrepl();

private:
  RbtGPGenomePtr RwSelect() const;
  RbtGPGenomePtr RkSelect() const;
  RbtGPGenomePtr TSelect(double) const;
  RbtGPGenomeList pop;    // population: array to pointers of RbtGPGenomes
  RbtGPGenomeList newpop; // array where the new individuals created are stored
  int popsize;            // size of population
  int nrepl;              // Number of new individuals that get
                          // created at each generation
  double total;           // total addition of the score values
  double ave;             // average
  double stdev;           // standard deviation
  double c;               // Sigma Truncation Multiplier
  double *psum;           // stores the partial sums of the fitness values. Used
                          // when selecting individuals from the population
  RbtRand &m_rand;        // keep a reference to the singleton
                          // random number generator
  int bestInd;            // keeps the index of the best individual
  RbtGPFitnessFunctionPtr ff;
  RbtReturnTypeArray ittrain, sfttrain;
};
typedef SmartPtr<RbtGPPopulation> RbtGPPopulationPtr;
typedef vector<RbtGPPopulationPtr> RbtGPPopulationList;
typedef RbtGPPopulationList::iterator RbtGPPopulationListIter;

// Compare class to compare different populations. Use to find the
// genome with best scoring function between all the populations
class PopCmp {
public:
  bool operator()(const RbtGPPopulationPtr &x,
                  const RbtGPPopulationPtr &y) const {
    return (x->Best()->GetFitness() < y->Best()->GetFitness());
  }
};

#endif //_RbtGPPopulation_H_
