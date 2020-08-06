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

#include "rxdock/Error.h"
#include "rxdock/geneticprogram/GPFFGold.h"
#include "rxdock/geneticprogram/GPFitnessFunction.h"
#include "rxdock/geneticprogram/GPGenome.h"

namespace rxdock {

namespace geneticprogram {

class GPPopulation {
public:
  static std::string _CT;
  ///////////////////
  // Constructors
  ///////////////////
  GPPopulation(int, int, GPFitnessFunctionPtr, ReturnTypeArray &,
               ReturnTypeArray &);
  void SetTrainingSets(ReturnTypeArray &it, ReturnTypeArray sf) {
    ittrain = it;
    sfttrain = sf;
  }

  ///////////////////
  // Destructors
  ///////////////////
  GPPopulation(const GPPopulation &);
  virtual ~GPPopulation();
  void Initialise(double, bool);
  void ScaleFitness();
  void Eval(double, bool);
  GPGenomePtr Select(std::string) const;

  void SelectionUpdate(std::string);

  void GAstep(std::string, double, double, double, double, double, bool);
  void EPstep(std::string, double, double, double, double, double, bool);
  GPGenomePtr Best() const;
  std::ostream &Print(std::ostream &) const;
  friend std::ostream &operator<<(std::ostream &, const GPPopulation &);
  void Swap(SmartPtr<GPPopulation>, int);
  void QSort(GPGenomeList &);
  void QSort();
  void MergePops();
  static bool Gen_eq(GPGenome *, GPGenome *);
  int GetSize();
  int GetNrepl();

private:
  GPGenomePtr RwSelect() const;
  GPGenomePtr RkSelect() const;
  GPGenomePtr TSelect(double) const;
  GPGenomeList pop;    // population: array to pointers of GPGenomes
  GPGenomeList newpop; // array where the new individuals created are stored
  int popsize;         // size of population
  int nrepl;           // Number of new individuals that get
                       // created at each generation
  double total;        // total addition of the score values
  double ave;          // average
  double stdev;        // standard deviation
  double c;            // Sigma Truncation Multiplier
  double *psum;        // stores the partial sums of the fitness values. Used
                       // when selecting individuals from the population
  Rand &m_rand;        // keep a reference to the singleton
                       // random number generator
  int bestInd;         // keeps the index of the best individual
  GPFitnessFunctionPtr ff;
  ReturnTypeArray ittrain, sfttrain;
};
typedef SmartPtr<GPPopulation> GPPopulationPtr;
typedef std::vector<GPPopulationPtr> GPPopulationList;
typedef GPPopulationList::iterator GPPopulationListIter;

// Compare class to compare different populations. Use to find the
// genome with best scoring function between all the populations
class PopCmp {
public:
  bool operator()(const GPPopulationPtr &x, const GPPopulationPtr &y) const {
    return (x->Best()->GetFitness() < y->Best()->GetFitness());
  }
};

} // namespace geneticprogram

} // namespace rxdock

#endif //_GPPopulation_H_
