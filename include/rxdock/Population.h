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

// Manages a population of genomes.
// Main public method (GAstep) performs one iteration of a GA.
#ifndef _RBTPOPULATION_H_
#define _RBTPOPULATION_H_

#include "rxdock/Error.h"
#include "rxdock/Genome.h"

#include <nlohmann/json.hpp>

using json = nlohmann::json;

namespace rxdock {

class BaseSF; // forward definition

class Population {
public:
  static const std::string _CT;
  // Sole constructor to create a randomised genome population of a fixed size.
  // pChr is the seed chromosome to clone to create each genome.
  // size is the population size to create.
  // pSF is the scoring function used to rank the genomes.
  // Notes:
  // 1) Each genome is randomised by calling the Randomise() method
  // of the underlying chromosome.
  // 2) Scoring function scores are calculated.
  // 3) Genomes are sorted by score.
  // 4) Roulette wheel fitness values are calculated.
  // 5) Model coords are updated to match the fittest chromosome
  // An BadArgument error is thrown if size is <=0, or if pChr or pSF is
  // null.
  RBTDLL_EXPORT Population(ChromElement *pChr, int size, BaseSF *pSF);
  virtual ~Population();

  friend void to_json(json &j, const Population &population);
  friend void from_json(const json &j, Population &population);

  // Gets the maximum size of the population as defined in the constructor.
  int GetMaxSize() const { return m_size; }
  // Gets the actual size of the population (may be < GetMaxSize())
  int GetActualSize() const { return m_pop.size(); }
  // Gets the best genome in the sorted population (element zero).
  RBTDLL_EXPORT GenomePtr Best() const;
  // Gets the average raw score across entire population
  double GetScoreMean() const { return m_scoreMean; }
  // Gets the raw score variance across entire population
  double GetScoreVariance() const { return m_scoreVariance; }
  // Gets the vector of genomes in the population.
  const GenomeList &GetGenomeList() const { return m_pop; }

  // Gets the scoring function used for ranking genomes
  BaseSF *GetSF() const { return m_pSF; }
  // Sets the scoring function used for ranking genomes
  // Forces a recalculation of the genome scores, and the roulette wheel fitness
  // values SetSF should be called whenever the scoring function parameters have
  // changed e.g. in between GA stages. An BadArgument error is thrown if pSF
  // is null. Model coords are updated to match the fittest chromosome
  void SetSF(BaseSF *pSF);

  // Main method for performing a GA iteration
  RBTDLL_EXPORT void
  GAstep(int nReplicates,    // Number of new genomes to create in the iteration
         double relStepSize, // Relative step size for chromosome mutations
         double equalityThreshold, // Equality threshold for genomes
         double pcross,            // Probability of crossover
         bool xovermut, // if true, perform cauchy mutation following crossover
         bool cmutate   // true=cauchy mutations, false=regular mutations
  );
  RBTDLL_EXPORT GenomePtr RouletteWheelSelect() const;

  void Print(std::ostream &) const;
  friend std::ostream &operator<<(std::ostream &, const Population &);

private:
  // Merges the new individuals created into the main population
  // Duplicate genomes are removed (based on equality of chromosome elements,
  // not scores)
  void MergeNewPop(GenomeList &newPop, double equalityThreshold);
  void EvaluateRWFitness();
  Population(const Population &);            // Disable
  Population &operator=(const Population &); // Disable

  GenomeList m_pop;       // The population of genomes
  unsigned int m_size;    // The maximum size of the population
  double m_c;             // Sigma Truncation Multiplier
  BaseSF *m_pSF;          // The scoring function
  Rand &m_rand;           // reference to the singleton random number generator
  double m_scoreMean;     // the average raw score across all genomes
  double m_scoreVariance; // the variance of raw scores across all genomes
};

void to_json(json &j, const Population &population);
void from_json(const json &j, Population &population);

typedef SmartPtr<Population> PopulationPtr;

} // namespace rxdock

#endif //_RBTPOPULATION_H_
