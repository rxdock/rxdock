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

#include "rxdock/Population.h"
#include "rxdock/Debug.h"
#include "rxdock/DockingError.h"
#include <algorithm>

using namespace rxdock;

const std::string Population::_CT = "Population";

Population::Population(ChromElement *pChr, int size, BaseSF *pSF)
    : m_size(size), m_c(2.0), m_pSF(pSF), m_rand(GetRandInstance()),
      m_scoreMean(0.0), m_scoreVariance(0.0) {
  if (pChr == nullptr) {
    throw BadArgument(
        _WHERE_, "Null chromosome element passed to Population constructor");
  } else if (size <= 0) {
    throw BadArgument(_WHERE_, "Population size must be positive (non-zero)");
  }
  // Create a random population
  m_pop.reserve(m_size);
  for (unsigned int i = 0; i < m_size; ++i) {
    // The Genome constructor clones the chromosome to create an independent
    // copy
    GenomePtr genome = new Genome(pChr);
    genome->GetChrom()->Randomise();
    m_pop.push_back(genome);
  }
  // Calculate the scores and evaluate roulette wheel fitness
  SetSF(m_pSF);
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

Population::~Population() { _RBTOBJECTCOUNTER_DESTR_(_CT); }

// Sets the scoring function used for ranking genomes
void Population::SetSF(BaseSF *pSF) {
  if (pSF == nullptr) {
    throw BadArgument(_WHERE_, "Null scoring function passed to SetSF");
  }
  m_pSF = pSF;
  for (GenomeListIter iter = m_pop.begin(); iter != m_pop.end(); ++iter) {
    (*iter)->SetScore(m_pSF);
  }
  std::stable_sort(m_pop.begin(), m_pop.end(), GenomeCmp_Score());
  EvaluateRWFitness();
}

void Population::GAstep(int nReplicates, double relStepSize,
                        double equalityThreshold, double pcross, bool xovermut,
                        bool cmutate) {
  if (nReplicates <= 0) {
    throw BadArgument(_WHERE_, "nReplicates must be positive (non-zero)");
  }
  GenomeList newPop;
  newPop.reserve(nReplicates);
  for (int i = 0; i < nReplicates / 2; i++) {
    GenomePtr mother = RouletteWheelSelect();
    GenomePtr father = RouletteWheelSelect();
    // Check that mother and father are not the same genome
    // The check is on the pointers, not that the chromosomes are near-equal
    // If we repeatedly get the same genomes selected, this must mean
    // the population lacks diversity
    int j = 0;
    while (father == mother) {
      father = RouletteWheelSelect();
      if (j > 100)
        throw DockingError(_WHERE_,
                           "Population failure - not enough diversity");
      j++;
    }
    GenomePtr child1 = new Genome(*mother);
    GenomePtr child2 = new Genome(*father);
    // Crossover
    if (m_rand.GetRandom01() < pcross) {
      Crossover(father->GetChrom(), mother->GetChrom(), child1->GetChrom(),
                child2->GetChrom());
      // Cauchy mutation following crossover
      if (xovermut) {
        child1->GetChrom()->CauchyMutate(0.0, relStepSize);
        child2->GetChrom()->CauchyMutate(0.0, relStepSize);
      }
    }
    // Mutation
    else {
      // Cauchy mutation
      if (cmutate) {
        child1->GetChrom()->CauchyMutate(0.0, relStepSize);
        child2->GetChrom()->CauchyMutate(0.0, relStepSize);
      }
      // Regular mutation
      else {
        child1->GetChrom()->Mutate(relStepSize);
        child2->GetChrom()->Mutate(relStepSize);
      }
    }
    newPop.push_back(child1);
    newPop.push_back(child2);
  }
  // check if one more is needed (odd nReplicates).
  if (nReplicates % 2) {
    GenomePtr mother = RouletteWheelSelect();
    GenomePtr child = new Genome(*mother);
    child->GetChrom()->CauchyMutate(0.0, relStepSize);
    newPop.push_back(child);
  }
  MergeNewPop(newPop, equalityThreshold);
  EvaluateRWFitness();
}

GenomePtr Population::Best() const {
  return m_pop.empty() ? GenomePtr() : m_pop.front();
}

void Population::Print(std::ostream &s) const {
  s << m_pop.size() << std::endl;
  for (GenomeListConstIter iter = m_pop.begin(); iter != m_pop.end(); ++iter) {
    (*iter)->Print(s);
  }
}

std::ostream &operator<<(std::ostream &s, const Population &p) {
  p.Print(s);
  return s;
}

GenomePtr Population::RouletteWheelSelect() const {
  double cutoff = m_rand.GetRandom01();
  int size = m_pop.size();
  int lower = 0;
  int upper = size - 1;
  while (upper >= lower) {
    int i = lower + (upper - lower) / 2;
    if (m_pop[i]->GetRWFitness() > cutoff)
      upper = i - 1;
    else
      lower = i + 1;
  }
  // make sure lower is a number between 0 and size - 1
  lower = std::min(size - 1, lower);
  lower = std::max(0, lower);
  return (m_pop[lower]);
}

void Population::MergeNewPop(GenomeList &newPop, double equalityThreshold) {
  // Assume newPop needs scoring and sorting
  for (GenomeListIter iter = newPop.begin(); iter != newPop.end(); ++iter) {
    (*iter)->SetScore(m_pSF);
  }
  std::stable_sort(newPop.begin(), newPop.end(), GenomeCmp_Score());

  GenomeList mergedPop;
  mergedPop.reserve(m_pop.size() + newPop.size());
  // Merge pops by score
  std::merge(m_pop.begin(), m_pop.end(), newPop.begin(), newPop.end(),
             std::back_inserter(mergedPop), GenomeCmp_Score());
  // Remove neighbouring duplicates by equality of chromosome element values
  GenomeListIter end = std::unique(mergedPop.begin(), mergedPop.end(),
                                   isGenome_eq(equalityThreshold));
  mergedPop.erase(end, mergedPop.end());
  m_pop.clear();
  end = (mergedPop.size() > m_size) ? (mergedPop.begin() + m_size)
                                    : mergedPop.end();
  std::copy(mergedPop.begin(), end, back_inserter(m_pop));
}

void Population::EvaluateRWFitness() {
  // Determine mean and variance of true scores
  double sum(0.0);
  double sumSq(0.0);
  for (GenomeListConstIter iter = m_pop.begin(); iter != m_pop.end(); ++iter) {
    double score = (*iter)->GetScore();
    sum += score;
    sumSq += score * score;
  }
  double popSize = static_cast<double>(m_pop.size());
  m_scoreMean = sum / popSize;
  m_scoreVariance = (sumSq / popSize) - (m_scoreMean * m_scoreMean);
  double sigma = std::sqrt(m_scoreVariance);
  // calculate scaled fitness values using sigma truncation
  // Goldberg page 124
  double offset = m_scoreMean - m_c * sigma;
  double partialSum = 0.0;
  // First set the unnormalised fitness values
  for (GenomeListIter iter = m_pop.begin(); iter != m_pop.end(); ++iter) {
    partialSum = (*iter)->SetRWFitness(offset, partialSum);
  }
  // Now normalise so that fitness values run from 0 to 1, ready for
  // roulette wheel selection
  for (GenomeListIter iter = m_pop.begin(); iter != m_pop.end(); ++iter) {
    (*iter)->NormaliseRWFitness(partialSum);
  }
}
