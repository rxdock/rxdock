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

// Genome for roulette wheel selection.
// Manages a chromosome element, an associated raw score from a scoring function
// and a scaled fitness value for roulette wheel selection
#ifndef _RBT_GENOME_H_
#define _RBT_GENOME_H_

#include "rxdock/ChromElement.h"

namespace rxdock {

class BaseSF;

class Genome {
public:
  static std::string _CT;
  // Sole constructor accepting an existing chromosome.
  // A clone of the existing chromosome is stored,
  // and is deleted in the Genome destructor
  Genome(ChromElement *pChr);
  // Copy constructor and Assignment operator.
  // The underlying chromosome is cloned during the copy.
  Genome(const Genome &);
  Genome &operator=(const Genome &);
  Genome *clone() const;
  virtual ~Genome();

  // Gets the underlying chromosome
  ChromElement *GetChrom() const { return m_chrom; }

  // Sets the raw score from the scoring function.
  // Model coordinates are updated to reflect the current chromosome element
  // values. Scores are negated so that higher positive scores are better. pSF
  // is a pointer to a scoring function object. If pSF is null, a zero score is
  // set.
  void SetScore(BaseSF *pSF);
  // Gets the stored raw score (without re-evaluation of the scoring function).
  double GetScore() const { return m_score; }

  // ROULETTE WHEEL SELECTION METHODS
  // The following three methods should be used only in the context of
  // genome lists sorted into ascending order by score.
  // See Population::EvaluateRWFitness for example of use.
  //
  // Sets the fitness value used for roulette wheel selection
  // sigmaOffset is the calculated offset for scaling using sigma truncation
  //(Goldberg page 124)
  // partialSum is the partial sum of all previous fitness values in the sorted
  // genome list.
  double SetRWFitness(double sigmaOffset, double partialSum);
  // Normalises the roulette wheel fitness values to lie between 0 and 1
  // total is the total sum of absolute fitness values returned by the last call
  // to SetRWFitness
  void NormaliseRWFitness(double total);
  // Gets the fitness value for roulette wheel selection
  double GetRWFitness() const { return m_RWFitness; }

  // Tests equality based on chromosome element values.
  // Does not take score into account.
  bool Equals(const Genome &g, double threshold) const {
    return m_chrom->Equals(*(g.m_chrom), threshold);
  }
  friend bool operator==(const Genome &g1, const Genome &g2) {
    return g1.Equals(g2, ChromElement::_THRESHOLD);
  }

  void Print(std::ostream &) const;
  friend std::ostream &operator<<(std::ostream &s, const Genome &g);

private:
  Genome(); // Default constructor disabled

private:
  /////////////////////
  // Private data
  /////////////////////
  ChromElement *m_chrom;
  double m_score;     // raw value of the scoring function
  double m_RWFitness; // scaled value of the raw score suitable for use
                      // with roulette wheel selection
};

// Useful typedefs
typedef SmartPtr<Genome> GenomePtr;        // Smart pointer
typedef std::vector<GenomePtr> GenomeList; // Vector of smart pointers
typedef GenomeList::iterator GenomeListIter;
typedef GenomeList::const_iterator GenomeListConstIter;

// Compare class to compare different genomes inside a population
// This is used to sort the genomes in a population depending
// of the value of their scoring function.
class GenomeCmp_Score {
public:
  bool operator()(const Genome *pG1, const Genome *pG2) const {
    return (pG1->GetScore() > pG2->GetScore());
  }
};

// STL binary predicate to compare two Genome* for equality.
// Tests equality based on chromosome element values.
// Does not take score into account.
// Uses operator== which ultimately calls ChromElement::Equals()
class isGenome_eq : public std::function<bool(Genome *, Genome *)> {
private:
  double m_threshold; // equality threshold
public:
  explicit isGenome_eq(double threshold) : m_threshold(threshold) {}
  bool operator()(const Genome *pG1, const Genome *pG2) const {
    return pG1->Equals(*pG2, m_threshold);
  }
};

} // namespace rxdock

#endif //_RBT_GENOME_H_
