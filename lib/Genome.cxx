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

#include "Genome.h"
#include "BaseSF.h"
#include "ChromElement.h"

using namespace rxdock;

std::string Genome::_CT("Genome");

Genome::Genome(ChromElement *pChr)
    : m_chrom(pChr->clone()), m_score(0.0), m_RWFitness(0.0) {
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

Genome::Genome(const Genome &g)
    : m_chrom((g.m_chrom)->clone()), m_score(g.m_score),
      m_RWFitness(g.m_RWFitness) {
  _RBTOBJECTCOUNTER_COPYCONSTR_(_CT);
}

Genome &Genome::operator=(const Genome &g) {
  if (&g != this) {
    m_chrom = (g.m_chrom)->clone();
    m_score = g.m_score;
    m_RWFitness = g.m_RWFitness;
  }
  return *this;
}

Genome::~Genome() {
  delete m_chrom;
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

Genome *Genome::clone() const { return new Genome(*this); }

void Genome::SetScore(BaseSF *pSF) {
  if (pSF != nullptr) {
    m_chrom->SyncToModel();
    m_score = -pSF->Score();
  } else {
    m_score = 0.0;
  }
  SetRWFitness(0.0, 0.0);
}

double Genome::SetRWFitness(double sigmaOffset, double partialSum) {
  // Apply sigma truncation to the raw score
  m_RWFitness = std::max(0.0, GetScore() - sigmaOffset);
  // The fitness value we store is the partial sum from previous Genome
  // elements Not normalised at this point.
  m_RWFitness += partialSum;
  return m_RWFitness;
}

void Genome::NormaliseRWFitness(double total) {
  if (total > 0.0) {
    m_RWFitness /= total;
  }
}

void Genome::Print(std::ostream &s) const {
  s << *m_chrom << std::endl;
  s << "Score: " << GetScore() << "; RWFitness: " << GetRWFitness()
    << std::endl;
}

std::ostream &operator<<(std::ostream &s, const Genome &g) {
  g.Print(s);
  return s;
}
