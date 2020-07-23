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

// GPGenome Class. GPGenome represents a genome in the CGP

#ifndef _RBT_GPGENOME_H_
#define _RBT_GPGENOME_H_

#include "RbtGPChromosome.h"
#include "RbtGPTypes.h"
#include "RbtRand.h"
// #include <cmath>

namespace rxdock {

namespace geneticprogram {

class RbtGPGenome {
public:
  static std::string _CT;
  ///////////////////
  // Constructors
  ///////////////////
  RbtGPGenome();
  RbtGPGenome(const RbtGPGenome &);
  RbtGPGenome(std::string);
  RbtGPGenome(std::istream &);

  RbtGPGenome &operator=(const RbtGPGenome &);

  ///////////////////
  // Destructor
  //////////////////
  virtual ~RbtGPGenome();
  static void SetStructure(int, int, int, int, int, int, int, int);
  static int GetNIP() { return npi; }
  static void SetNIP(int n) { npi = n; }
  static int GetNIF() { return nfi; }
  static int GetNN() { return nn; }
  static int GetNO() { return no; }
  static int GetNSFI() { return nsfi; }
  static void SetNSFI(int n) { nsfi = n; }
  RbtGPChromosomePtr GetChrom() const { return (new RbtGPChromosome(*chrom)); }
  void Initialise();
  void Mutate(double);
  void SetFitness(double f) { fitness = f; }
  double GetFitness() const { return fitness; }
  void UniformCrossover(const RbtGPGenome &, const RbtGPGenome &);
  void Crossover(RbtGPGenome &);
  friend std::ostream &operator<<(std::ostream &s, const RbtGPGenome &p);
  std::ostream &Print(std::ostream &) const;

protected:
private:
  void MutateGene(int);

  static int npi, nfi, nsfi, no, nn, nf, nr, nc, l;
  // defines structures and constraints of each genome
  RbtRand &m_rand;
  double fitness;           // scaled score. For now same than score
  RbtGPChromosomePtr chrom; // list of integers that represent a chrom
};

// Useful typedefs
typedef SmartPtr<RbtGPGenome> RbtGPGenomePtr;        // Smart pointer
typedef std::vector<RbtGPGenomePtr> RbtGPGenomeList; // Vector of smart pointers
typedef RbtGPGenomeList::iterator RbtGPGenomeListIter;
typedef RbtGPGenomeList::const_iterator RbtGPGenomeListConstIter;

std::ostream &operator<<(std::ostream &s, const RbtGPGenome &p);

} // namespace geneticprogram

} // namespace rxdock

#endif //_RbtGPGenome_H_
