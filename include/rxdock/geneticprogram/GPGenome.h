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

#include "rxdock/Rand.h"
#include "rxdock/geneticprogram/GPChromosome.h"
#include "rxdock/geneticprogram/GPTypes.h"
// #include <cmath>

namespace rxdock {

namespace geneticprogram {

class GPGenome {
public:
  static std::string _CT;
  ///////////////////
  // Constructors
  ///////////////////
  GPGenome();
  GPGenome(const GPGenome &);
  GPGenome(std::string);
  GPGenome(std::istream &);

  GPGenome &operator=(const GPGenome &);

  ///////////////////
  // Destructor
  //////////////////
  virtual ~GPGenome();
  static void SetStructure(int, int, int, int, int, int, int, int);
  static int GetNIP() { return npi; }
  static void SetNIP(int n) { npi = n; }
  static int GetNIF() { return nfi; }
  static int GetNN() { return nn; }
  static int GetNO() { return no; }
  static int GetNSFI() { return nsfi; }
  static void SetNSFI(int n) { nsfi = n; }
  GPChromosomePtr GetChrom() const { return (new GPChromosome(*chrom)); }
  void Initialise();
  void Mutate(double);
  void SetFitness(double f) { fitness = f; }
  double GetFitness() const { return fitness; }
  void UniformCrossover(const GPGenome &, const GPGenome &);
  void Crossover(GPGenome &);
  friend std::ostream &operator<<(std::ostream &s, const GPGenome &p);
  std::ostream &Print(std::ostream &) const;

protected:
private:
  void MutateGene(int);

  static int npi, nfi, nsfi, no, nn, nf, nr, nc, l;
  // defines structures and constraints of each genome
  Rand &m_rand;
  double fitness;        // scaled score. For now same than score
  GPChromosomePtr chrom; // list of integers that represent a chrom
};

// Useful typedefs
typedef SmartPtr<GPGenome> GPGenomePtr;        // Smart pointer
typedef std::vector<GPGenomePtr> GPGenomeList; // Vector of smart pointers
typedef GPGenomeList::iterator GPGenomeListIter;
typedef GPGenomeList::const_iterator GPGenomeListConstIter;

std::ostream &operator<<(std::ostream &s, const GPGenome &p);

} // namespace geneticprogram

} // namespace rxdock

#endif //_GPGenome_H_
