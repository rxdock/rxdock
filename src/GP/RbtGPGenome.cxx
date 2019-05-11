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

#include "RbtGPGenome.h"
#include "RbtDebug.h"
#include "RbtGPParser.h"
#include <fstream>
#include <sstream>

std::string RbtGPGenome::_CT("RbtGPGenome");
int RbtGPGenome::npi;
int RbtGPGenome::nfi;
int RbtGPGenome::nsfi;
int RbtGPGenome::no;
int RbtGPGenome::nn;
int RbtGPGenome::nf;
int RbtGPGenome::nr;
int RbtGPGenome::nc;
int RbtGPGenome::l;

// Constructors
RbtGPGenome::RbtGPGenome() : m_rand(Rbt::GetRbtRand()) {
  chrom = new RbtGPChromosome(npi, nfi, nn, no, nr, nc);
  fitness = 0.0;
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtGPGenome::RbtGPGenome(const RbtGPGenome &g) : m_rand(Rbt::GetRbtRand()) {
  chrom = new RbtGPChromosome(*(g.chrom));
  fitness = g.fitness;
  _RBTOBJECTCOUNTER_COPYCONSTR_(_CT);
}

RbtGPGenome::RbtGPGenome(std::istream &in) : m_rand(Rbt::GetRbtRand()) {
  long int seed;
  in >> seed;
  m_rand.Seed(seed);
  // Get structure
  in >> npi >> nfi >> nsfi >> no >> nf >> nr >> nc >> l;
  nn = nr * nc; // number of nodes
  chrom = new RbtGPChromosome(npi, nfi, nn, no, nr, nc);
  in >> *chrom;
  fitness = 0.0;
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtGPGenome::RbtGPGenome(std::string str) : m_rand(Rbt::GetRbtRand()) {
  std::istringstream ist(str);
  // Get structure
  ist >> npi >> nfi >> nsfi >> no >> nf >> nr >> nc >> l;
  nn = nr * nc; // number of nodes
  chrom = new RbtGPChromosome(npi, nfi, nn, no, nr, nc);
  ist >> *chrom;
  fitness = 0.0;
}

RbtGPGenome &RbtGPGenome::operator=(const RbtGPGenome &g) {
  if (this != &g) {
    *chrom = *(g.chrom);
    fitness = g.fitness;
  }
  return *this;
}

///////////////////
// Destructor
//////////////////
RbtGPGenome::~RbtGPGenome() { _RBTOBJECTCOUNTER_DESTR_(_CT); }

void RbtGPGenome::SetStructure(int inpi, int infi, int insfi, int ino, int inf,
                               int inr, int inc, int il) {
  npi = inpi;   // number of program inputs
  nfi = infi;   // number of inputs per function
  nsfi = insfi; // number of inputs needed to calculate the SF
  no = ino;     // number of outputs of the program.
                // Assumption: functions return one output
  nf = inf;     // number of functions
  nr = inr;     // number of rows
  nc = inc;     // number of columns
  nn = nr * nc; // number of nodes
  l = il;       // connectivity level, i.e., how many previous columns
                // of cells may have their outputs connected to a
                // node in the current column
}

void RbtGPGenome::Initialise() {
  int s = chrom->size();
  int g = 0, max, min;
  for (int i = 0; i < nc; i++) {
    for (int j = 0; j < nr; j++) {
      for (int k = 0; k < nfi; k++) {
        max = npi + i * nr;
        if (i < l)
          (*chrom)[g++] = m_rand.GetRandomInt(max);
        else {
          min = npi + (i - l) * nr;
          (*chrom)[g++] = m_rand.GetRandomInt(max - min) + min;
        }
      }
      // node's function
      int f = m_rand.GetRandomInt(nf + 1);
      if (f == nf) // constant Int
      {
        int a, b;
        a = m_rand.GetRandomInt(200) - 100;
        b = m_rand.GetRandomInt(10) - 5;
        double r =
            (a / 10.0) * std::pow(10.0, b); // m_rand.GetRandom01() * 20 - 10;
        chrom->SetConstant(r, g);
      }
      (*chrom)[g++] = f;
    }
  }
  // outputs of program
  min = npi + (nc - l) * nr;
  max = npi + nc * nr;
  for (int i = 0; i < no; i++)
    (*chrom)[g++] = m_rand.GetRandomInt(max - min) + min;
  fitness = 0.0;
}

void RbtGPGenome::MutateGene(int i) {
  int min, max, column;
  // is it part of the nodes?
  if (i < (nr * nc * (nfi + 1))) {
    // is it an input of a function?
    if (((i + 1) % (nfi + 1)) != 0) {
      column = i / (nr * (nfi + 1));
      max = npi + column * nr;
      if (column < l)
        (*chrom)[i] = m_rand.GetRandomInt(max);
      else {
        min = npi + (column - l) * nr;
        (*chrom)[i] = m_rand.GetRandomInt(max - min) + min;
      }
    } else // it is a function
    {
      int f = m_rand.GetRandomInt(nf + 1); // m_rand.GetRandomInt(nf);
      if (f == nf)                         // constant Int
      {
        int a, b;
        a = m_rand.GetRandomInt(200) - 100;
        b = m_rand.GetRandomInt(10) - 5;
        double r =
            (a / 10.0) * std::pow(10.0, b); // m_rand.GetRandom01() * 20 - 10;
        chrom->SetConstant(r, i);
      } else if ((*chrom)[i] == nf) {
        chrom->ResetConstant(i);
      }
      (*chrom)[i] = f;
    }
  } else // it is an output
  {
    min = npi + (nc - l) * nr;
    max = npi + nc * nr;
    (*chrom)[i] = m_rand.GetRandomInt(max - min) + min;
  }
}

void RbtGPGenome::Mutate(double mutRate) {
  for (int i = 0; i < chrom->size(); i++)
    if (m_rand.GetRandom01() < mutRate) {
      MutateGene(i);
    }
}

void RbtGPGenome::UniformCrossover(const RbtGPGenome &mom,
                                   const RbtGPGenome &dad) {
  int coin;
  for (int i = 0; i < mom.chrom->size(); i++) {
    coin = m_rand.GetRandomInt(2);
    if (coin == 0)
      (*chrom)[i] = (*mom.chrom)[i];
    else
      (*chrom)[i] = (*dad.chrom)[i];
  }
}

void RbtGPGenome::Crossover(RbtGPGenome &) {}

std::ostream &RbtGPGenome::Print(std::ostream &s) const {
  int g = 0;
  // Print structure
  s << npi << " " << nfi << " " << nsfi << " " << no << " " << nf << " " << nr
    << " " << nc << " " << l << std::endl;
  s << *chrom << std::endl;
  RbtGPParser p(npi, nfi, 0, no);
  std::ifstream nstr("descnames", std::ios_base::in);
  s << p.PrintParse(nstr, chrom, true, false);
  return s;
}

std::ostream &operator<<(std::ostream &s, const RbtGPGenome &g) {
  return g.Print(s);
}
