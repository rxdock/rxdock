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

#include "GPChromosome.h"
#include "Cell.h"
#include "Command.h"

#include <loguru.hpp>

#include <fstream>

using namespace rxdock;
using namespace rxdock::geneticprogram;

std::string GPChromosome::_CT("GPChromosome");

// Constructors
GPChromosome::GPChromosome(int npi, int nfi, int nn, int no, int nr, int nc) {
  chrom = std::vector<int>(nn * (nfi + 1) + no);
  nProgramInputs = npi;
  nFunctionsInputs = nfi;
  int ncells = nn + npi; // number of nodes + inputs
  nProgramOutputs = no;
  nColumns = nc;
  nRows = nr;
  cells = CellList(ncells);
  for (int i = 0; i < ncells; i++)
    cells[i] = new Cell();
}
GPChromosome::GPChromosome(const GPChromosome &c) {
  chrom = std::vector<int>(c.chrom.size());
  cells = CellList(c.cells.size());
  for (unsigned int i = 0; i < cells.size(); i++)
    cells[i] = c.cells[i];
  for (unsigned int i = 0; i < chrom.size(); i++)
    chrom[i] = c.chrom[i];
  nProgramInputs = c.nProgramInputs;
  nFunctionsInputs = c.nFunctionsInputs;
  nProgramOutputs = c.nProgramOutputs;
}

///////////////////
// Destructor
//////////////////
GPChromosome::~GPChromosome() {}

void GPChromosome::Clear() {
  for (unsigned int i = 0; i < cells.size(); i++)
    cells[i]->Clear();
  Command::Clear();
}
GPChromosome &GPChromosome::operator=(const GPChromosome &c) {
  for (unsigned int i = 0; i < cells.size(); i++)
    *(cells[i]) = *(c.cells[i]);
  for (unsigned int i = 0; i < chrom.size(); i++)
    chrom[i] = c.chrom[i];
  nProgramInputs = c.nProgramInputs;
  nFunctionsInputs = c.nFunctionsInputs;
  nProgramOutputs = c.nProgramOutputs;
  return *this;
}

std::istream &rxdock::geneticprogram::operator>>(std::istream &s,
                                                 GPChromosome &p) {
  for (unsigned int i = 0; i < p.chrom.size(); i++)
    s >> p.chrom[i];
  int cell;
  ReturnType value;
  while (!s.eof()) {
    s >> cell >> value;
    p.Cells(cell)->SetConstant(value);
  }
  return s;
}

std::ostream &GPChromosome::Print(std::ostream &s) const {
  int g = 0;
  for (int i = 0; i < nColumns; i++) {
    for (int j = 0; j < nRows; j++) {
      for (int k = 0; k <= nFunctionsInputs; k++)
        s << chrom[g++] << " ";
      s << "\t";
    }
    s << std::endl;
  }
  for (int i = 0; i < nProgramOutputs; i++)
    s << chrom[g++] << " ";
  s << std::endl;
  for (unsigned int i = 0; i < cells.size(); i++)
    if (cells[i]->Constant())
      s << i << "\t" << cells[i]->GetResult() << std::endl;
  return s;
}

std::ostream &rxdock::geneticprogram::operator<<(std::ostream &s,
                                                 const GPChromosome &p) {
  p.Print(s);
  return s;
}

void GPChromosome::SetConstant(ReturnType cte, int idx) {
  int cell = idx / (nFunctionsInputs + 1) + nProgramInputs;
  LOG_F(1, "cell={}", cell);
  cells[cell]->SetConstant(cte);
}

void GPChromosome::ResetConstant(int idx) {
  int cell = idx / (nFunctionsInputs + 1) + nProgramInputs;
  LOG_F(1, "cell={}", cell);
  cells[cell]->ResetConstant();
  if (cells[cell]->Evaluated())
    cells[cell]->ResetConstant();
}
