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

#include "RbtGPChromosome.h"
#include "RbtCell.h"
#include "RbtCommand.h"

#include <fstream>

std::string RbtGPChromosome::_CT("RbtGPChromosome");

// Constructors
RbtGPChromosome::RbtGPChromosome(int npi, int nfi, int nn, int no, int nr,
                                 int nc) {
  chrom = RbtIntList(nn * (nfi + 1) + no);
  int i = 0, n;
  nProgramInputs = npi;
  nFunctionsInputs = nfi;
  int ncells = nn + npi; // number of nodes + inputs
  nProgramOutputs = no;
  nColumns = nc;
  nRows = nr;
  cells = RbtCellList(ncells);
  for (int i = 0; i < ncells; i++)
    cells[i] = new RbtCell();
}
RbtGPChromosome::RbtGPChromosome(const RbtGPChromosome &c) {
  chrom = RbtIntList(c.chrom.size());
  cells = RbtCellList(c.cells.size());
  for (int i = 0; i < cells.size(); i++)
    cells[i] = c.cells[i];
  for (int i = 0; i < chrom.size(); i++)
    chrom[i] = c.chrom[i];
  nProgramInputs = c.nProgramInputs;
  nFunctionsInputs = c.nFunctionsInputs;
  nProgramOutputs = c.nProgramOutputs;
}

///////////////////
// Destructor
//////////////////
RbtGPChromosome::~RbtGPChromosome() {}

void RbtGPChromosome::Clear() {
  for (int i = 0; i < cells.size(); i++)
    cells[i]->Clear();
  RbtCommand::Clear();
}
RbtGPChromosome &RbtGPChromosome::operator=(const RbtGPChromosome &c) {
  for (int i = 0; i < cells.size(); i++)
    *(cells[i]) = *(c.cells[i]);
  for (int i = 0; i < chrom.size(); i++)
    chrom[i] = c.chrom[i];
  nProgramInputs = c.nProgramInputs;
  nFunctionsInputs = c.nFunctionsInputs;
  nProgramOutputs = c.nProgramOutputs;
}

std::istream &operator>>(std::istream &s, RbtGPChromosome &p) {
  for (int i = 0; i < p.chrom.size(); i++)
    s >> p.chrom[i];
  int cell;
  RbtReturnType value;
  while (!s.eof()) {
    s >> cell >> value;
    p.Cells(cell)->SetConstant(value);
  }
}

std::ostream &RbtGPChromosome::Print(std::ostream &s) const {
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
  for (int i = 0; i < cells.size(); i++)
    if (cells[i]->Constant())
      s << i << "\t" << cells[i]->GetResult() << std::endl;
  return s;
}

std::ostream &operator<<(std::ostream &s, const RbtGPChromosome &p) {
  p.Print(s);
  return s;
}

void RbtGPChromosome::SetConstant(RbtReturnType cte, int idx) {
  int cell = idx / (nFunctionsInputs + 1) + nProgramInputs;
  //    std::cout << cell << std::endl;
  cells[cell]->SetConstant(cte);
}

void RbtGPChromosome::ResetConstant(int idx) {
  int cell = idx / (nFunctionsInputs + 1) + nProgramInputs;
  //    std::cout << cell << std::endl;
  cells[cell]->ResetConstant();
  if (cells[cell]->Evaluated())
    cells[cell]->ResetConstant();
}
