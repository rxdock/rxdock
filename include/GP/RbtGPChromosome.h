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

// GP chromosome. Keeps the list of integers that represent a chromosome and
// a list of cells that represent a program.

#ifndef _RBT_GPCHROMOSOME_H
#define _RBT_GPCHROMOSOME_H

#include "RbtCell.h"
#include "RbtCommands.h"
#include "RbtGPTypes.h"

class RbtGPChromosome {
public:
  static std::string _CT;
  ///////////////////
  // Constructors
  ///////////////////
  RbtGPChromosome(int, int, int, int, int, int);
  RbtGPChromosome(const RbtGPChromosome &);

  ///////////////////
  // Destructor
  //////////////////
  virtual ~RbtGPChromosome();

  void Clear();
  std::ostream &Print(std::ostream &) const;
  friend std::ostream &operator<<(std::ostream &s, const RbtGPChromosome &p);
  friend std::istream &operator>>(std::istream &s, RbtGPChromosome &p);
  const int &operator[](int idx) const { return chrom[idx]; }
  int &operator[](int idx) { return chrom[idx]; }
  RbtGPChromosome &operator=(const RbtGPChromosome &c);
  int size() const { return chrom.size(); }
  RbtCellPtr Cells(int idx) const { return (cells[idx]); }
  void SetConstant(RbtReturnType, int);
  void ResetConstant(int);
  int GetStartingCell() { return (chrom[chrom.size() - 1]); }
  RbtReturnType GetRCte(int ncell) { return cells[ncell]->GetResult(); }
  RbtCommands GetCommand(int ncell) {
    if (ncell < nProgramInputs)
      throw RbtError(_WHERE_, "This cell does not have a command");
    RbtCommands comm(chrom[(nFunctionsInputs + 1) * (ncell - nProgramInputs) +
                           nFunctionsInputs]);
    return comm;
  }
  int GetArgument(int ncell, int narg) {
    if (ncell < nProgramInputs)
      throw RbtError(_WHERE_, "This cell does not have arguments");
    return (chrom[(nFunctionsInputs + 1) * (ncell - nProgramInputs) + narg]);
  }
  bool IsProgramInput(int ncell) { return (ncell < nProgramInputs); }

private:
  /////////////////////
  // Private data
  /////////////////////
  std::vector<int> chrom;
  RbtCellList cells;
  int nProgramInputs, nFunctionsInputs, nProgramOutputs, nRows, nColumns;
};

// Useful typedefs
typedef SmartPtr<RbtGPChromosome> RbtGPChromosomePtr; // Smart pointer

#endif //_RBTGPCHROMOSOME_H_
