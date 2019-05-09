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

// GPParser Class. Parses a genome from a cartesian GP

#ifndef _RBT_GPPARSER_H_
#define _RBT_GPPARSER_H_

#include "RbtCell.h"
#include "RbtCommand.h"
#include "RbtContainers.h"
#include "RbtGPChromosome.h"
#include "RbtGPTypes.h"

class RbtGPParser {
public:
  static std::string _CT;
  ///////////////////
  // Constructors
  ///////////////////
  RbtGPParser(int, int, int, int);
  RbtGPParser(const RbtGPParser &);

  ///////////////////
  // Destructor
  //////////////////
  virtual ~RbtGPParser();

  void Clear();
  std::ostream &Print(std::ostream &) const;
  friend std::ostream &operator<<(std::ostream &s, const RbtGPParser &p);
  RbtReturnTypeList Parse(RbtGPChromosomePtr, RbtReturnTypeList &);
  std::string PrintParse(std::istream &, RbtGPChromosomePtr, bool, bool);
  bool IsConstantInt(int i) { return (i == CINT); }
  bool IsConstantFloat(int i) { return (i == CFLOAT); }

private:
  /////////////////////
  // Private data
  /////////////////////
  RbtReturnType Eval(RbtGPChromosomePtr, int);
  RbtReturnType Parse1Output(RbtGPChromosomePtr, int);
  std::string PrintEval(RbtGPChromosomePtr, int, bool, bool);
  std::string PrintParse1Output(RbtGPChromosomePtr, int, bool, bool);
  // RbtGPChromosomePtr chrom;
  // RbtCellList cells;
  RbtCommandList commands;
  int nProgramInputs, nFunctionsInputs, nProgramOutputs;
  int ADD, SUB, MUL, DIV, IF, CINT, CFLOAT;
};

// Useful typedefs
typedef SmartPtr<RbtGPParser> RbtGPParserPtr; // Smart pointer

#endif //_RbtGPParser_H_
