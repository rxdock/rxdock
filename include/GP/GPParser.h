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

#include "Cell.h"
#include "Command.h"
#include "GPChromosome.h"
#include "GPTypes.h"

namespace rxdock {

namespace geneticprogram {

class GPParser {
public:
  static std::string _CT;
  ///////////////////
  // Constructors
  ///////////////////
  GPParser(int, int, int, int);
  GPParser(const GPParser &);

  ///////////////////
  // Destructor
  //////////////////
  virtual ~GPParser();

  void Clear();
  std::ostream &Print(std::ostream &) const;
  friend std::ostream &operator<<(std::ostream &s, const GPParser &p);
  ReturnTypeList Parse(GPChromosomePtr, ReturnTypeList &);
  std::string PrintParse(std::istream &, GPChromosomePtr, bool, bool);
  bool IsConstantInt(int i) { return (i == CINT); }
  bool IsConstantFloat(int i) { return (i == CFLOAT); }

private:
  /////////////////////
  // Private data
  /////////////////////
  ReturnType Eval(GPChromosomePtr, int);
  ReturnType Parse1Output(GPChromosomePtr, int);
  std::string PrintEval(GPChromosomePtr, int, bool, bool);
  std::string PrintParse1Output(GPChromosomePtr, int, bool, bool);
  // GPChromosomePtr chrom;
  // CellList cells;
  CommandList commands;
  int nProgramInputs, nFunctionsInputs, nProgramOutputs;
  int ADD, SUB, MUL, DIV, IF, CINT, CFLOAT;
};

// Useful typedefs
typedef SmartPtr<GPParser> GPParserPtr; // Smart pointer

} // namespace geneticprogram

} // namespace rxdock

#endif //_GPParser_H_
