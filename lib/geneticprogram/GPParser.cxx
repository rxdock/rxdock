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

#include "rxdock/geneticprogram/GPParser.h"
#include "rxdock/geneticprogram/Cell.h"
#include "rxdock/geneticprogram/Command.h"
#include <fstream>
#include <sstream>

using namespace rxdock;
using namespace rxdock::geneticprogram;

const std::string GPParser::_CT = "GPParser";
int Command::ntabs = 0;
bool Command::inside = false;

// Constructors
GPParser::GPParser(int npi, int nfi, int nn, int no) {
  nProgramInputs = npi;
  nFunctionsInputs = nfi;
  nProgramOutputs = no;
  commands.push_back(new AddCommand());
  commands.push_back(new SubCommand());
  commands.push_back(new MulCommand());
  commands.push_back(new DivCommand());
  commands.push_back(new IfCommand());
  IF = 4;
  commands.push_back(new LogCommand());
  commands.push_back(new ExpCommand());
}

GPParser::GPParser(const GPParser &p) {
  nProgramInputs = p.nProgramInputs;
  nFunctionsInputs = p.nFunctionsInputs;
  nProgramOutputs = p.nProgramOutputs;
  commands.push_back(new AddCommand());
  commands.push_back(new SubCommand());
  commands.push_back(new MulCommand());
  commands.push_back(new DivCommand());
  commands.push_back(new IfCommand());
  IF = 4;
  commands.push_back(new LogCommand());
  commands.push_back(new ExpCommand());
}

///////////////////
// Destructor
//////////////////
GPParser::~GPParser() {}

ReturnType GPParser::Eval(GPChromosomePtr chrom, int n) {
  if (!(chrom->Cells(n)->Evaluated())) {
    ReturnType o = Parse1Output(chrom, n);
    chrom->Cells(n)->SetResult(o);
    return o;
  } else
    return chrom->Cells(n)->GetResult();
}

ReturnType GPParser::Parse1Output(GPChromosomePtr chrom, int output) {
  int startCell = (nFunctionsInputs + 1) * (output - nProgramInputs);
  int ncomm = (*chrom)[startCell + nFunctionsInputs];
  int nargs = commands[ncomm]->GetNArgs();
  ReturnType *args = new ReturnType[nargs];
  // I had this two loops as one before, but because is a
  // recursive function, they MUST be separated, otherwise the command
  // can get the arguments from another cell
  for (int i = 0; i < nargs; i++)
    args[i] = Eval(chrom, (*chrom)[startCell + i]);
  for (int i = 0; i < nargs; i++)
    commands[ncomm]->SetArg(i, args[i]);
  ReturnType f2 = commands[ncomm]->Execute();
  delete[] args;
  return (commands[ncomm]->Execute());
}

ReturnTypeList GPParser::Parse(GPChromosomePtr chrom, ReturnTypeList &inputs) {
  // Give the values of the inputs to the corresponding cells
  for (int i = 0; i < nProgramInputs; i++)
    chrom->Cells(i)->SetResult(*(inputs[i]));
  // for each output, find the node where its output is taken
  // from, and parse it
  int size = chrom->size();
  ReturnTypeList outputs;
  for (int i = 0; i < nProgramOutputs; i++) {
    ReturnType o = Eval(chrom, (*chrom)[size - nProgramOutputs + i]);
    outputs.push_back(new ReturnType(o));
  }
  return outputs;
}

std::string GPParser::PrintEval(GPChromosomePtr chrom, int n, bool numbers,
                                bool ins) {
  if (!(chrom->Cells(n)->Named())) {
    std::string o(PrintParse1Output(chrom, n, numbers, ins));
    std::ostringstream nm;
    if (numbers) {
      nm << "(";
      nm << n << " ";
      nm << o << ") ";
    } else
      nm << o;
    chrom->Cells(n)->SetName(nm.str());
    return nm.str();
  } else
    return chrom->Cells(n)->GetName();
}

std::string GPParser::PrintParse1Output(GPChromosomePtr chrom, int output,
                                        bool numbers, bool ins) {
  int startCell = (nFunctionsInputs + 1) * (output - nProgramInputs);
  int ncomm = (*chrom)[startCell + nFunctionsInputs];
  if (ncomm == IF)
    Command::ntabs++;
  int nargs = commands[ncomm]->GetNArgs();
  std::vector<std::string> fs(nargs);
  for (int i = 0; i < commands[ncomm]->GetNArgs(); i++)
    if ((ncomm == IF) && (i == 0)) {
      //           Command::inside = true;
      fs[i] = PrintEval(chrom, (*chrom)[startCell + i], numbers, true);
      //            Command::inside = false;
    } else
      fs[i] = PrintEval(chrom, (*chrom)[startCell + i], numbers, ins);
  for (int i = 0; i < commands[ncomm]->GetNArgs(); i++)
    commands[ncomm]->SetNameArg(i, fs[i]);
  std::string c;
  Command::inside = ins;
  if (ncomm == IF) {
    c = commands[ncomm]->GetName();
    Command::ntabs--;
  } else
    c = commands[ncomm]->GetName();
  return (c);
}

std::string GPParser::PrintParse(std::istream &istr, GPChromosomePtr chrom,
                                 bool numbers, bool ins) {
  // Give the names of the inputs to the corresponding cells
  std::string name;
  int node;
  while (istr >> node) {
    istr.get();
    getline(istr, name);
    chrom->Cells(node)->SetName(name + " ");
  }
  // for each output, find the node where its output is taken
  // from, and parse it
  int size = chrom->size();
  ReturnTypeList outputs;
  std::string f;
  for (int i = 0; i < nProgramOutputs; i++)
    f += PrintEval(chrom, (*chrom)[size - nProgramOutputs + i], numbers, ins) +
         "\n";
  return f;
}

// std::ostream& Print(std::ostream&) const;
// friend std::ostream& operator<<(std::ostream& s, const GPParser &p);
