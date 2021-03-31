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

#include "rxdock/CellTokenIter.h"

using namespace rxdock;

const std::string CellTokenIter::_CT = "CellTokenIter";

///////////////////
// Constructors
///////////////////
CellTokenIter::CellTokenIter(const geneticprogram::GPChromosomePtr c,
                             ContextPtr contextp)
    : chrom(c) {
  int startCell = chrom->GetStartingCell();
  if (chrom->IsProgramInput(startCell)) {
    current = new Token(contextp->GetVble(startCell));
  } else {
    Commands comm = chrom->GetCommand(startCell);
    if (comm.IsRCte()) {
      double val = chrom->GetRCte(startCell);
      contextp->Assign(startCell, val);
      current = new Token(contextp->GetVble(startCell));
    } else {
      current = new Token(comm);
      int nargs = comm.GetNArgs();
      for (int i = nargs - 1; i >= 0; i--)
        cells.push(chrom->GetArgument(startCell, i));
    }
  }
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

CellTokenIter::CellTokenIter(const CellTokenIter &ti)
    : chrom(ti.chrom), current(ti.current), cells(ti.cells) {
  _RBTOBJECTCOUNTER_COPYCONSTR_(_CT);
}

CellTokenIter::~CellTokenIter() { _RBTOBJECTCOUNTER_DESTR_(_CT); }

void CellTokenIter::Next(ContextPtr contextp) {
  if (cells.size() > 0) {
    int cell = cells.top();
    cells.pop();
    if (chrom->IsProgramInput(cell)) {
      current = new Token(contextp->GetVble(cell));
    } else {
      Commands comm = chrom->GetCommand(cell);
      if (comm.IsRCte()) {
        double val = chrom->GetRCte(cell);
        contextp->Assign(cell, val);
        current = new Token(contextp->GetVble(cell));
      } else {
        current = new Token(comm);
        int nargs = comm.GetNArgs();
        for (int i = nargs - 1; i >= 0; i--)
          cells.push(chrom->GetArgument(cell, i));
      }
    }
  }
}

TokenPtr CellTokenIter::Current() { return current; }

void CellTokenIter::copy(const CellTokenIter &ti) {
  chrom = ti.chrom;
  current = ti.current;
  cells = ti.cells;
}
