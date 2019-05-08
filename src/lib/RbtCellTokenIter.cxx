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

#include "RbtCellTokenIter.h"

std::string RbtCellTokenIter::_CT("RbtCellTokenIter");

///////////////////
// Constructors
///////////////////
RbtCellTokenIter::RbtCellTokenIter(const RbtGPChromosomePtr c,
                                   RbtContextPtr contextp)
    : chrom(c) {
  int startCell = chrom->GetStartingCell();
  if (chrom->IsProgramInput(startCell)) {
    current = new RbtToken(contextp->GetVble(startCell));
  } else {
    RbtCommands comm = chrom->GetCommand(startCell);
    if (comm.IsRCte()) {
      double val = chrom->GetRCte(startCell);
      contextp->Assign(startCell, val);
      current = new RbtToken(contextp->GetVble(startCell));
    } else {
      current = new RbtToken(comm);
      int nargs = comm.GetNArgs();
      for (int i = nargs - 1; i >= 0; i--)
        cells.push(chrom->GetArgument(startCell, i));
    }
  }
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtCellTokenIter::RbtCellTokenIter(const RbtCellTokenIter &ti)
    : chrom(ti.chrom), cells(ti.cells), current(ti.current) {
  _RBTOBJECTCOUNTER_COPYCONSTR_(_CT);
}

RbtCellTokenIter::~RbtCellTokenIter() { _RBTOBJECTCOUNTER_DESTR_(_CT); }

void RbtCellTokenIter::Next(RbtContextPtr contextp) {
  if (cells.size() > 0) {
    int cell = cells.top();
    cells.pop();
    if (chrom->IsProgramInput(cell)) {
      current = new RbtToken(contextp->GetVble(cell));
    } else {
      RbtCommands comm = chrom->GetCommand(cell);
      if (comm.IsRCte()) {
        double val = chrom->GetRCte(cell);
        contextp->Assign(cell, val);
        current = new RbtToken(contextp->GetVble(cell));
      } else {
        current = new RbtToken(comm);
        int nargs = comm.GetNArgs();
        for (int i = nargs - 1; i >= 0; i--)
          cells.push(chrom->GetArgument(cell, i));
      }
    }
  }
}

RbtTokenPtr RbtCellTokenIter::Current() { return current; }

void RbtCellTokenIter::copy(const RbtCellTokenIter &ti) {
  chrom = ti.chrom;
  current = ti.current;
  cells = ti.cells;
}
