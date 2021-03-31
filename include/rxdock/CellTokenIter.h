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

// Defines the interface for a token iterator to be used by the parser

#ifndef _RBT_CELLTOKENITER_H_
#define _RBT_CELLTOKENITER_H_

#include "rxdock/Config.h"
#include "rxdock/Context.h"
#include "rxdock/Token.h"
#include "rxdock/geneticprogram/GPChromosome.h"
#include <stack>

#include "rxdock/TokenIter.h"

namespace rxdock {

class CellTokenIter : public TokenIter {
public:
  static const std::string _CT;
  ///////////////////
  // Constructors
  ///////////////////
  CellTokenIter(const CellTokenIter &);
  CellTokenIter(const geneticprogram::GPChromosomePtr, ContextPtr);

  ///////////////////
  // Destructor
  //////////////////
  virtual ~CellTokenIter();

  virtual void copy(const CellTokenIter &);
  // void First();
  void Next(ContextPtr);
  // Bool IsDone() const;
  TokenPtr Current();

  ////////////////////
  // Private methods
  ////////////////////

private:
  CellTokenIter(); // Default constructor disabled
  geneticprogram::GPChromosomePtr chrom;
  TokenPtr current;
  std::stack<int> cells;
  ContextPtr contextp;
};

// Useful typedefs
typedef SmartPtr<CellTokenIter> CellTokenIterPtr; // Smart pointer
typedef std::vector<CellTokenIterPtr>
    CellTokenIterList; // Vector of smart pointers
typedef CellTokenIterList::iterator CellTokenIterListIter;
typedef CellTokenIterList::const_iterator CellTokenIterListConstIter;

} // namespace rxdock

#endif //_CellTokenIter_H_
