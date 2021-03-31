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

// Abstract StringTokenIter class. Defines the interface for a token iterator
// to be used by the parser.

#ifndef _RBT_STRINGTOKENITER_H_
#define _RBT_STRINGTOKENITER_H_

#include "rxdock/Config.h"
#include "rxdock/Context.h"
#include "rxdock/TokenIter.h"
#include "rxdock/geneticprogram/GPChromosome.h"
#include <fstream>
#include <stack>

namespace rxdock {

typedef SmartPtr<std::istream> istreamPtr;
class StringTokenIter : public TokenIter {
public:
  static const std::string _CT;
  ///////////////////
  // Constructors
  ///////////////////
  StringTokenIter(const StringTokenIter &);
  StringTokenIter(const istreamPtr, ContextPtr);

  ///////////////////
  // Destructor
  //////////////////
  virtual ~StringTokenIter();

  void copy(const StringTokenIter &);
  // void First();
  void Next(ContextPtr);
  // Bool IsDone() const;
  TokenPtr Current();

  //  Bool operator==(const StringTokenIter&) const;
  //  std::ostream& Print(std::ostream&) const;
  //  friend std::ostream& operator<<(std::ostream& s, const StringTokenIter
  //  &p);

  ////////////////////
  // Private methods
  ////////////////////

private:
  StringTokenIter(); // Default constructor disabled
  TokenPtr translate(std::string);

  istreamPtr filep;
  TokenPtr current;
  std::string strtok;
  ContextPtr contextp;
};

// Useful typedefs
typedef SmartPtr<StringTokenIter> StringTokenIterPtr; // Smart pointer
typedef std::vector<StringTokenIterPtr>
    StringTokenIterList; // Vector of smart pointers
typedef StringTokenIterList::iterator StringTokenIterListIter;
typedef StringTokenIterList::const_iterator StringTokenIterListConstIter;

} // namespace rxdock

#endif //_StringTokenIter_H_
