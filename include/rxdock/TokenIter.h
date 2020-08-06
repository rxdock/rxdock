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

// Abstract TokenIter class. Defines the interface for a token iterator
// to be used by the parser.

#ifndef _RBT_TOKENITER_H_
#define _RBT_TOKENITER_H_

#include "Context.h"
#include "Token.h"

namespace rxdock {

class TokenIter {
public:
  //  virtual void copy(const TokenIter &) = 0;
  virtual void Next(ContextPtr) = 0;
  virtual TokenPtr Current() = 0;
  virtual ~TokenIter() {}

  //  virtual TokenPtr IsDone() const = 0;
  //  virtual void First() = 0;
  //  Bool operator==(const TokenIter&) const;
  //  std::ostream& Print(std::ostream&) const;
  //  friend std::ostream& operator<<(std::ostream& s, const TokenIter &p);
};

// Useful typedefs
typedef SmartPtr<TokenIter> TokenIterPtr;        // Smart pointer
typedef std::vector<TokenIterPtr> TokenIterList; // Vector of smart pointers
typedef TokenIterList::iterator TokenIterListIter;
typedef TokenIterList::const_iterator TokenIterListConstIter;

} // namespace rxdock

#endif //_TokenIter_H_
