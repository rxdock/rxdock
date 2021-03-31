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

// Token class. Defines a token of the grammar to be used by the parser.

#ifndef _RBT_TOKEN_H_
#define _RBT_TOKEN_H_

#include "rxdock/Commands.h"
#include "rxdock/Config.h"
#include "rxdock/Vble.h"

namespace rxdock {

class Token {
public:
  static const std::string _CT;
  ///////////////////
  // Constructors
  ///////////////////

  Token(const Vble &);
  Token(Commands);
  Token(const Token &);

  ///////////////////
  // Destructor
  //////////////////
  virtual ~Token();
  const Vble &GetVble() const;
  //  void SetVbleNumber(Int);

  //  void copy(const Token &);
  bool IsVble();
  bool IsLog();
  bool IsExp();
  bool IsAdd();
  bool IsSub();
  bool IsMul();
  bool IsDiv();
  bool IsAnd();
  bool IsIf();

  //  Bool operator==(const Token&) const;
  //  std::ostream& Print(std::ostream&) const;
  //  friend std::ostream& operator<<(std::ostream& s, const Token &p);

  ////////////////////
  // Private methods
  ////////////////////

private:
  Token(); // Default constructor disabled
  Commands comm;
  const Vble &vble;
  bool isvble;
};

// Useful typedefs
typedef SmartPtr<Token> TokenPtr;        // Smart pointer
typedef std::vector<TokenPtr> TokenList; // Vector of smart pointers
typedef TokenList::iterator TokenListIter;
typedef TokenList::const_iterator TokenListConstIter;

// Compare class to compare different genomes inside a population
// This is used to sort the genomes in a population depending
// of the value of their scoring function.

} // namespace rxdock

#endif //_Token_H_
