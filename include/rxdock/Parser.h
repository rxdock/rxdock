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

// Parser class. Takes a list of tokens and parsers them to create a
// Filter Expression

#ifndef _RBT_PARSER_H_
#define _RBT_PARSER_H_

#include "rxdock/Context.h"
#include "rxdock/FilterExpression.h"
#include "rxdock/Token.h"
#include "rxdock/TokenIter.h"

namespace rxdock {

class FilterExpression;

class Parser {
public:
  static const std::string _CT;
  ///////////////////
  // Constructors
  ///////////////////
  Parser(const Parser &);
  Parser();
  FilterExpressionPtr Parse(TokenIterPtr, ContextPtr);

  ///////////////////
  // Destructor
  //////////////////
  virtual ~Parser();

  ////////////////////
  // Private methods
  ////////////////////

private:
};

// Useful typedefs
typedef SmartPtr<Parser> ParserPtr;        // Smart pointer
typedef std::vector<ParserPtr> ParserList; // Vector of smart pointers
typedef ParserList::iterator ParserListIter;
typedef ParserList::const_iterator ParserListConstIter;

// Compare class to compare different genomes inside a population
// This is used to sort the genomes in a population depending
// of the value of their scoring function.

} // namespace rxdock

#endif //_Parser_H_
