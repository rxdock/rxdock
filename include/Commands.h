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

// Defines the interface for a token to be used by the parser

#ifndef _RBT_COMMANDS_H_
#define _RBT_COMMANDS_H_

namespace rxdock {

class Commands {
public:
  enum CommName {
    ADD = 0,
    SUB = 1,
    MUL = 2,
    DIV = 3,
    IF = 4,
    LOG = 5,
    EXP = 6,
    RCTE = 7,
    AND = 8,
    NULLC = -1
  };

  ///////////////////
  // Constructors
  ///////////////////
  Commands(const Commands &c) : name(c.name), nargs(c.nargs) {}
  inline Commands(int ncomm) {
    name = CommName(ncomm);
    switch (name) {
    case ADD:
    case SUB:
    case MUL:
    case DIV:
    case AND:
      nargs = 2;
      break;
    case IF:
      nargs = 3;
      break;
    case LOG:
    case EXP:
      nargs = 1;
      break;
    case RCTE:
      nargs = 0;
      break;
    case NULLC:
      nargs = 0;
      break;
    }
  }

  inline CommName GetName() const { return name; }

  inline int GetNArgs() const { return nargs; }

  inline bool IsAdd() const { return (name == ADD); }
  inline bool IsSub() const { return (name == SUB); }
  inline bool IsMul() const { return (name == MUL); }
  inline bool IsDiv() const { return (name == DIV); }
  inline bool IsIf() const { return (name == IF); }
  inline bool IsLog() const { return (name == LOG); }
  inline bool IsExp() const { return (name == EXP); }
  inline bool IsRCte() const { return (name == RCTE); }
  inline bool IsAnd() const { return (name == AND); }

  ///////////////////
  // Destructor
  //////////////////
  virtual ~Commands() {}

  virtual void copy(const Commands &c) {
    name = c.name;
    nargs = c.nargs;
  }

  ////////////////////
  // Private methods
  ////////////////////

private:
  Commands(); // Default constructor disabled
  CommName name;
  int nargs;
};

// Useful typedefs
// typedef SmartPtr<Commands> CommandsPtr;  //Smart pointer
// typedef std::vector<CommandsPtr> CommandsList;//Vector of smart
// pointers typedef CommandsList::iterator CommandsListIter; typedef
// CommandsList::const_iterator CommandsListConstIter;

// Compare class to compare different genomes inside a population
// This is used to sort the genomes in a population depending
// of the value of their scoring function.

} // namespace rxdock

#endif //_Commands_H_
