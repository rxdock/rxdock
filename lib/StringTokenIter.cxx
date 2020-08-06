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

#include "StringTokenIter.h"
#include "Commands.h"
#include <cerrno>
#include <cstdlib>

using namespace rxdock;

std::string StringTokenIter::_CT("StringTokenIter");

///////////////////
// Constructors
///////////////////
StringTokenIter::StringTokenIter(const istreamPtr fn, ContextPtr co)
    : filep(fn), contextp(co) {
  (std::istream &)(*filep) >> strtok;
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

StringTokenIter::StringTokenIter(const StringTokenIter &ti)
    : filep(ti.filep), current(ti.current), strtok(ti.strtok),
      contextp(ti.contextp) {
  _RBTOBJECTCOUNTER_COPYCONSTR_(_CT);
}

StringTokenIter::~StringTokenIter() { _RBTOBJECTCOUNTER_DESTR_(_CT); }

void StringTokenIter::Next(ContextPtr) { (*filep) >> strtok; }

TokenPtr StringTokenIter::Current() {
  current = translate(strtok);
  return current;
}

void StringTokenIter::copy(const StringTokenIter &ti) {
  strtok = ti.strtok;
  current = ti.current;
  filep = ti.filep;
  contextp = ti.contextp;
}

TokenPtr StringTokenIter::translate(std::string s) {
  if (s.length() == 0)
    throw Error(_WHERE_, "Missing token, can't translate the expression");
  if (s == "+")
    return new Token(Commands(Commands::ADD));
  if (s == "-")
    return new Token(Commands(Commands::SUB));
  if (s == "*")
    return new Token(Commands(Commands::MUL));
  if (s == "/")
    return new Token(Commands(Commands::DIV));
  if (s == "if") {
    return new Token(Commands(Commands::IF));
  }
  if (s == "log")
    return new Token(Commands(Commands::LOG));
  if (s == "exp")
    return new Token(Commands(Commands::EXP));
  if (s == "and")
    return new Token(Commands(Commands::AND));
  if ((std::string::npos != s.find("SCORE")) ||
      (std::string::npos != s.find("SITE")) ||
      (std::string::npos != s.find("LIG"))) {
    contextp->Assign(s, 0.0);
    return new Token(contextp->GetVble(s));
  } else // I assume is a double. If it is not, need to send an Error
  {
    char *error;
    errno = 0;
    double val = std::strtod(s.c_str(), &error);
    if (!errno && !*error) // This checks for errors
    {
      contextp->Assign(s, val);
      return new Token(contextp->GetVble(s));
    }
  }
  throw Error(_WHERE_, "Can't read " + s);
}
