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

#include "Token.h"
#include "Commands.h"
#include "Debug.h"

using namespace rxdock;

static const Vble default_vble;

std::string Token::_CT("Token");

///////////////////
// Constructors
///////////////////
Token::Token(const Vble &v) : comm(-1), vble(v), isvble(true) {
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

Token::Token(Commands c) : comm(c), vble(default_vble), isvble(false) {
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

Token::Token(const Token &t) : comm(t.comm), vble(t.vble), isvble(t.isvble) {
  _RBTOBJECTCOUNTER_COPYCONSTR_(_CT);
}

///////////////////
// Destructor
//////////////////
Token::~Token() { _RBTOBJECTCOUNTER_DESTR_(_CT); }

/*void Token::copy(const Token &t)
{
    isvble = t.isvble;
    comm = t.comm;
    vble = t.vble;
}*/

const Vble &Token::GetVble() const {
  if (!isvble)
    throw Error(_WHERE_, "the token is not a vble");
  return vble;
}
/* void Token::SetVbleNumber(Int vn)
{
    vble.SetValue(vn);
}*/

bool Token::IsVble() { return (isvble); }
bool Token::IsLog() { return (!isvble && comm.IsLog()); }
bool Token::IsExp() { return (!isvble && comm.IsExp()); }
bool Token::IsAdd() { return (!isvble && comm.IsAdd()); }
bool Token::IsSub() { return (!isvble && comm.IsSub()); }
bool Token::IsMul() { return (!isvble && comm.IsMul()); }
bool Token::IsDiv() { return (!isvble && comm.IsDiv()); }
bool Token::IsAnd() { return (!isvble && comm.IsAnd()); }
bool Token::IsIf() { return (!isvble && comm.IsIf()); }
