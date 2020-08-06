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

// Command Class. Command represents a function or input variable of the CGP

#ifndef _RBT_COMMAND_H_
#define _RBT_COMMAND_H_

#include "rxdock/geneticprogram/GPTypes.h"

#include <cmath>

namespace rxdock {

class Command {
public:
  virtual ReturnType Execute() = 0;
  virtual int GetNArgs() = 0;
  virtual std::string GetName() = 0;
  virtual void SetArg(int i, ReturnType f) = 0;       // { *(Arg[i]) = f;};
  virtual void SetNameArg(int i, std::string &n) = 0; // { NameArg[i] = n;};
  static void Clear() { ntabs = 0; }

  ///////////////////
  // Destructor
  //////////////////
  virtual ~Command() {}
  static int ntabs;
  static bool inside;

protected:
  ///////////////////
  // Constructors
  ///////////////////
  Command() {} // Default constructor disabled
  std::string tabs(int n) {
    std::string st = "";
    for (int i = 0; i < n; i++)
      st += "  ";
    return st;
  }

private:
  //    friend std::ostream& operator<<(std::ostream& s, const Command &p);
};

// Useful typedefs
typedef SmartPtr<Command> CommandPtr;        // Smart pointer
typedef std::vector<CommandPtr> CommandList; // Vector of smart pointers
typedef CommandList::iterator CommandListIter;
typedef CommandList::const_iterator CommandListConstIter;
class AddCommand : public Command {
public:
  int GetNArgs() { return 2; }
  std::string GetName() { return ("(" + NameArg[0] + "+ " + NameArg[1] + ")"); }
  AddCommand() {} //{Arg = ReturnType(2);}
  void SetArg(int i, ReturnType f) { Arg[i] = f; }
  void SetNameArg(int i, std::string &n) { NameArg[i] = n; }
  ReturnType Execute() { return Arg[0] + Arg[1]; }

private:
  ReturnType Arg[2];
  std::string NameArg[2];
};

class SubCommand : public Command {
public:
  int GetNArgs() { return 2; }
  std::string GetName() { return ("(" + NameArg[0] + "- " + NameArg[1] + ")"); }
  void SetArg(int i, ReturnType f) { Arg[i] = f; }
  void SetNameArg(int i, std::string &n) { NameArg[i] = n; }
  SubCommand() {}
  ReturnType Execute() { return Arg[0] - Arg[1]; }

private:
  ReturnType Arg[2];
  std::string NameArg[2];
};

class MulCommand : public Command {
public:
  int GetNArgs() { return 2; }
  std::string GetName() { return ("(" + NameArg[0] + "* " + NameArg[1] + ")"); }
  void SetArg(int i, ReturnType f) { Arg[i] = f; }
  void SetNameArg(int i, std::string &n) { NameArg[i] = n; }
  MulCommand() {}
  ReturnType Execute() { return Arg[0] * Arg[1]; }

private:
  ReturnType Arg[2];
  std::string NameArg[2];
};

class DivCommand : public Command {
public:
  int GetNArgs() { return 2; }
  std::string GetName() {
    return ("(" + NameArg[0] + "div " + NameArg[1] + ")");
  }
  DivCommand() {}
  void SetArg(int i, ReturnType f) { Arg[i] = f; }
  void SetNameArg(int i, std::string &n) { NameArg[i] = n; }
  ReturnType Execute() {
    if (std::fabs(Arg[1]) < 0.000001)
      return Arg[0];
    return Arg[0] / Arg[1];
  }

private:
  ReturnType Arg[2];
  std::string NameArg[2];
};

class IfCommand : public Command {
public:
  int GetNArgs() { return 3; }
  std::string GetName() {
    //	return ("\niff(" + NameArg[0] + "," + NameArg[1] + "," +
    //			NameArg[2] + ")");};
    if (inside)
      return ("iff(" + NameArg[0] + " , " + NameArg[1] + " , " + NameArg[2] +
              ")");
    return ("\n" + tabs(ntabs - 1) + "if " + NameArg[0] + "> 0 then\n" +
            tabs(ntabs) + NameArg[1] + "\n" + tabs(ntabs - 1) + "else \n" +
            tabs(ntabs) + NameArg[2] + "\n" + tabs(ntabs - 1) + "end\n");
  }
  IfCommand() {}
  void SetArg(int i, ReturnType f) { Arg[i] = f; }
  void SetNameArg(int i, std::string &n) { NameArg[i] = n; }
  ReturnType Execute() {
    if (Arg[0] > 0.0)
      return Arg[1];
    return Arg[2];
  }

private:
  ReturnType Arg[3];
  std::string NameArg[3];
};

class LogCommand : public Command {
public:
  int GetNArgs() { return 1; }
  std::string GetName() { return ("log(" + NameArg[0] + ")"); }
  LogCommand() {}
  void SetArg(int i, ReturnType f) { Arg[i] = f; }
  void SetNameArg(int i, std::string &n) { NameArg[i] = n; }
  ReturnType Execute() {
    if (std::fabs(Arg[0]) < 0.000001)
      return 0;
    return (std::log(std::fabs(Arg[0])));
  }

private:
  ReturnType Arg[1];
  std::string NameArg[1];
};

class ExpCommand : public Command {
public:
  int GetNArgs() { return 1; }
  std::string GetName() { return ("exp(" + NameArg[0] + ")"); }
  ExpCommand() {}
  void SetArg(int i, ReturnType f) { Arg[i] = f; }
  void SetNameArg(int i, std::string &n) { NameArg[i] = n; }
  ReturnType Execute() {
    if (Arg[0] > 200)
      return std::exp(200);
    if (Arg[0] < -200)
      return 0;
    return std::exp(Arg[0]);
  }

private:
  ReturnType Arg[2];
  std::string NameArg[2];
};

} // namespace rxdock

#endif //_Command_H_
