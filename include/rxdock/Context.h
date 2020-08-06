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

// Stores the information that is global to the interpreter
// in this case, values of ctes, names of variables,...
// It basically contains a map of variables, class Vble is also defined here

#ifndef _RBT_CONTEXT_H_
#define _RBT_CONTEXT_H_

#include "BaseSF.h"
#include "DockingSite.h"
#include "GPTypes.h"
#include "Model.h"
#include "Vble.h"
#include <fstream>

namespace rxdock {

typedef std::map<std::string, VblePtr> StringVbleMap; // Map of Vbles
typedef std::map<int, VblePtr> IntVbleMap;            // Map of Vbles
typedef StringVbleMap::iterator StringVbleMapIter;
typedef IntVbleMap::iterator IntVbleMapIter;

class Context {
public:
  static std::string _CT;
  ///////////////////
  // Constructors
  ///////////////////
  Context(const Context &c);
  Context(); // default constructor disabled
             ///////////////////
             // Destructors
             ///////////////////
  virtual ~Context();
  virtual void Assign(std::string, ReturnType) = 0;
  virtual void Assign(int, ReturnType) = 0;
  virtual const Vble &GetVble(int) = 0;
  virtual const Vble &GetVble(std::string) = 0;
  virtual void SetVble(int key, const Vble &v) = 0;
  //    virtual String GetName(Int)=0;
  //    virtual ReturnType GetValue(Int)=0;
  //    virtual String GetName(String)=0;
  //    virtual ReturnType GetValue(String)=0;
};

class CellContext : public Context {
public:
  // static String _CT;
  CellContext(std::ifstream &ifile);
  CellContext();
  CellContext(const CellContext &c);
  virtual ~CellContext();
  void Assign(int key, ReturnType val) {
    IntVbleMapIter it = vm.find(key);
    if (it != vm.end())
      vm[key]->SetValue(val);
    else {
      vm[key] = new Vble(std::to_string(val), val);
    }
  }
  void Assign(std::string s, ReturnType val) {
    throw Error(_WHERE_, "This is not a string context");
  }
  //    String GetName(Int key) { return vm[key].GetName();}
  //    ReturnType GetValue(Int key) {return vm[key].GetValue();}
  //    String GetName(String){return "";}
  //    ReturnType GetValue(String){return 0.0;}
  const Vble &GetVble(int key) { return *(vm[key]); }
  void SetVble(int key, const Vble &v) { *(vm[key]) = v; }
  const Vble &GetVble(std::string key) {
    throw Error(_WHERE_, "This is not a string context");
  }
  //    void Clear();

private:
  IntVbleMap vm;
  int ninputs;
};

class StringContext : public Context {

public:
  StringContext();
  StringContext(SmartPtr<std::ifstream> ifile);
  StringContext(const StringContext &c);
  virtual ~StringContext();
  void Assign(std::string key, ReturnType val) {
    StringVbleMapIter it = vm.find(key);
    if (it != vm.end())
      vm[key]->SetValue(val);
    else {
      vm[key] = new Vble(key, val);
    }
  }
  void Assign(int i, ReturnType val) {
    throw Error(_WHERE_, "This is not a cell context");
  }

  //    String GetName(Int){return "";}
  //    ReturnType GetValue(Int){return 0.0;}
  //    String GetName(String key) { return vm[key].GetName();}
  //    ReturnType GetValue(String key) {return vm[key].GetValue();}
  double Get(ModelPtr lig, std::string name);
  double Get(ModelPtr rec, DockingSitePtr site, std::string name);
  double Get(BaseSF *spSF, std::string name, ModelPtr lig);
  const Vble &GetVble(std::string key) { return *(vm[key]); }
  const Vble &GetVble(int key) {
    throw Error(_WHERE_, "This is not a cell context");
  }
  void SetVble(int key, const Vble &v) { *(vm[""]) = v; }
  void UpdateLigs(ModelPtr lig);
  void UpdateSite(ModelPtr rec, DockingSitePtr site);
  void UpdateScores(BaseSF *spSF, ModelPtr lig);

private:
  StringVbleMap vm;
};

// Useful typedefs
typedef SmartPtr<CellContext> CellContextPtr;        // Smart pointer
typedef std::vector<CellContextPtr> CellContextList; // Vector of smart
                                                     // pointers
typedef CellContextList::iterator CellContextListIter;
typedef CellContextList::const_iterator CellContextListConstIter;
typedef SmartPtr<StringContext> StringContextPtr; // Smart pointer
typedef std::vector<StringContextPtr>
    StringContextList; // Vector of smart pointers
typedef StringContextList::iterator StringContextListIter;
typedef StringContextList::const_iterator StringContextListConstIter;
typedef SmartPtr<Context> ContextPtr;        // Smart pointer
typedef std::vector<ContextPtr> ContextList; // Vector of smart pointers
typedef ContextList::iterator ContextListIter;
typedef ContextList::const_iterator ContextListConstIter;

} // namespace rxdock

#endif //_Context_H_
