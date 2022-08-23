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

#ifndef _RBTSETUPSASF_H_
#define _RBTSETUPSASF_H_

#include "rxdock/BaseInterSF.h"

namespace rxdock {

// after Hasel, Hendrickson and Still
struct SimpleSolvationParameters {
  double p;
  double r;
};

class SetupSASF : public BaseInterSF {
  AtomList theLigandList;   // ligand typing
  AtomList theReceptorList; // receptor typing

public:
  static const std::string _CT;

  SetupSASF(const std::string &strName = "setup-sa");
  ~SetupSASF();

protected:
  virtual void SetupReceptor();
  virtual void SetupLigand();
  virtual void SetupScore();
  virtual double RawScore() const;

  void SetupReceptorSATypes(void);
  void SetupLigandSATypes(void);

private:
  friend void to_json(json &j, const SetupSASF &sasf) {
    json atomList;
    for (const auto &cIter : sasf.theLigandList) {
      json atom = *cIter;
      atomList.push_back(atom);
    }
    json atomList2;
    for (const auto &cIter : sasf.theReceptorList) {
      json atom = *cIter;
      atomList2.push_back(atom);
    }

    j = json{{"ligand-list", atomList}, {"receptor-list", atomList2}};
  };

  friend void from_json(const json &j, SetupSASF &sasf) {
    sasf.theLigandList.clear();
    sasf.theLigandList.reserve(j.at("atoms").size());
    for (auto &atom : j.at("atoms")) {
      AtomPtr spAtom = AtomPtr(new Atom(atom));
      sasf.theLigandList.push_back(spAtom);
    }

    sasf.theReceptorList.clear();
    sasf.theReceptorList.reserve(j.at("atoms").size());
    for (auto &atom : j.at("atoms")) {
      AtomPtr spAtom = AtomPtr(new Atom(atom));
      sasf.theReceptorList.push_back(spAtom);
    }
  };
};

} // namespace rxdock

#endif //_RBTSETUPSASF_H_
