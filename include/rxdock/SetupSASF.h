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
  static std::string _CT;

  SetupSASF(const std::string &strName = "SETUP_SA");
  ~SetupSASF();

protected:
  virtual void SetupReceptor();
  virtual void SetupLigand();
  virtual void SetupScore();
  virtual double RawScore() const;

  void SetupReceptorSATypes(void);
  void SetupLigandSATypes(void);
};

} // namespace rxdock

#endif //_RBTSETUPSASF_H_
