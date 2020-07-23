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

// Simple scoring function which estimates the loss of rotational entropy of the
// ligand from the number of rotable bonds.

#ifndef _RBTROTSF_H_
#define _RBTROTSF_H_

#include "RbtBaseInterSF.h"

namespace rxdock {

class RbtRotSF : public RbtBaseInterSF {
public:
  // Static data member for class type
  static std::string _CT;
  // Parameter names
  // Boolean controlling whether to include bonds to NH3+ in rotable bond count
  static std::string _INCNH3;
  // Boolean controlling whether to include bonds to OH in  rotable bond count
  static std::string _INCOH;

  RbtRotSF(const std::string &strName = "ROT");
  virtual ~RbtRotSF();

protected:
  virtual void SetupReceptor();
  virtual void SetupLigand();
  virtual void SetupScore();
  virtual double RawScore() const;
  void ParameterUpdated(const std::string &strName);

private:
  int nRot;
  bool bIncNH3;
  bool bIncOH;
};

} // namespace rxdock

#endif //_RBTROTSF_H_
