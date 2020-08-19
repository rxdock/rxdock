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

// Target intramolecular dihedral scoring function
// for flexible receptors
#ifndef _RBTDIHEDRALTARGETSF_H_
#define _RBTDIHEDRALTARGETSF_H_

#include "rxdock/BaseInterSF.h"
#include "rxdock/DihedralSF.h"

namespace rxdock {

class DihedralTargetSF : public BaseInterSF, public DihedralSF {
public:
  // Class type string
  static std::string _CT;
  // Parameter names

  DihedralTargetSF(const std::string &strName = "dihedral");
  virtual ~DihedralTargetSF();

protected:
  virtual void SetupReceptor();
  virtual void SetupLigand();
  virtual void SetupScore();
  virtual double RawScore() const;

  // Clear the dihedral list
  // As we are not using smart pointers, there is some memory management to do
  void ClearReceptor();

private:
  DihedralList m_dihList;
};

} // namespace rxdock

#endif //_RBTDIHEDRALTARGETSF_H_
