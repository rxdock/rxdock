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

// Ligand intramolecular dihedral scoring function
#ifndef _RBTDIHEDRALINTRASF_H_
#define _RBTDIHEDRALINTRASF_H_

#include "rxdock/BaseIntraSF.h"
#include "rxdock/DihedralSF.h"

namespace rxdock {

class DihedralIntraSF : public BaseIntraSF, public DihedralSF {
public:
  // Class type string
  static std::string _CT;
  // Parameter names

  DihedralIntraSF(const std::string &strName = "dihedral");
  virtual ~DihedralIntraSF();

protected:
  virtual void SetupScore();
  virtual double RawScore() const;

  // Clear the dihedral list
  // As we are not using smart pointers, there is some memory management to do
  void ClearModel();

private:
  DihedralList m_dihList;
};

} // namespace rxdock

#endif //_RBTDIHEDRALINTRASF_H_
