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

// Generic nmr restraint scoring function

#ifndef _RBTNMRSF_H_
#define _RBTNMRSF_H_

#include "rxdock/BaseIdxSF.h"
#include "rxdock/BaseInterSF.h"
#include "rxdock/NoeRestraint.h"

namespace rxdock {

class NmrSF : public BaseInterSF, public BaseIdxSF {
public:
  // Class type string
  static std::string _CT;
  // Parameter names
  static std::string _FILENAME; // Nmr restraint file name
  static std::string
      _QUADRATIC; // True = quadratic penalty function; false = linear

  NmrSF(const std::string &strName = "NMR");
  virtual ~NmrSF();

protected:
  virtual void SetupReceptor();
  virtual void SetupLigand();
  virtual void SetupScore();
  virtual double RawScore() const;
  void ParameterUpdated(const std::string &strName);

private:
  double NoeDistance(const NoeRestraintAtoms &noe) const;
  double StdDistance(const StdRestraintAtoms &std) const;

  bool m_bQuadratic; // synchronised with QUADRATIC named parameter
  NonBondedGridPtr m_spGrid;
  AtomList m_ligAtomList;          // All ligand atoms
  NoeRestraintAtomsList m_noeList; // List of all NOE interactions
  StdRestraintAtomsList m_stdList; // List of all STD interactions
};

} // namespace rxdock

#endif //_RBTNMRSF_H_
