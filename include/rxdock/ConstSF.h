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

// Scoring function to penalise ligand and solvent binding
// Ligand binding RawScore is always 1
// Solvent binding RawScore is:
// SOLVENT_PENALTY * the number of enabled solvent models
#ifndef _RBTCONSTSF_H_
#define _RBTCONSTSF_H_

#include "rxdock/BaseInterSF.h"

namespace rxdock {

class ConstSF : public BaseInterSF {
public:
  // Static data member for class type (i.e. "ConstSF")
  static std::string _CT;
  static std::string _SOLVENT_PENALTY;

  ConstSF(const std::string &strName = "CONST");
  virtual ~ConstSF();

  virtual void ScoreMap(StringVariantMap &scoreMap) const;

protected:
  virtual void SetupReceptor() {}
  virtual void SetupLigand() {}
  virtual void SetupScore() {}
  virtual double RawScore() const;
  void ParameterUpdated(const std::string &strName);

private:
  // The original constant score for ligand binding
  double InterScore() const;
  // The solvent binding penalty
  double SystemScore() const;
  double m_solventPenalty;
};

} // namespace rxdock

#endif //_RBTCONSTSF_H_
