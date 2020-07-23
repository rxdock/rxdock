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

// Pharmacophore distance restraint scoring function

#ifndef _RBTPHARMASF_H_
#define _RBTPHARMASF_H_

#include "RbtBaseInterSF.h"
#include "RbtBaseMolecularFileSink.h"
#include "RbtConstraint.h"

namespace rxdock {

class RbtPharmaSF : public RbtBaseInterSF {
public:
  // Class type string
  static std::string _CT;
  // Parameter names
  static std::string _CONSTRAINTS_FILE;
  static std::string _OPTIONAL_FILE;
  static std::string _NOPT;
  static std::string _WRITE_ERRORS;

  RbtPharmaSF(const std::string &strName = "PHARMA");
  virtual ~RbtPharmaSF();
  // Override RbtBaseSF::ScoreMap to provide additional raw descriptors
  virtual void ScoreMap(RbtStringVariantMap &scoreMap) const;

protected:
  virtual void SetupReceptor();
  virtual void SetupLigand();
  virtual void SetupScore();
  virtual double RawScore() const;
  // DM 25 Oct 2000 - track changes to parameter values in local data members
  // ParameterUpdated is invoked by RbtParamHandler::SetParameter
  void ParameterUpdated(const std::string &strName);

private:
  RbtConstraintList m_constrList;
  RbtConstraintList m_optList;
  int m_nopt;
  RbtMolecularFileSinkPtr m_spErrorFile;
  bool m_bWriteErrors;
  // Keep track of individual constraint scores for ScoreMap
  mutable std::vector<double> m_conScores; // Mandatory constraint scores
  mutable std::vector<double> m_optScores; // Optional constraint scores
};

} // namespace rxdock

#endif //_RBTPHARMASF_H_
