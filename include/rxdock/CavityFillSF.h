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

// Post-filtering scoring function for calculating percentage cavity unfilled
//(mode=0) or percentage ligand unbound (mode=1)

#ifndef _RBTCAVITYFILLSF_H_
#define _RBTCAVITYFILLSF_H_

#include "rxdock/BaseInterSF.h"
#include "rxdock/FFTGrid.h"

#include <nlohmann/json.hpp>

using json = nlohmann::json;

namespace rxdock {

class CavityFillSF : public BaseInterSF {
public:
  // Class type string
  static const std::string _CT;
  // Parameter names
  CavityFillSF(const std::string &strName = "cavity-fill");
  virtual ~CavityFillSF();
 
  friend void to_json(json &j, const CavityFillSF &sf);
  friend void from_json(const json &j, CavityFillSF &sf);

protected:
  virtual void SetupReceptor();
  virtual void SetupLigand();
  virtual void SetupScore();
  virtual double RawScore() const;
  void ParameterUpdated(const std::string &strName);

private:
  FFTGridPtr m_spGrid;
  AtomList m_ligAtomList;
};

void to_json(json &j, const CavityFillSF &sf);
void from_json(const json &j, CavityFillSF &sf);

} // namespace rxdock

#endif //_RBTCAVITYFILLSF_H_
