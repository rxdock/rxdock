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

// Tethered atoms distance restraint scoring function

#ifndef _RBTTETHERSF_H_
#define _RBTTETHERSF_H_

#include "rxdock/BaseInterSF.h"

#include <nlohmann/json.hpp>

using json = nlohmann::json;

namespace rxdock {

class TetherSF : public BaseInterSF {
public:
  // Class type string
  static const std::string _CT;
  // Parameter names
  static const std::string _REFERENCE_FILE;

  TetherSF(const std::string &strName = "tether");
  virtual ~TetherSF();

  friend void to_json(json &j, const TetherSF &tetsf);
  friend void from_json(const json &j, TetherSF &tetsf);

protected:
  virtual void SetupReceptor();
  virtual void SetupLigand();
  virtual void SetupScore();
  virtual double RawScore() const;
  // DM 25 Oct 2000 - track changes to parameter values in local data members
  // ParameterUpdated is invoked by ParamHandler::SetParameter
  void ParameterUpdated(const std::string &strName);

private:
  std::vector<int> ReadTetherAtoms(std::vector<std::string> &);
  AtomList m_ligAtomList;
  std::vector<int> m_tetherAtomList;
  CoordList m_tetherCoords;
};

void to_json(json &j, const TetherSF &tetsf);
void from_json(const json &j, TetherSF &tetsf);

} // namespace rxdock

#endif //_RBTTETHERSF_H_
