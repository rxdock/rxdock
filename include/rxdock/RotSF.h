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

#include "rxdock/BaseInterSF.h"

#include <nlohmann/json.hpp>

using json = nlohmann::json;

namespace rxdock {

class RotSF : public BaseInterSF {
public:
  // Static data member for class type
  static const std::string _CT;
  // Parameter names
  // Boolean controlling whether to include bonds to NH3+ in rotable bond count
  static const std::string _INCNH3;
  // Boolean controlling whether to include bonds to OH in  rotable bond count
  static const std::string _INCOH;

  RotSF(const std::string &strName = "rot");
  virtual ~RotSF();

  friend void to_json(json &j, const RotSF &rotsf);
  friend void from_json(const json &j, RotSF &rotsf);

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

void to_json(json &j, const RotSF &rotsf);
void from_json(const json &j, RotSF &rotsf);

} // namespace rxdock

#endif //_RBTROTSF_H_
