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

// Ligand intramolecular scoring function for vdw potential

#ifndef _RBTVDWINTRASF_H_
#define _RBTVDWINTRASF_H_

#include "rxdock/BaseIntraSF.h"
#include "rxdock/VdwSF.h"

#include <nlohmann/json.hpp>

using json = nlohmann::json;

namespace rxdock {

class VdwIntraSF : public BaseIntraSF, public VdwSF {
public:
  // Class type string
  static const std::string _CT;

  RBTDLL_EXPORT VdwIntraSF(const std::string &strName = "vdw");
  virtual ~VdwIntraSF();

  friend void to_json(json &j, const VdwIntraSF &vdwIntraSF);
  friend void from_json(const json &j, VdwIntraSF &vdwIntraSF);

  // Request Handling method
  // Handles the Partition request
  virtual void HandleRequest(RequestPtr spRequest);

protected:
  virtual void SetupScore();
  virtual double RawScore() const;

  // DM 25 Oct 2000 - track changes to parameter values in local data members
  // ParameterUpdated is invoked by ParamHandler::SetParameter
  void ParameterUpdated(const std::string &strName);

private:
  AtomRListList m_vdwIntns; // The full list of vdW interactions
  AtomRListList
      m_prtIntns; // The partitioned interactions (within partition distance)
  AtomRList m_ligAtomList;
};

void to_json(json &j, const VdwIntraSF &vdwIntraSF);
void from_json(const json &j, VdwIntraSF &vdwIntraSF);

} // namespace rxdock

#endif //_RBTVDWINTRASF_H_
