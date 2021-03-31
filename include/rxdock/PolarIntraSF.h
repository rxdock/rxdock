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

// Ligand intramolecular scoring function for all attractive polar
// interactions (HBD,HBA,metal,guanidinium carbon)

#ifndef _RBTPOLARINTRASF_H_
#define _RBTPOLARINTRASF_H_

#include "rxdock/BaseIntraSF.h"
#include "rxdock/InteractionGrid.h"
#include "rxdock/PolarSF.h"

namespace rxdock {

class PolarIntraSF : public BaseIntraSF, public PolarSF {
public:
  // Class type string
  static const std::string _CT;
  // Parameter names
  static const std::string _ATTR;

  PolarIntraSF(const std::string &strName = "polar");
  virtual ~PolarIntraSF();

protected:
  virtual void SetupScore();
  virtual double RawScore() const;

  // Request Handling method
  // Handles the Partition request
  virtual void HandleRequest(RequestPtr spRequest);

  // Clear the model lists
  // As we are not using smart pointers, there is some memory management to do
  void ClearModel();

  // DM 25 Oct 2000 - track changes to parameter values in local data members
  // ParameterUpdated is invoked by ParamHandler::SetParameter
  void ParameterUpdated(const std::string &strName);

private:
  InteractionCenterList m_posList;
  InteractionCenterList m_negList;
  InteractionListMap m_intns;
  InteractionListMap m_prtIntns;
  bool m_bAttr;
};

} // namespace rxdock

#endif //_RBTPOLARINTRASF_H_
