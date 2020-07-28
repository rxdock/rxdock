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

// Precalculated-grid-based intermolecular vdw scoring function

#ifndef _RBTVDWGRIDSF_H_
#define _RBTVDWGRIDSF_H_

#include "BaseInterSF.h"
#include "RealGrid.h"

namespace rxdock {

class VdwGridSF : public BaseInterSF {
public:
  // Class type string
  static std::string _CT;

  RBTDLL_EXPORT static std::string &GetCt();

  // Parameter names
  static std::string _GRID;     // Suffix for grid filename
  static std::string _SMOOTHED; // Controls whether to smooth the grid values

  VdwGridSF(const std::string &strName = "VDW");
  virtual ~VdwGridSF();

protected:
  virtual void SetupReceptor();
  virtual void SetupLigand();
  virtual void SetupSolvent();
  virtual void SetupScore();
  virtual double RawScore() const;
  // DM 25 Oct 2000 - track changes to parameter values in local data members
  // ParameterUpdated is invoked by ParamHandler::SetParameter
  void ParameterUpdated(const std::string &strName);

private:
  // Read grids from input stream
  void ReadGrids(json vdwGrids);

  RealGridList m_grids;
  AtomRList m_ligAtomList;
  TriposAtomTypeList m_ligAtomTypes;
  bool m_bSmoothed;
};

} // namespace rxdock

#endif //_RBTVDWGRIDSF_H_
