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

#include "RbtBaseInterSF.h"
#include "RbtRealGrid.h"

namespace rxdock {

class RbtVdwGridSF : public RbtBaseInterSF {
public:
  // Class type string
  static std::string _CT;

  RBTDLL_EXPORT static std::string &GetCt();

  // Parameter names
  static std::string _GRID;     // Suffix for grid filename
  static std::string _SMOOTHED; // Controls whether to smooth the grid values

  RbtVdwGridSF(const std::string &strName = "VDW");
  virtual ~RbtVdwGridSF();

protected:
  virtual void SetupReceptor();
  virtual void SetupLigand();
  virtual void SetupSolvent();
  virtual void SetupScore();
  virtual double RawScore() const;
  // DM 25 Oct 2000 - track changes to parameter values in local data members
  // ParameterUpdated is invoked by RbtParamHandler::SetParameter
  void ParameterUpdated(const std::string &strName);

private:
  // Read grids from input stream
  void ReadGrids(std::istream &istr);

  RbtRealGridList m_grids;
  RbtAtomRList m_ligAtomList;
  RbtTriposAtomTypeList m_ligAtomTypes;
  bool m_bSmoothed;
};

} // namespace rxdock

#endif //_RBTVDWGRIDSF_H_
