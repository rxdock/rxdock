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

// Penalty function for keeping ligand and solvent within the docking volume
#ifndef _RBTCAVITYGRIDSF_H_
#define _RBTCAVITYGRIDSF_H_

#include "rxdock/BaseInterSF.h"
#include "rxdock/FlexDataVisitor.h"
#include "rxdock/RealGrid.h"

namespace rxdock {

class CavityGridSF : public BaseInterSF {
public:
  // FlexDataVisitor class to generate list of movable non-hydrogen atoms
  // i.e. those atoms that have the potential to move outside of the cavity
  class RBTDLL_EXPORT HeavyAtomFactory : public FlexDataVisitor {
  public:
    HeavyAtomFactory() {}
    HeavyAtomFactory(ModelList modelList);
    virtual void VisitReceptorFlexData(ReceptorFlexData *);
    virtual void VisitLigandFlexData(LigandFlexData *);
    virtual void VisitSolventFlexData(SolventFlexData *);
    AtomRList GetAtomList() const { return m_atomList; }

  private:
    AtomRList m_atomList;
  };

  // Class type string
  static std::string _CT;
  // Parameter names
  static std::string _RMAX;
  static std::string
      _QUADRATIC; // True = quadratic penalty function; false = linear

  CavityGridSF(const std::string &strName = "cavity");
  virtual ~CavityGridSF();

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
  RealGridPtr m_spGrid;
  double m_maxDist;     // Max distance of any grid point from the cavity
  AtomRList m_atomList; // combined list of all movable heavy atoms
  double m_rMax;
  bool m_bQuadratic; // synchronised with QUADRATIC named parameter
};

} // namespace rxdock

#endif //_RBTCAVITYGRIDSF_H_
