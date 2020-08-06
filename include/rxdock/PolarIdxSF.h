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

// Indexed-grid-based intermolecular scoring function for all attractive polar
// interactions (HBD,HBA,metal,guanidinium carbon)

#ifndef _RBTPOLARIDXSF_H_
#define _RBTPOLARIDXSF_H_

#include "rxdock/BaseIdxSF.h"
#include "rxdock/BaseInterSF.h"
#include "rxdock/PolarSF.h"

namespace rxdock {

class PolarIdxSF : public BaseInterSF, public BaseIdxSF, public PolarSF {
public:
  // Class type string
  static std::string _CT;
  // Parameter names
  static std::string _INCR;
  static std::string _ATTR;
  // DM 12 Jun 2002 - score thresholds used for counting attractive and
  // repulsive interactions
  static std::string _THRESHOLD_POS;
  static std::string _THRESHOLD_NEG;

  RBTDLL_EXPORT static std::string &GetIncr();

  RBTDLL_EXPORT PolarIdxSF(const std::string &strName = "POLAR");
  virtual ~PolarIdxSF();

  // Override BaseSF::ScoreMap to provide additional raw descriptors
  virtual void ScoreMap(StringVariantMap &scoreMap) const;

protected:
  virtual void SetupReceptor();
  virtual void SetupLigand();
  virtual void SetupSolvent();
  virtual void SetupScore();
  virtual double RawScore() const;

  // Clear the receptor and ligand grids and lists respectively
  // As we are not using smart pointers, there is some memory management to do
  void ClearReceptor();
  void ClearLigand();
  void ClearSolvent();
  // Helper function for above
  void DeleteList(InteractionCenterList &icList);

  // DM 25 Oct 2000 - track changes to parameter values in local data members
  // ParameterUpdated is invoked by ParamHandler::SetParameter
  void ParameterUpdated(const std::string &strName);

private:
  double ReceptorScore() const;
  double SolventScore() const;
  double InterScore() const;
  double ReceptorSolventScore() const;
  double LigandSolventScore() const;

  double InterScore(const InteractionCenterList &posList,
                    const InteractionCenterList &negList, bool bCount) const;
  InteractionGridPtr m_spPosGrid;
  InteractionGridPtr m_spNegGrid;
  InteractionCenterList m_recepPosList;
  InteractionCenterList m_recepNegList;
  InteractionCenterList m_flexRecPosList;
  InteractionCenterList m_flexRecNegList;
  InteractionListMap m_flexRecIntns;
  InteractionListMap m_flexRecPrtIntns; // partitioned interactions

  InteractionCenterList m_ligPosList;
  InteractionCenterList m_ligNegList;

  InteractionCenterList m_solventPosList;
  InteractionCenterList m_solventNegList;
  InteractionListMap m_solventIntns;

  bool m_bAttr;
  bool m_bFlexRec;
  bool m_bSolvent;
  // DM 12 Jun 2002 - keep track of number of ligand atoms involved in non-zero
  // polar interactions
  mutable int m_nPos; //#positive centers with non-zero scores
  mutable int m_nNeg; //#negative centers with non-zero scores
  double m_posThreshold;
  double m_negThreshold;
};

} // namespace rxdock

#endif //_RBTPOLARIDXSF_H_
