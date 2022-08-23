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

// Base implementation class for all attractive and repulsive polar scoring fns
// Provides methods for: creating interaction centers, calculating interaction
// scores

#ifndef _RBTPOLARSF_H_
#define _RBTPOLARSF_H_

#include "rxdock/AnnotationHandler.h"
#include "rxdock/BaseSF.h"
#include "rxdock/InteractionGrid.h"

#include <nlohmann/json.hpp>

using json = nlohmann::json;

namespace rxdock {

class PolarSF : public virtual BaseSF, public virtual AnnotationHandler {
public:
  // Class type string
  static const std::string _CT;
  // Parameter names
  static const std::string _INCR;
  static const std::string _R12FACTOR;
  static const std::string _R12INCR;
  static const std::string _DR12MIN;
  static const std::string _DR12MAX;
  static const std::string _A1;
  static const std::string _DA1MIN;
  static const std::string _DA1MAX;
  static const std::string _A2;
  static const std::string _DA2MIN;
  static const std::string _DA2MAX;
  static const std::string _INCMETAL;
  static const std::string _INCHBD;
  static const std::string _INCHBA;
  static const std::string _INCGUAN;
  static const std::string _GUAN_PLANE;
  static const std::string _ABS_DR12;
  static const std::string _LP_OSP2;
  static const std::string _LP_PHI;
  static const std::string _LP_DPHIMIN;
  static const std::string _LP_DPHIMAX;
  static const std::string _LP_DTHETAMIN;
  static const std::string _LP_DTHETAMAX;

  virtual ~PolarSF();

  friend void to_json(json &j, const PolarSF &polarSF);
  friend void from_json(const json &j, PolarSF &polarSF);

protected:
  PolarSF();

  InteractionCenterList
  CreateAcceptorInteractionCenters(const AtomList &atomList) const;
  InteractionCenterList
  CreateDonorInteractionCenters(const AtomList &atomList) const;

  // Index the intramolecular interactions between two lists
  void BuildIntraMap(const InteractionCenterList &ICList1,
                     const InteractionCenterList &ICList2,
                     InteractionListMap &intns) const;
  // Index the intramolecular interactions within a single list
  void BuildIntraMap(const InteractionCenterList &ICList,
                     InteractionListMap &intns) const;

  double IntraScore(const InteractionCenterList &posList,
                    const InteractionCenterList &negList,
                    const InteractionListMap &prtIntns, bool attr) const;
  void Partition(const InteractionCenterList &posList,
                 const InteractionCenterList &negList,
                 const InteractionListMap &intns, InteractionListMap &prtIntns,
                 double dist = 0.0) const;

  // Generic scoring function params
  struct f1prms {
    double R0, DRMin, DRMax, slope;
    f1prms(double R, double DMin, double DMax)
        : R0(R), DRMin(DMin), DRMax(DMax), slope(1.0 / (DMax - DMin)) {}
  };

  inline f1prms GetRprms() const { return f1prms(0.0, m_DR12Min, m_DR12Max); }
  inline f1prms GetA1prms() const { return f1prms(m_A1, m_DA1Min, m_DA1Max); }
  inline f1prms GetA2prms() const { return f1prms(m_A2, m_DA2Min, m_DA2Max); }

  double PolarScore(const InteractionCenter *intn,
                    const InteractionCenterList &intnList, const f1prms &Rprms,
                    const f1prms &A1prms, const f1prms &A2prms) const;

  // As this has a virtual base class we need a separate OwnParameterUpdated
  // which can be called by concrete subclass ParameterUpdated methods
  // See Stroustrup C++ 3rd edition, p395, on programming virtual base classes
  void OwnParameterUpdated(const std::string &strName);

private:
  // Generic scoring function primitive
  inline double f1(double DR, const f1prms &prms) const {
    return (DR > prms.DRMax)   ? 0.0
           : (DR > prms.DRMin) ? 1.0 - prms.slope * (DR - prms.DRMin)
                               : 1.0;
  }

  void UpdateLPprms();

  // DM 25 Oct 2000 - heavily used params
  double m_R12Factor;
  double m_R12Incr;
  double m_DR12Min;
  double m_DR12Max;
  double m_A1;
  double m_DA1Min;
  double m_DA1Max;
  double m_A2;
  double m_DA2Min;
  double m_DA2Max;
  bool m_bAbsDR12;
  double m_LP_PHI;
  double m_LP_DPHIMin;
  double m_LP_DPHIMax;
  double m_LP_DTHETAMin;
  double m_LP_DTHETAMax;
  f1prms m_PHI_lp_prms;
  f1prms m_PHI_plane_prms;
  f1prms m_THETAprms;
};

void to_json(json &j, const PolarSF &polarSF);
void from_json(const json &j, PolarSF &polarSF);

} // namespace rxdock

#endif //_RBTPOLARSF_H_
