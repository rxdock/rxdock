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

#include "RbtAnnotationHandler.h"
#include "RbtBaseSF.h"
#include "RbtInteractionGrid.h"

namespace rxdock {

class RbtPolarSF : public virtual RbtBaseSF,
                   public virtual RbtAnnotationHandler {
public:
  // Class type string
  static std::string _CT;
  // Parameter names
  static std::string _INCR;
  static std::string _R12FACTOR;
  static std::string _R12INCR;
  static std::string _DR12MIN;
  static std::string _DR12MAX;
  static std::string _A1;
  static std::string _DA1MIN;
  static std::string _DA1MAX;
  static std::string _A2;
  static std::string _DA2MIN;
  static std::string _DA2MAX;
  static std::string _INCMETAL;
  static std::string _INCHBD;
  static std::string _INCHBA;
  static std::string _INCGUAN;
  static std::string _GUAN_PLANE;
  static std::string _ABS_DR12;
  static std::string _LP_OSP2;
  static std::string _LP_PHI;
  static std::string _LP_DPHIMIN;
  static std::string _LP_DPHIMAX;
  static std::string _LP_DTHETAMIN;
  static std::string _LP_DTHETAMAX;

  virtual ~RbtPolarSF();

protected:
  RbtPolarSF();

  RbtInteractionCenterList
  CreateAcceptorInteractionCenters(const RbtAtomList &atomList) const;
  RbtInteractionCenterList
  CreateDonorInteractionCenters(const RbtAtomList &atomList) const;

  // Index the intramolecular interactions between two lists
  void BuildIntraMap(const RbtInteractionCenterList &ICList1,
                     const RbtInteractionCenterList &ICList2,
                     RbtInteractionListMap &intns) const;
  // Index the intramolecular interactions within a single list
  void BuildIntraMap(const RbtInteractionCenterList &ICList,
                     RbtInteractionListMap &intns) const;

  double IntraScore(const RbtInteractionCenterList &posList,
                    const RbtInteractionCenterList &negList,
                    const RbtInteractionListMap &prtIntns, bool attr) const;
  void Partition(const RbtInteractionCenterList &posList,
                 const RbtInteractionCenterList &negList,
                 const RbtInteractionListMap &intns,
                 RbtInteractionListMap &prtIntns, double dist = 0.0) const;

  // Generic scoring function params
  struct f1prms {
    double R0, DRMin, DRMax, slope;
    f1prms(double R, double DMin, double DMax)
        : R0(R), DRMin(DMin), DRMax(DMax), slope(1.0 / (DMax - DMin)) {}
  };

  inline f1prms GetRprms() const { return f1prms(0.0, m_DR12Min, m_DR12Max); }
  inline f1prms GetA1prms() const { return f1prms(m_A1, m_DA1Min, m_DA1Max); }
  inline f1prms GetA2prms() const { return f1prms(m_A2, m_DA2Min, m_DA2Max); }

  double PolarScore(const RbtInteractionCenter *intn,
                    const RbtInteractionCenterList &intnList,
                    const f1prms &Rprms, const f1prms &A1prms,
                    const f1prms &A2prms) const;

  // As this has a virtual base class we need a separate OwnParameterUpdated
  // which can be called by concrete subclass ParameterUpdated methods
  // See Stroustrup C++ 3rd edition, p395, on programming virtual base classes
  void OwnParameterUpdated(const std::string &strName);

private:
  // Generic scoring function primitive
  inline double f1(double DR, const f1prms &prms) const {
    return (DR > prms.DRMax)
               ? 0.0
               : (DR > prms.DRMin) ? 1.0 - prms.slope * (DR - prms.DRMin) : 1.0;
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

} // namespace rxdock

#endif //_RBTPOLARSF_H_
