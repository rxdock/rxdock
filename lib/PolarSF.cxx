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

#include "rxdock/PolarSF.h"
#include "rxdock/Plane.h"
#include "rxdock/WorkSpace.h"

#include <loguru.hpp>

#include <functional>

using namespace rxdock;

// Static data members
const std::string PolarSF::_CT = "PolarSF";
const std::string PolarSF::_R12FACTOR = "R12FACTOR";
const std::string PolarSF::_R12INCR = "R12INCR";
const std::string PolarSF::_DR12MIN = "DR12MIN";
const std::string PolarSF::_DR12MAX = "DR12MAX";
const std::string PolarSF::_A1 = "A1";
const std::string PolarSF::_DA1MIN = "DA1MIN";
const std::string PolarSF::_DA1MAX = "DA1MAX";
const std::string PolarSF::_A2 = "A2";
const std::string PolarSF::_DA2MIN = "DA2MIN";
const std::string PolarSF::_DA2MAX = "DA2MAX";
const std::string PolarSF::_INCMETAL = "INCMETAL";
const std::string PolarSF::_INCHBD = "INCHBD";
const std::string PolarSF::_INCHBA = "INCHBA";
const std::string PolarSF::_INCGUAN = "INCGUAN";
const std::string PolarSF::_GUAN_PLANE = "GUAN_PLANE";
const std::string PolarSF::_ABS_DR12 = "ABS_DR12";
const std::string PolarSF::_LP_OSP2 = "LP_OSP2";
const std::string PolarSF::_LP_PHI = "LP_PHI";
const std::string PolarSF::_LP_DPHIMIN = "LP_DPHIMIN";
const std::string PolarSF::_LP_DPHIMAX = "LP_DPHIMAX";
const std::string PolarSF::_LP_DTHETAMIN = "LP_DTHETAMIN";
const std::string PolarSF::_LP_DTHETAMAX = "LP_DTHETAMAX";

PolarSF::PolarSF()
    : m_R12Factor(1.0), m_R12Incr(0.6), m_DR12Min(0.25), m_DR12Max(0.6),
      m_A1(180.0), m_DA1Min(30.0), m_DA1Max(80.0), m_A2(150.0), m_DA2Min(30.0),
      m_DA2Max(70.0), m_bAbsDR12(true), m_LP_PHI(45.0), m_LP_DPHIMin(15.0),
      m_LP_DPHIMax(30.0), m_LP_DTHETAMin(20.0), m_LP_DTHETAMax(60.0),
      m_PHI_lp_prms(m_LP_PHI, m_LP_DPHIMin, m_LP_DPHIMax),
      m_PHI_plane_prms(0.0, m_LP_PHI + m_LP_DPHIMin, m_LP_PHI + m_LP_DPHIMax),
      m_THETAprms(0.0, m_LP_DTHETAMin, m_LP_DTHETAMax) {
  LOG_F(2, "PolarSF default constructor");
  // Add parameters
  AddParameter(_R12FACTOR, m_R12Factor);
  AddParameter(_R12INCR, m_R12Incr);
  AddParameter(_DR12MIN, m_DR12Min);
  AddParameter(_DR12MAX, m_DR12Max);
  AddParameter(_A1, m_A1);
  AddParameter(_DA1MIN, m_DA1Min);
  AddParameter(_DA1MAX, m_DA1Max);
  AddParameter(_A2, m_A2);
  AddParameter(_DA2MIN, m_DA2Min);
  AddParameter(_DA2MAX, m_DA2Max);
  AddParameter(_INCMETAL, true);
  AddParameter(_INCHBD, true);
  AddParameter(_INCHBA, true);
  AddParameter(_INCGUAN, true);
  AddParameter(_GUAN_PLANE, true);
  AddParameter(_ABS_DR12, m_bAbsDR12);
  AddParameter(_LP_OSP2, false);
  AddParameter(_LP_PHI, m_LP_PHI);
  AddParameter(_LP_DPHIMIN, m_LP_DPHIMin);
  AddParameter(_LP_DPHIMAX, m_LP_DPHIMax);
  AddParameter(_LP_DTHETAMIN, m_LP_DTHETAMin);
  AddParameter(_LP_DTHETAMAX, m_LP_DTHETAMax);
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

PolarSF::~PolarSF() {
  LOG_F(2, "PolarSF destructor");
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

InteractionCenterList
PolarSF::CreateDonorInteractionCenters(const AtomList &atomList) const {
  const bool bIncHBD = GetParameter(_INCHBD);
  const bool bIncMetal = GetParameter(_INCMETAL);
  const bool bIncGuan = GetParameter(_INCGUAN);
  const bool bGuanPlane = GetParameter(_GUAN_PLANE);

  InteractionCenterList intnList;

  // H-Bond donors
  if (bIncHBD) {
    AtomList posList = GetAtomListWithPredicate(atomList, isAtomHBondDonor());
    for (AtomListConstIter iter = posList.begin(); iter != posList.end();
         iter++) {
      AtomList parentList = GetBondedAtomList(*iter);
      // Create interaction center from donor H and donor parent
      intnList.push_back(new InteractionCenter(*iter, parentList.front()));
    }
  }

  // Metals
  if (bIncMetal) {
    AtomList posList = GetAtomListWithPredicate(atomList, isAtomMetal());
    for (AtomListConstIter iter = posList.begin(); iter != posList.end();
         iter++) {
      // Create interaction center from metal ion only
      intnList.push_back(new InteractionCenter(*iter));
    }
  }

  // Guanidinium carbons
  if (bIncGuan) {
    AtomList posList =
        GetAtomListWithPredicate(atomList, isAtomGuanidiniumCarbon());
    for (AtomListConstIter iter = posList.begin(); iter != posList.end();
         iter++) {
      AtomList parentList = GetBondedAtomList(*iter);
      // Two options: for attractive interactions (HBA..GUAN) the SF depends on
      // the angle between the acceptor and the normal to the plane of the
      // guanidinium, so we define an interaction center with 3 atoms
      if (bGuanPlane) {
        // Create interaction center from guanidinium carbon and first 2 bonded
        // atoms
        intnList.push_back(
            new InteractionCenter(*iter, parentList[0], parentList[1]));
      }
      // For repulsive interactions (HBD..GUAN or GUAN..GUAN) the SF depends
      // only on distance, so we define an interaction center with just 1 atom
      else {
        // Create interaction center from guanidinium carbon only
        intnList.push_back(new InteractionCenter(*iter));
      }
    }
  }

  return intnList;
}

InteractionCenterList
PolarSF::CreateAcceptorInteractionCenters(const AtomList &atomList) const {
  const bool bIncHBA = GetParameter(_INCHBA);
  const bool bLP = GetParameter(_LP_OSP2);

  InteractionCenterList intnList;
  isAtomicNo_eq bIsC(6);
  isAtomicNo_eq bIsN(7);
  isAtomicNo_eq bIsO(8);
  isAtomRNA bIsRNA;
  isAtomAnionic bIsAnionic;

  // Negative polar (HBA)
  if (bIncHBA) {
    AtomList negList =
        GetAtomListWithPredicate(atomList, isAtomHBondAcceptor());
    for (AtomListConstIter iter = negList.begin(); iter != negList.end();
         iter++) {
      AtomList parentList = GetBondedAtomList(*iter);

      // Unconnected acceptor - ignore
      if (parentList.empty())
        continue;

      // Singly connected acceptor
      if (parentList.size() == 1) {
        AtomPtr spAcceptorParent = parentList.front();
        // Detect terminal oxygen bonded to carbon or nitrogen here
        // and define an interaction center including lone pair geometry if
        // required
        if (bLP && bIsO(*iter) &&
            (bIsC(spAcceptorParent) || bIsN(spAcceptorParent))) {
          // Get the atoms bonded to the acceptor parent (not including the
          // acceptor itself)
          AtomList grandParentList =
              GetAtomListWithPredicate(GetBondedAtomList(spAcceptorParent),
                                       std::bind(std::not2(isAtomPtr_eq()),
                                                 std::placeholders::_1, *iter));
          if (!grandParentList.empty()) {
            AtomPtr spGrandParent = grandParentList.front();
            // If O is -ve charged (likely carboxylate), or is in RNA, or is in
            // nitro group define interaction center with LONEPAIR geometry
            // (directionality and in-plane of LP)
            if (bIsAnionic(*iter) || bIsRNA(*iter) || bIsN(spAcceptorParent)) {
              intnList.push_back(
                  new InteractionCenter(*iter, spAcceptorParent, spGrandParent,
                                        InteractionCenter::LONEPAIR));
              LOG_F(1, "LONEPAIR acceptor: {}", (*iter)->GetFullAtomName());
            }
            // else define interaction center with PLANE geometry (in-plane of
            // LP, broader directionality)
            else {
              intnList.push_back(
                  new InteractionCenter(*iter, spAcceptorParent, spGrandParent,
                                        InteractionCenter::PLANE));
              LOG_F(1, "PLANE acceptor: {}", (*iter)->GetFullAtomName());
            }
          } else {
            intnList.push_back(new InteractionCenter(*iter, spAcceptorParent));
          }
        } else {
          intnList.push_back(new InteractionCenter(*iter, spAcceptorParent));
        }
      } else {
        // define pseudo atom for acceptor parent if more than one parent
        // pseudo atoms stored permanently in receptor model
        PseudoAtomPtr spPseudoAtom =
            (*iter)->GetModelPtr()->AddPseudoAtom(parentList);
        AtomPtr spAcceptorParent = spPseudoAtom;
        intnList.push_back(new InteractionCenter(*iter, spAcceptorParent));
      }
    }
  }
  return intnList;
}

void PolarSF::BuildIntraMap(const InteractionCenterList &ICList1,
                            const InteractionCenterList &ICList2,
                            InteractionListMap &intns) const {
  SelectInteractionCenter selectIC(true);
  isInteractionCenterSelected isSelected;
  bool bSingleList = ICList2.empty(); // If true, then index the flexible
                                      // interactions within ICList1
  for (InteractionCenterListConstIter iter = ICList1.begin();
       iter != ICList1.end(); iter++) {
    // 1. First select all atoms of all the interaction centers in the list
    //(Needed to deal with situations where the interaction centers do not all
    // belong to the same model)
    std::for_each(ICList1.begin(), ICList1.end(), selectIC);
    // We use the ID of Atom1 of the interaction center (IC) as the index into
    // the map but we check all the constituent atoms of the IC to assess
    // whether the interaction should be deemed "flexible" Reason is that some
    // interactions have an angular dependence: the primary distance may be
    // rigid, but the angular dependence may not be
    Atom *pAtom = (*iter)->GetAtom1Ptr();
    unsigned int id = pAtom->GetAtomId() - 1;
    Assert<Assertion>(id < intns.size());
    // 2. Deselect all atoms in the same model as the current atom
    pAtom->GetModelPtr()->SetAtomSelectionFlags(false);
    // 3. Set the selection flags for all atoms in the model whose distance can
    // vary to *any* of the atoms in the IC
    AtomRList atomList = (*iter)->GetAtomList();
    std::for_each(atomList.begin(), atomList.end(), SelectFlexAtoms());
    // Now we retrieve the ICs from the 2nd list whose constituent atoms have
    // been selected
    if (bSingleList) {
      std::copy_if(iter + 1, ICList1.end(), std::back_inserter(intns[id]),
                   isSelected);
    } else {
      std::copy_if(ICList2.begin(), ICList2.end(),
                   std::back_inserter(intns[id]), isSelected);
    }
    LOG_F(1, "{}: Intns to {}", GetFullName(), pAtom->GetFullAtomName());
    for (InteractionCenterListConstIter i = intns[id].begin();
         i != intns[id].end(); i++) {
      LOG_F(1, "    {}", (*i)->GetAtom1Ptr()->GetFullAtomName());
    }
  }
}

// Convenience method for above.
// Indexes the flexible interactions within a single interaction center list.
void PolarSF::BuildIntraMap(const InteractionCenterList &ICList,
                            InteractionListMap &intns) const {
  InteractionCenterList emptyList;
  BuildIntraMap(ICList, emptyList, intns);
}

double PolarSF::IntraScore(const InteractionCenterList &posList,
                           const InteractionCenterList &negList,
                           const InteractionListMap &intns, bool attr) const {
  double score = 0.0; // Total score

  PolarSF::f1prms Rprms = GetRprms();   // Distance params
  PolarSF::f1prms A1prms = GetA1prms(); // Donor angle params
  PolarSF::f1prms A2prms = GetA2prms(); // Acceptor angle params

  double s(0.0); // Partial scores

  if (attr) {
    // Pos-neg
    for (InteractionCenterListConstIter iter = posList.begin();
         iter != posList.end(); iter++) {
      Atom *pAtom = (*iter)->GetAtom1Ptr();
      int id = pAtom->GetAtomId() - 1;
      // NO CHECK ON ID IN RANGE
      s = PolarScore(*iter, intns[id], Rprms, A1prms, A2prms);
      if (std::fabs(s) > 0) {
        LOG_F(1, "{}: pos-neg {}; s={}", GetFullName(),
              pAtom->GetFullAtomName(), s);
      }
      score += pAtom->GetUser1Value() * s;
    }
    // Neg-pos
    for (InteractionCenterListConstIter iter = negList.begin();
         iter != negList.end(); iter++) {
      Atom *pAtom = (*iter)->GetAtom1Ptr();
      int id = pAtom->GetAtomId() - 1;
      // NO CHECK ON ID IN RANGE
      s = PolarScore(*iter, intns[id], Rprms, A2prms, A1prms);
      if (std::fabs(s) > 0) {
        LOG_F(1, "{}: neg-pos {}; s={}", GetFullName(),
              pAtom->GetFullAtomName(), s);
      }
      score += pAtom->GetUser1Value() * s;
    }
  } else {
    // pos-pos
    for (InteractionCenterListConstIter iter = posList.begin();
         iter != posList.end(); iter++) {
      Atom *pAtom = (*iter)->GetAtom1Ptr();
      int id = pAtom->GetAtomId() - 1;
      // NO CHECK ON ID IN RANGE
      s = PolarScore(*iter, intns[id], Rprms, A1prms, A1prms);
      if (std::fabs(s) > 0) {
        LOG_F(1, "{}: pos-pos {}; s={}", GetFullName(),
              pAtom->GetFullAtomName(), s);
      }
      score += pAtom->GetUser1Value() * s;
    }
    // neg-neg
    for (InteractionCenterListConstIter iter = negList.begin();
         iter != negList.end(); iter++) {
      Atom *pAtom = (*iter)->GetAtom1Ptr();
      int id = pAtom->GetAtomId() - 1;
      // NO CHECK ON ID IN RANGE
      s = PolarScore(*iter, intns[id], Rprms, A2prms, A2prms);
      if (std::fabs(s) > 0) {
        LOG_F(1, "{}: neg-neg {}; s={}", GetFullName(),
              pAtom->GetFullAtomName(), s);
      }
      score += pAtom->GetUser1Value() * s;
    }
  }
  return score;
}

void PolarSF::Partition(const InteractionCenterList &posList,
                        const InteractionCenterList &negList,
                        const InteractionListMap &intns,
                        InteractionListMap &prtIntns, double dist) const {
  // Clear the existing partitioned lists
  for (InteractionListMapIter iter = prtIntns.begin(); iter != prtIntns.end();
       iter++)
    (*iter).clear();
  LOG_F(1, "Partition (dist={})", dist);

  for (InteractionCenterListConstIter iter = posList.begin();
       iter != posList.end(); iter++) {
    Atom *pAtom = (*iter)->GetAtom1Ptr();
    int id = pAtom->GetAtomId() - 1;
    if (dist > 0.0) {
      isInteractionD_lt bInRange(*iter, dist);
      // This copies all interactions which are within dist A of atom
      std::copy_if(intns[id].begin(), intns[id].end(),
                   std::back_inserter(prtIntns[id]), bInRange);
    } else {
      // Remove partition - simply copy from the master list of vdW interactions
      prtIntns[id] = intns[id];
    }
    LOG_F(1, "{}: #intn={}: #prtn={}",
          (*iter)->GetAtom1Ptr()->GetFullAtomName(), intns[id].size(),
          prtIntns[id].size());
  }

  for (InteractionCenterListConstIter iter = negList.begin();
       iter != negList.end(); iter++) {
    Atom *pAtom = (*iter)->GetAtom1Ptr();
    int id = pAtom->GetAtomId() - 1;
    if (dist > 0.0) {
      isInteractionD_lt bInRange(*iter, dist);
      // This copies all interactions which are within dist A of atom
      std::copy_if(intns[id].begin(), intns[id].end(),
                   std::back_inserter(prtIntns[id]), bInRange);
    } else {
      // Remove partition - simply copy from the master list of vdW interactions
      prtIntns[id] = intns[id];
    }
    LOG_F(1, "{}: #intn={}: #prtn={}",
          (*iter)->GetAtom1Ptr()->GetFullAtomName(), intns[id].size(),
          prtIntns[id].size());
  }
}

double PolarSF::PolarScore(const InteractionCenter *pIC1,
                           const InteractionCenterList &IC2List,
                           const f1prms &Rprms, const f1prms &A1prms,
                           const f1prms &A2prms) const {
  double s(0.0);
  if (IC2List.empty()) {
    return s;
  }

  Coord nullCoord;
  bool bAnnotate = isAnnotationEnabled();
  // std::string strName = GetName();

  Atom *pAtom1_1 = pIC1->GetAtom1Ptr();
  // DM 7 June 2006 - check for disabled interaction centre 1
  if (!pAtom1_1->GetEnabled()) {
    return s;
  }
  const Coord &cAtom1_1 = pAtom1_1->GetCoords();
  Atom *pAtom1_2 = pIC1->GetAtom2Ptr();
  Atom *pAtom1_3 = pIC1->GetAtom3Ptr();
  InteractionCenter::eLP eLP1 = pIC1->LP();
  bool bAngle1 =
      ((pAtom1_2 != nullptr) && (pAtom1_3 == nullptr)) ? true : false;
  bool bPlane1 = ((pAtom1_2 != nullptr) && (pAtom1_3 != nullptr) &&
                  (eLP1 == InteractionCenter::NONE))
                     ? true
                     : false;
  bool bLP1 = ((pAtom1_2 != nullptr) && (pAtom1_3 != nullptr) &&
               (eLP1 != InteractionCenter::NONE))
                  ? true
                  : false;
  const f1prms &PHI1prms =
      (eLP1 == InteractionCenter::LONEPAIR) ? m_PHI_lp_prms : m_PHI_plane_prms;
  const Coord &cAtom1_2 =
      (bAngle1 || bPlane1 || bLP1) ? pAtom1_2->GetCoords() : nullCoord;
  const Coord &cAtom1_3 = (bPlane1 || bLP1) ? pAtom1_3->GetCoords() : nullCoord;
  Plane pl1 = (bPlane1 || bLP1) ? Plane(cAtom1_1, cAtom1_2, cAtom1_3) : Plane();
  double radius1 = pAtom1_1->GetVdwRadius();

  for (InteractionCenterListConstIter IC2Iter = IC2List.begin();
       IC2Iter != IC2List.end(); IC2Iter++) {
    Atom *pAtom2_1 = (*IC2Iter)->GetAtom1Ptr();
    // if (pAtom1_1 == pAtom2_1) continue;//check for self-interactions
    if (!pAtom2_1->GetEnabled())
      continue; // check for disabled interaction centre 2
    const Coord &cAtom2_1 = pAtom2_1->GetCoords();
    Atom *pAtom2_2 = (*IC2Iter)->GetAtom2Ptr();
    Atom *pAtom2_3 = (*IC2Iter)->GetAtom3Ptr();
    InteractionCenter::eLP eLP2 = (*IC2Iter)->LP();
    bool bAngle2 =
        ((pAtom2_2 != nullptr) && (pAtom2_3 == nullptr)) ? true : false;
    bool bPlane2 = ((pAtom2_2 != nullptr) && (pAtom2_3 != nullptr) &&
                    (eLP2 == InteractionCenter::NONE))
                       ? true
                       : false;
    bool bLP2 = ((pAtom2_2 != nullptr) && (pAtom2_3 != nullptr) &&
                 (eLP2 != InteractionCenter::NONE))
                    ? true
                    : false;
    double radius2 = pAtom2_1->GetVdwRadius();
    double R12 = m_R12Factor * (radius1 + radius2) + m_R12Incr;
    Vector v12 = cAtom1_1 - cAtom2_1;
    double R = v12.Length();
    double DR = R - R12;
    double f = m_bAbsDR12 ? f1(std::fabs(DR), Rprms) : f1(DR, Rprms);
    if (f > 0.0) {
      LOG_F(1, "R12: {}, {}, {}, {}, f={}", pAtom1_1->GetFullAtomName(),
            pAtom2_1->GetFullAtomName(), R12, DR, f);
      // For guanidinium bPlane2 interacting with C=O lone pair bLP1, we want to
      // use the regular angular dependence
      if (bAngle1 || (bPlane2 && bLP1)) {
        double DA1 = Angle(cAtom1_2, cAtom1_1, cAtom2_1) - A1prms.R0;
        f *= f1(std::fabs(DA1), A1prms);
        LOG_F(1, "A1: {}, {}, f={}", A1prms.R0, DA1, f);
      } else if (bPlane1) {
        double A =
            std::acos(-std::fabs(Dot(v12.Unit(), pl1.VNorm()))) * 180.0 / M_PI;
        double DA1 = A - A1prms.R0;
        f *= f1(std::fabs(DA1), A1prms);
        LOG_F(1, "Pl1: {}, {}, f={}", A1prms.R0, DA1, f);
      }
      // Lone pair geometry dependence around Osp2
      else if (bLP1) {
        // Perpendicular distance from donor H to plane of lone pairs
        double dPerp = DistanceFromPointToPlane(cAtom2_1, pl1);
        // Coordinate of donor H projected into plane of lone pairs
        Coord cPerp = cAtom2_1 - dPerp * pl1.VNorm();
        LOG_F(1, "dPerp = {} check cPerp = {}", dPerp,
              DistanceFromPointToPlane(cPerp, pl1));
        // Can calculate std::sin(theta) directly from the two distances we have
        // already
        double theta = std::asin(dPerp / R) * 180.0 / M_PI;
        f *= f1(std::fabs(theta), m_THETAprms);
        if (f > 0.0) {
          double phi = 180.0 - Angle(cPerp, cAtom1_1, cAtom1_2);
          double Dphi = phi - PHI1prms.R0;
          f *= f1(std::fabs(Dphi), PHI1prms);
          LOG_F(1, "LP1: {}, {}, {}, f={}", theta, PHI1prms.R0, Dphi, f);
        }
      }
      if (f > 0.0) {
        // For guanidinium bPlane1 interacting with C=O lone pair bLP2, we want
        // to use the regular angular dependence
        if (bAngle2 || (bPlane1 && bLP2)) {
          const Coord &cAtom2_2 = pAtom2_2->GetCoords();
          double DA2 = Angle(cAtom1_1, cAtom2_1, cAtom2_2) - A2prms.R0;
          f *= f1(std::fabs(DA2), A2prms);
          LOG_F(1, "A2: {}, {}, f={}", A2prms.R0, DA2, f);
        } else if (bPlane2) {
          const Coord &cAtom2_2 = pAtom2_2->GetCoords();
          const Coord &cAtom2_3 = pAtom2_3->GetCoords();
          Plane pl2 = Plane(cAtom2_1, cAtom2_2, cAtom2_3);
          double A = std::acos(-std::fabs(Dot(v12.Unit(), pl2.VNorm()))) *
                     180.0 / M_PI;
          double DA2 = A - A2prms.R0;
          f *= f1(std::fabs(DA2), A2prms);
          LOG_F(1, "Pl2: {}, {}, f={}", A2prms.R0, DA2, f);
        } else if (bLP2) {
          const Coord &cAtom2_2 = pAtom2_2->GetCoords();
          const Coord &cAtom2_3 = pAtom2_3->GetCoords();
          const f1prms &PHI2prms = (eLP2 == InteractionCenter::LONEPAIR)
                                       ? m_PHI_lp_prms
                                       : m_PHI_plane_prms;
          Plane pl2 = Plane(cAtom2_1, cAtom2_2, cAtom2_3);
          double dPerp = DistanceFromPointToPlane(cAtom1_1, pl2);
          Coord cPerp = cAtom1_1 - dPerp * pl2.VNorm();
          LOG_F(1, "dPerp = {} check cPerp = {}", dPerp,
                DistanceFromPointToPlane(cPerp, pl2));
          double theta = std::asin(dPerp / R) * 180.0 / M_PI;
          f *= f1(std::fabs(theta), m_THETAprms);
          if (f > 0.0) {
            double phi = 180.0 - Angle(cPerp, cAtom2_1, cAtom2_2);
            double Dphi = phi - PHI2prms.R0;
            f *= f1(std::fabs(Dphi), PHI2prms);
            LOG_F(1, "LP1: {}, {}, {}, f={}", theta, PHI2prms.R0, Dphi, f);
          }
        }
        if (f > 0.0) {
          LOG_F(1, "{} {} {} {} {}", pAtom1_1->GetFullAtomName(),
                pAtom1_1->GetUser1Value(), pAtom2_1->GetFullAtomName(),
                pAtom2_1->GetUser1Value(), f);
          s += pAtom2_1->GetUser1Value() * f;
          if (bAnnotate) {
            // DM 6 Feb 2003. Include the charge and receptor density scaling
            // factors in the polar annotation score so that the score truly
            // reflects the strength of the interaction, as scored by rDock This
            // also distinguishes attractive (-ve) and repulsive (+ve)
            // annotations
            AnnotationPtr spAnnotation(new Annotation(
                pAtom1_1, pAtom2_1, R,
                f * pAtom1_1->GetUser1Value() * pAtom2_1->GetUser1Value()));
            AddAnnotation(spAnnotation);
          }
        }
      }
    }
  }
  return s;
}

// As this has a virtual base class we need a separate OwnParameterUpdated
// which can be called by concrete subclass ParameterUpdated methods
// See Stroustrup C++ 3rd edition, p395, on programming virtual base classes
void PolarSF::OwnParameterUpdated(const std::string &strName) {
  // DM 25 Oct 2000 - heavily used params
  if (strName == _R12FACTOR) {
    m_R12Factor = GetParameter(_R12FACTOR);
  } else if (strName == _R12INCR) {
    m_R12Incr = GetParameter(_R12INCR);
  } else if (strName == _DR12MIN) {
    m_DR12Min = GetParameter(_DR12MIN);
  } else if (strName == _DR12MAX) {
    m_DR12Max = GetParameter(_DR12MAX);
  } else if (strName == _A1) {
    m_A1 = GetParameter(_A1);
  } else if (strName == _DA1MIN) {
    m_DA1Min = GetParameter(_DA1MIN);
  } else if (strName == _DA1MAX) {
    m_DA1Max = GetParameter(_DA1MAX);
  } else if (strName == _A2) {
    m_A2 = GetParameter(_A2);
  } else if (strName == _DA2MIN) {
    m_DA2Min = GetParameter(_DA2MIN);
  } else if (strName == _DA2MAX) {
    m_DA2Max = GetParameter(_DA2MAX);
  } else if (strName == _ABS_DR12) {
    m_bAbsDR12 = GetParameter(_ABS_DR12);
  } else if (strName == _LP_PHI) {
    m_LP_PHI = GetParameter(_LP_PHI);
    UpdateLPprms();
  } else if (strName == _LP_DPHIMIN) {
    m_LP_DPHIMin = GetParameter(_LP_DPHIMIN);
    UpdateLPprms();
  } else if (strName == _LP_DPHIMAX) {
    m_LP_DPHIMax = GetParameter(_LP_DPHIMAX);
    UpdateLPprms();
  } else if (strName == _LP_DTHETAMIN) {
    m_LP_DTHETAMin = GetParameter(_LP_DTHETAMIN);
    UpdateLPprms();
  } else if (strName == _LP_DTHETAMAX) {
    m_LP_DTHETAMax = GetParameter(_LP_DTHETAMAX);
    UpdateLPprms();
  }
}

void PolarSF::UpdateLPprms() {
  m_PHI_lp_prms = f1prms(m_LP_PHI, m_LP_DPHIMin, m_LP_DPHIMax);
  m_PHI_plane_prms =
      f1prms(0.0, m_LP_PHI + m_LP_DPHIMin, m_LP_PHI + m_LP_DPHIMax);
  m_THETAprms = f1prms(0.0, m_LP_DTHETAMin, m_LP_DTHETAMax);
}
