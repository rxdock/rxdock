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

#include "VdwSF.h"
#include "Model.h"

#include <loguru.hpp>

using namespace rxdock;

// Static data members
std::string VdwSF::_CT("VdwSF");
std::string VdwSF::_USE_4_8("USE_4_8");
std::string VdwSF::_USE_TRIPOS("USE_TRIPOS");
std::string VdwSF::_RMAX("RMAX");
std::string VdwSF::_ECUT("ECUT");
std::string VdwSF::_E0("E0");

std::string &VdwSF::GetEcut() { return _ECUT; }

VdwSF::VdwSF()
    : m_use_4_8(true), m_use_tripos(false), m_rmax(1.5), m_ecut(1.0),
      m_e0(1.5) {
  LOG_F(2, "VdwSF default constructor");
  // Add parameters
  AddParameter(_USE_4_8, m_use_4_8);
  AddParameter(_USE_TRIPOS, m_use_tripos);
  AddParameter(_RMAX, m_rmax);
  AddParameter(_ECUT, m_ecut);
  AddParameter(_E0, m_e0);
  m_spVdwSource = ParameterFileSourcePtr(
      new ParameterFileSource(GetDataFileName("data/sf", "Tripos52_vdw.prm")));
  Setup();
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

VdwSF::~VdwSF() {
  LOG_F(2, "VdwSF default constructor");
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

// As this has a virtual base class we need a separate OwnParameterUpdated
// which can be called by concrete subclass ParameterUpdated methods
// See Stroustrup C++ 3rd edition, p395, on programming virtual base classes
void VdwSF::OwnParameterUpdated(const std::string &strName) {
  // DM 25 Oct 2000 - heavily used params
  if (strName == _USE_4_8) {
    m_use_4_8 = GetParameter(_USE_4_8);
    Setup();
  } else if (strName == _USE_TRIPOS) {
    m_use_tripos = GetParameter(_USE_TRIPOS);
    Setup();
  } else if (strName == _RMAX) {
    m_rmax = GetParameter(_RMAX);
    Setup();
  } else if (strName == _ECUT) {
    m_ecut = GetParameter(_ECUT);
    SetupCloseRange();
  } else if (strName == _E0) {
    m_e0 = GetParameter(_E0);
    SetupCloseRange();
  }
}

// Used by subclasses to calculate vdW potential between pAtom and all atoms in
// atomList
double VdwSF::VdwScore(const Atom *pAtom, const AtomRList &atomList) const {
  double score = 0.0;
  if (atomList.empty()) {
    return score;
  }

  const Coord &c1 = pAtom->GetCoords();
  // Get the iterator into the appropriate row of the vdw table for this atom
  // type
  TriposAtomType::eType type1 = pAtom->GetTriposType();
  VdwTableConstIter iter1 = m_vdwTable.begin() + type1;

  // 4-8 potential, never annotated
  if (m_use_4_8) {
    for (AtomRListConstIter iter = atomList.begin(); iter != atomList.end();
         iter++) {
      const Coord &c2 = (*iter)->GetCoords();
      double R_sq = Length2(c1, c2); // Distance squared
      TriposAtomType::eType type2 = (*iter)->GetTriposType();
      // iter2 points to the vdw params for this atom type pair
      VdwRowConstIter iter2 = (*iter1).begin() + type2;
      double s = f4_8(R_sq, *iter2);
      // XB NOTE: Apply weight here
      // double wxb = (*iter)->GetReweight();
      // s *=  wxb;
      // LOG_F(1, "4-8vdw {}  weight: {} score {} score*w {}",
      //       (*iter)->GetName(), (*iter)->GetReweight(), s,
      //       s * (*iter)->GetReweight());
      // XB END MODIFICATIONS
      score += s;
    }
  }
  // 6-12 with annotation
  else if (isAnnotationEnabled()) {
    for (AtomRListConstIter iter = atomList.begin(); iter != atomList.end();
         iter++) {
      const Coord &c2 = (*iter)->GetCoords();
      double R_sq = Length2(c1, c2); // Distance squared
      TriposAtomType::eType type2 = (*iter)->GetTriposType();
      // iter2 points to the vdw params for this atom type pair
      VdwRowConstIter iter2 = (*iter1).begin() + type2;
      double s = f6_12(R_sq, *iter2);
      // XB NOTE: Apply weight here
      // double wxb = (*iter)->GetReweight();
      // s *=  wxb;
      // LOG_F(1, "6-12vdw {}  weight: {} score {} score*w {}",
      //       (*iter)->GetName(), (*iter)->GetReweight(), s,
      //       s * (*iter)->GetReweight());
      // XB END MODIFICATIONS
      if (s != 0.0) {
        score += s;
        AnnotationPtr spAnnotation(
            new Annotation(pAtom, *iter, std::sqrt(R_sq), s));
        AddAnnotation(spAnnotation);
      }
    }
  }
  // 6-12 without annotation
  else {
    for (AtomRListConstIter iter = atomList.begin(); iter != atomList.end();
         iter++) {
      const Coord &c2 = (*iter)->GetCoords();
      double R_sq = Length2(c1, c2); // Distance squared
      TriposAtomType::eType type2 = (*iter)->GetTriposType();
      // iter2 points to the vdw params for this atom type pair
      VdwRowConstIter iter2 = (*iter1).begin() + type2;
      double s = f6_12(R_sq, *iter2);
      // XB NOTE: Apply weight here
      // double wxb = (*iter)->GetReweight();
      // s *=  wxb;
      // LOG_F(1, "6-12vdw {}  weight: {} score {} score*w {}",
      //       (*iter)->GetName(), (*iter)->GetReweight(), s,
      //       s * (*iter)->GetReweight());
      // XB END MODIFICATIONS
      score += s;
    }
  }
  return score;
}

// As above, but score is calculated only between enabled atoms
double VdwSF::VdwScoreEnabledOnly(const Atom *pAtom,
                                  const AtomRList &atomList) const {
  double score = 0.0;
  if (!pAtom->GetEnabled() || atomList.empty()) {
    return score;
  }

  const Coord &c1 = pAtom->GetCoords();
  // Get the iterator into the appropriate row of the vdw table for this atom
  // type
  TriposAtomType::eType type1 = pAtom->GetTriposType();
  VdwTableConstIter iter1 = m_vdwTable.begin() + type1;

  // 4-8 potential, never annotated
  if (m_use_4_8) {
    for (AtomRListConstIter iter = atomList.begin(); iter != atomList.end();
         iter++) {
      if ((*iter)->GetEnabled()) {
        const Coord &c2 = (*iter)->GetCoords();
        double R_sq = Length2(c1, c2); // Distance squared
        TriposAtomType::eType type2 = (*iter)->GetTriposType();
        // iter2 points to the vdw params for this atom type pair
        VdwRowConstIter iter2 = (*iter1).begin() + type2;
        double s = f4_8(R_sq, *iter2);
        score += s;
      }
    }
  }
  // 6-12 with annotation
  else if (isAnnotationEnabled()) {
    for (AtomRListConstIter iter = atomList.begin(); iter != atomList.end();
         iter++) {
      if ((*iter)->GetEnabled()) {
        const Coord &c2 = (*iter)->GetCoords();
        double R_sq = Length2(c1, c2); // Distance squared
        TriposAtomType::eType type2 = (*iter)->GetTriposType();
        // iter2 points to the vdw params for this atom type pair
        VdwRowConstIter iter2 = (*iter1).begin() + type2;
        double s = f6_12(R_sq, *iter2);
        if (s != 0.0) {
          score += s;
          AnnotationPtr spAnnotation(
              new Annotation(pAtom, *iter, std::sqrt(R_sq), s));
          AddAnnotation(spAnnotation);
        }
      }
    }
  }
  // 6-12 without annotation
  else {
    for (AtomRListConstIter iter = atomList.begin(); iter != atomList.end();
         iter++) {
      if ((*iter)->GetEnabled()) {
        const Coord &c2 = (*iter)->GetCoords();
        double R_sq = Length2(c1, c2); // Distance squared
        TriposAtomType::eType type2 = (*iter)->GetTriposType();
        // iter2 points to the vdw params for this atom type pair
        VdwRowConstIter iter2 = (*iter1).begin() + type2;
        double s = f6_12(R_sq, *iter2);
        score += s;
      }
    }
  }
  return score;
}

// XB This is the old  VdwScore, without reweighting factors
// Double VdwSF::VdwScoreIntra(const Atom* pAtom, const AtomRList&
// atomList) const {
//  Double score = 0.0;
//  if (atomList.empty()) {
//    return score;
//  }
//
//  const Coord& c1 = pAtom->GetCoords();
//  //Get the iterator into the appropriate row of the vdw table for this atom
//  type TriposAtomType::eType type1 = pAtom->GetTriposType();
//  VdwTableConstIter iter1 = m_vdwTable.begin() + type1;
//
//  //4-8 potential, never annotated
//  if (m_use_4_8) {
//    for (AtomRListConstIter iter = atomList.begin(); iter !=
//    atomList.end(); iter++) {
//      const Coord& c2 = (*iter)->GetCoords();
//      Double R_sq = Length2(c1,c2);//Distance squared
//      TriposAtomType::eType type2 = (*iter)->GetTriposType();
//      //iter2 points to the vdw params for this atom type pair
//      VdwRowConstIter iter2 = (*iter1).begin() + type2;
//      Double s = f4_8(R_sq,*iter2);
//      score += s;
//    }
//  }
//  //6-12 with annotation
//  else if (isAnnotationEnabled()) {
//    for (AtomRListConstIter iter = atomList.begin(); iter !=
//    atomList.end(); iter++) {
//      const Coord& c2 = (*iter)->GetCoords();
//      Double R_sq = Length2(c1,c2);//Distance squared
//      TriposAtomType::eType type2 = (*iter)->GetTriposType();
//      //iter2 points to the vdw params for this atom type pair
//      VdwRowConstIter iter2 = (*iter1).begin() + type2;
//      Double s = f6_12(R_sq,*iter2);
//      if (s != 0.0) {
//	score += s;
//	AnnotationPtr spAnnotation(new
// Annotation(pAtom,*iter,std::sqrt(R_sq),s));
// AddAnnotation(spAnnotation);
//      }
//    }
//  }
//  //6-12 without annotation
//  else {
//    for (AtomRListConstIter iter = atomList.begin(); iter !=
//    atomList.end(); iter++) {
//      const Coord& c2 = (*iter)->GetCoords();
//      Double R_sq = Length2(c1,c2);//Distance squared
//      TriposAtomType::eType type2 = (*iter)->GetTriposType();
//      //iter2 points to the vdw params for this atom type pair
//      VdwRowConstIter iter2 = (*iter1).begin() + type2;
//      Double s = f6_12(R_sq,*iter2);
//     score += s;
//    }
//  }
//  return score;
//}
double VdwSF::MaxVdwRange(const Atom *pAtom) const {
  return m_maxRange[pAtom->GetTriposType()];
}

double VdwSF::MaxVdwRange(TriposAtomType::eType t) const {
  return m_maxRange[t];
}

// Initialise m_vdwTable with appropriate params for each atom type pair
void VdwSF::Setup() {
  LOG_F(2, "VdwSF::Setup");
  TriposAtomType triposType;
  LOG_F(1, "TYPE1,TYPE2,A,B,RMIN,KIJ,RMA");
  // Dummy read to force parsing of file, otherwise the first SetSection is
  // overridden
  std::vector<std::string> secList = m_spVdwSource->GetSectionList();
  std::string _R("R");
  std::string _K("K");
  std::string _IP("IP");
  std::string _POL("POL");
  std::string _ISHBD("isHBD");
  std::string _ISHBA("isHBA");
  m_vdwTable =
      VdwTable(TriposAtomType::MAXTYPES, VdwRow(TriposAtomType::MAXTYPES));
  m_maxRange = std::vector<double>(TriposAtomType::MAXTYPES, 0.0);
  for (int i = TriposAtomType::UNDEFINED; i < TriposAtomType::MAXTYPES; i++) {
    // Read the params for atom type i
    std::string stri = triposType.Type2Str(TriposAtomType::eType(i));
    m_spVdwSource->SetSection(stri);
    double Ri = m_spVdwSource->GetParameterValue(_R); // vdw radius
    double Ki = m_spVdwSource->GetParameterValue(_K); // Tripos 5.2 well
                                                      // depth
    bool hasIPi = m_spVdwSource->isParameterPresent(_IP);
    bool isHBDi = m_spVdwSource->isParameterPresent(_ISHBD);
    bool isHBAi = m_spVdwSource->isParameterPresent(_ISHBA);
    double Ii, alphai;
    if (hasIPi) {
      Ii = m_spVdwSource->GetParameterValue(_IP); // Ionisation potential, I
      alphai = m_spVdwSource->GetParameterValue(_POL); // Polarisability, alpha
    }
    // Symmetric matrix, so only need to generate one half (i.e.
    // m_vdwTable[i][j] = m_vdwTable[j][i])
    for (int j = i; j < TriposAtomType::MAXTYPES; j++) {
      // Read the params for atom type j
      std::string strj = triposType.Type2Str(TriposAtomType::eType(j));
      m_spVdwSource->SetSection(strj);
      double Rj = m_spVdwSource->GetParameterValue(_R); // vdw radius
      double Kj = m_spVdwSource->GetParameterValue(_K); // Tripos 5.2 well depth
      bool hasIPj = m_spVdwSource->isParameterPresent(_IP);
      bool isHBDj = m_spVdwSource->isParameterPresent(_ISHBD);
      bool isHBAj = m_spVdwSource->isParameterPresent(_ISHBA);
      double Ij, alphaj;
      if (hasIPj) {
        Ij = m_spVdwSource->GetParameterValue(_IP); // Ionisation potential, I
        alphaj = m_spVdwSource->GetParameterValue(_POL); // Polarisability,
                                                         // alpha
      }

      // Create the vdw params for this atom type pair
      vdwprms prms;
      prms.rmin = Ri + Rj; // rmin is the sum of Ri and Rj
      double rmax = prms.rmin * m_rmax;
      m_maxRange[i] = std::max(m_maxRange[i],
                               rmax); // Keep track of max range for atom type i
      m_maxRange[j] = std::max(m_maxRange[j],
                               rmax); // Keep track of max range for atom type j
      prms.rmax_sq = rmax * rmax;     // Max range**2
      // EITHER: Well depth is zero between donor Hs and acceptors
      if ((isHBDi && isHBAj) || (isHBDj && isHBAi)) {
        prms.kij = 0.0;
      }
      // OR, use Standard L-J combination rules for well depth
      // switch to this mode if either atom is missing IP values
      else if (m_use_tripos || !hasIPi || !hasIPj) {
        prms.kij = std::sqrt(Ki * Kj);
      }
      // OR, use GOLD rules based on ionisation potential and polarisability
      else {
        double D = 0.345 * Ii * Ij * alphai * alphaj / (Ii + Ij);
        double C = 0.5 * D * std::pow(prms.rmin, 6);
        prms.kij = D * D / (4.0 * C);
      }
      // Having got the well depth, we can now generate A and B for either 4-8
      // or 6-12
      double rmin_pwr =
          (m_use_4_8) ? std::pow(prms.rmin, 4) : std::pow(prms.rmin, 6);
      prms.A = prms.kij * rmin_pwr * rmin_pwr;
      prms.B = 2.0 * prms.kij * rmin_pwr;
      m_vdwTable[i][j] = prms;
      m_vdwTable[j][i] = prms;
      LOG_F(1, "{},{},{},{},{},{},{}",
            triposType.Type2Str(TriposAtomType::eType(i)),
            triposType.Type2Str(TriposAtomType::eType(j)), prms.A, prms.B,
            prms.rmin, prms.kij, std::sqrt(prms.rmax_sq));
    }
  }
  // Now we can regenerate the close range params
  SetupCloseRange();
  // Set the overall range automatically from the max range for each atom type
  double range = *(std::max_element(m_maxRange.begin(), m_maxRange.end()));
  SetRange(range);
  LOG_F(INFO, "Overall max range for scoring function = {}", range);
}

// Regenerate the short-range params only (called more frequently)
void VdwSF::SetupCloseRange() {
  LOG_F(2, "VdwSF::SetupCloseRange()");
  LOG_F(1, "TYPE1,TYPE2,RCUT,ECUT,SLOPE,E0");
  TriposAtomType triposType;
  double x = 1.0 + std::sqrt(1.0 + m_ecut);
  double p = (m_use_4_8) ? std::pow(x, 1.0 / 4.0) : std::pow(x, 1.0 / 6.0);
  double c = 1.0 / p;
  for (VdwTableIter iter1 = m_vdwTable.begin(); iter1 != m_vdwTable.end();
       iter1++) {
    for (VdwRowIter iter2 = (*iter1).begin(); iter2 != (*iter1).end();
         iter2++) {
      (*iter2).rcutoff_sq = std::pow((*iter2).rmin * c, 2);
      (*iter2).ecutoff = (*iter2).kij * m_ecut;
      (*iter2).e0 = (*iter2).ecutoff * m_e0;
      (*iter2).slope = ((*iter2).e0 - (*iter2).ecutoff) / (*iter2).rcutoff_sq;
      LOG_F(
          1, "{},{},{},{},{},{}",
          triposType.Type2Str(
              TriposAtomType::eType(iter1 - m_vdwTable.begin())),
          triposType.Type2Str(TriposAtomType::eType(iter2 - (*iter1).begin())),
          std::sqrt((*iter2).rcutoff_sq), (*iter2).ecutoff, (*iter2).slope,
          (*iter2).e0);
    }
  }
}

// Index the flexible interactions between two atom lists.
// Foreach atom in the first list, compile a list of the atoms in the second
// list whose distance can vary to that atom intns container should have been
// initialised to a vector of empty atom lists prior to calling BuildIntraMap
void VdwSF::BuildIntraMap(const AtomRList &atomList1,
                          const AtomRList &atomList2,
                          AtomRListList &intns) const {
  SelectAtom selectAtom(true);
  isAtomSelected isSelected;
  bool bSingleList = atomList2.empty(); // If true, then index the flexible
                                        // interactions within atomList1
  for (AtomRListConstIter iter = atomList1.begin(); iter != atomList1.end();
       iter++) {
    // 1. First select all atoms in the list
    std::for_each(atomList1.begin(), atomList1.end(), selectAtom);
    Atom *pAtom = (*iter);
    unsigned int id = pAtom->GetAtomId() - 1;
    Assert<Assertion>(id < intns.size());
    Model *pModel = pAtom->GetModelPtr();
    // 2. Now deselect all atoms in the same model as the current atom
    pModel->SetAtomSelectionFlags(false);
    // 3. Now select all the flexible interactions to this atom within the
    // current model (Atoms in different models are not affected, therefore
    // remain selected from step 1) This way, we deal correctly with situations
    // where the atoms in atomList may belong to different models
    pModel->SelectFlexAtoms(pAtom);
    if (bSingleList) {
      std::copy_if(iter + 1, atomList1.end(), std::back_inserter(intns[id]),
                   isSelected);
    } else {
      std::copy_if(atomList2.begin(), atomList2.end(),
                   std::back_inserter(intns[id]), isSelected);
    }
    LOG_F(1, "VdwSF::BuildIntraMap: {}; N = {}", pAtom->GetFullAtomName(),
          intns[id].size());
    LOG_IF_F(1, !intns[id].empty(), "{} to {}",
             intns[id].front()->GetFullAtomName(),
             intns[id].back()->GetFullAtomName());
  }
}

// Convenience method for above.
// Indexes the flexible interactions within a single atom list
void VdwSF::BuildIntraMap(const AtomRList &atomList,
                          AtomRListList &intns) const {
  AtomRList emptyList;
  BuildIntraMap(atomList, emptyList, intns);
}

void VdwSF::Partition(const AtomRList &atomList, const AtomRListList &intns,
                      AtomRListList &prtIntns, double dist) const {
  // Clear the existing partitioned lists
  for (AtomRListListIter iter = prtIntns.begin(); iter != prtIntns.end();
       iter++)
    (*iter).clear();
  LOG_F(1, "VdwSF::Partition: dist={}", dist);

  for (AtomRListConstIter iter = atomList.begin(); iter != atomList.end();
       iter++) {
    int id = (*iter)->GetAtomId() - 1;
    if (dist > 0.0) {
      isD_lt bInRange(*iter, dist);
      // This copies all interactions which are within dist A of atom
      std::copy_if(intns[id].begin(), intns[id].end(),
                   std::back_inserter(prtIntns[id]), bInRange);
    } else {
      // Remove partition - simply copy from the master list of vdW interactions
      prtIntns[id] = intns[id];
    }
    LOG_F(1, "VdwSF::Partition: {}; $vdw={}: #prtn={}",
          (*iter)->GetFullAtomName(), intns[id].size(), prtIntns[id].size());
  }
}
