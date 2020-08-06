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

#include "AromIdxSF.h"
#include "WorkSpace.h"

#include <loguru.hpp>

using namespace rxdock;

// Static data members
std::string AromIdxSF::_CT("AromIdxSF");
std::string AromIdxSF::_INCR("INCR");
std::string AromIdxSF::_R12("R12");
std::string AromIdxSF::_DR12MIN("DR12MIN");
std::string AromIdxSF::_DR12MAX("DR12MAX");
std::string AromIdxSF::_DAMIN("DAMIN");
std::string AromIdxSF::_DAMAX("DAMAX");
std::string AromIdxSF::_THRESHOLD("THRESHOLD");

AromIdxSF::AromIdxSF(const std::string &strName)
    : BaseSF(_CT, strName), m_R12(3.5), m_DR12Min(0.25), m_DR12Max(0.6),
      m_DAMin(20.0), m_DAMax(30.0), m_nArom(0), m_nGuan(0), m_threshold(0.25) {
  LOG_F(2, "AromIdxSF parameterised constructor");
  // Add parameters
  AddParameter(_INCR, 5.0);
  AddParameter(_R12, m_R12);
  AddParameter(_DR12MIN, m_DR12Min);
  AddParameter(_DR12MAX, m_DR12Max);
  AddParameter(_DAMIN, m_DAMin);
  AddParameter(_DAMAX, m_DAMax);
  AddParameter(_THRESHOLD, m_threshold);
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

AromIdxSF::~AromIdxSF() {
  LOG_F(2, "AromIdxSF destructor");
  ClearReceptor();
  ClearLigand();
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

// Override BaseSF::ScoreMap to provide additional raw descriptors
/*void AromIdxSF::ScoreMap(StringVariantMap& scoreMap) const {
  if (isEnabled()) {
    EnableAnnotations(true);
    ClearAnnotationList();
    Double rs = RawScore();
    EnableAnnotations(false);

    Double w = GetWeight();
    std::string name = GetFullName();
    Double s = w * rs;
    //DM 31 Mar 2003
    //New convention - report the raw scores and the weight in the SD file
    //rather than the weighted score and the raw score
    //Raw scores are less likely to change, whereas weights are quite volatile
    scoreMap[name] = rs;
    scoreMap[name+".w"] = w;
    scoreMap[name+".narom"] = m_nArom;
    scoreMap[name+".nguan"] = m_nGuan;
    //scoreMap[name+".ss"] = m_ss;
    //scoreMap[name+".sp"] = m_sp;
    //scoreMap[name+".pp"] = m_pp;
    std::vector<std::string> annList;
    RenderAnnotationList(GetName(),annList);
    scoreMap[AnnotationHandler::_ANNOTATION_FIELD] += annList;
    ClearAnnotationList();
    return s;
  }
  else {
    return 0.0;
  }
}*/

void AromIdxSF::SetupReceptor() {
  ClearReceptor();
  if (GetReceptor().Null())
    return;

  m_spAromGrid = CreateInteractionGrid();
  m_spGuanGrid = CreateInteractionGrid();

  double idxIncr = GetParameter(_INCR); // vdw Radius increment for indexing
  double maxError = GetMaxError();
  idxIncr += maxError;

  // First, index the receptor aromatic rings
  // Limit to those within RANGE of the docking site
  // Have to filter manually as GetReceptorAtomsInRange not appropriate for ring
  // centers
  AtomListList recepRingLists = GetReceptor()->GetRingAtomLists();
  DockingSitePtr spDS = GetWorkSpace()->GetDockingSite();

  int nCoords = GetReceptor()->GetNumSavedCoords() - 1;
  if (nCoords > 0) {
    for (AtomListListConstIter rIter = recepRingLists.begin();
         rIter != recepRingLists.end(); rIter++) {
      if (GetNumAtomsWithPredicate(*rIter, isPiAtom()) == (*rIter).size()) {
        PseudoAtomPtr spPseudoAtom = GetReceptor()->AddPseudoAtom(*rIter);
        AtomList recepPiAtoms = spPseudoAtom->GetAtomList();
        AtomPtr spRecepPiAtom1 = recepPiAtoms[0];
        AtomPtr spRecepPiAtom2 = recepPiAtoms[1];
        InteractionCenter *pIntnCenter(new InteractionCenter(
            spPseudoAtom, spRecepPiAtom1, spRecepPiAtom2));
        m_recepAromList.push_back(pIntnCenter);
      }
    }
    AtomList recepGuanList = GetAtomListWithPredicate(
        GetReceptor()->GetAtomList(), isAtomGuanidiniumCarbon());
    for (AtomListConstIter iter = recepGuanList.begin();
         iter != recepGuanList.end(); iter++) {
      AtomList recepPiAtoms = GetBondedAtomList(*iter);
      AtomPtr spRecepPiAtom1 = recepPiAtoms[0];
      AtomPtr spRecepPiAtom2 = recepPiAtoms[1];
      InteractionCenter *pIntnCenter(
          new InteractionCenter(*iter, spRecepPiAtom1, spRecepPiAtom2));
      m_recepGuanList.push_back(pIntnCenter); // Store the interaction center
    }

    for (int i = 1; i <= nCoords; i++) {
      LOG_F(1, "AromIdxSF setup indexing receptor coords #{}", i);
      GetReceptor()->RevertCoords(i);
      for (InteractionCenterListConstIter iter = m_recepAromList.begin();
           iter != m_recepAromList.end(); iter++) {
        m_spAromGrid->SetInteractionLists(*iter, idxIncr);
      }
      m_spAromGrid->UniqueInteractionLists();
      for (InteractionCenterListConstIter iter = m_recepGuanList.begin();
           iter != m_recepGuanList.end(); iter++) {
        m_spGuanGrid->SetInteractionLists(*iter, idxIncr);
      }
      m_spGuanGrid->UniqueInteractionLists();
    }
  } else {
    DockingSite::isAtomInRange bIsInRange(spDS->GetGrid(), 0.0,
                                          GetCorrectedRange());
    for (AtomListListConstIter rIter = recepRingLists.begin();
         rIter != recepRingLists.end(); rIter++) {
      // Check that all ring atoms are pi-atoms (crude test for aromaticity)
      if (GetNumAtomsWithPredicate(*rIter, isPiAtom()) == (*rIter).size()) {
        PseudoAtomPtr spPseudoAtom = GetReceptor()->AddPseudoAtom(*rIter);
        if (bIsInRange(spPseudoAtom)) {
          AtomList recepPiAtoms = spPseudoAtom->GetAtomList();
          AtomPtr spRecepPiAtom1 = recepPiAtoms[0];
          AtomPtr spRecepPiAtom2 = recepPiAtoms[1];
          InteractionCenter *pIntnCenter(new InteractionCenter(
              spPseudoAtom, spRecepPiAtom1, spRecepPiAtom2));
          m_recepAromList.push_back(pIntnCenter);
          m_spAromGrid->SetInteractionLists(pIntnCenter, idxIncr);
        }
      }
    }

    // Now the guanidinium groups (positively charged carbons)
    // Limit to those within RANGE of the docking site
    AtomList atomList = spDS->GetAtomList(GetReceptor()->GetAtomList(), 0.0,
                                          GetCorrectedRange());
    AtomList recepGuanList =
        GetAtomListWithPredicate(atomList, isAtomGuanidiniumCarbon());

    for (AtomListConstIter iter = recepGuanList.begin();
         iter != recepGuanList.end(); iter++) {
      AtomList recepPiAtoms = GetBondedAtomList(*iter);
      AtomPtr spRecepPiAtom1 = recepPiAtoms[0];
      AtomPtr spRecepPiAtom2 = recepPiAtoms[1];
      InteractionCenter *pIntnCenter(
          new InteractionCenter(*iter, spRecepPiAtom1, spRecepPiAtom2));
      m_recepGuanList.push_back(pIntnCenter); // Store the interaction center
      m_spGuanGrid->SetInteractionLists(pIntnCenter, idxIncr);
    }
  }
}

void AromIdxSF::SetupLigand() {
  ClearLigand();
  if (GetLigand().Null())
    return;

  // First, create the ligand aromatic interaction centers
  AtomListList ligRingLists = GetLigand()->GetRingAtomLists();
  for (AtomListListConstIter rIter = ligRingLists.begin();
       rIter != ligRingLists.end(); rIter++) {
    // Check that all ring atoms are pi-atoms (crude test for aromaticity)
    if (GetNumAtomsWithPredicate(*rIter, isPiAtom()) == (*rIter).size()) {
      PseudoAtomPtr spPseudoAtom = GetLigand()->AddPseudoAtom(*rIter);
      // AtomList ligPiAtoms = GetBondedAtomList(spPseudoAtom);
      AtomList ligPiAtoms = spPseudoAtom->GetAtomList();
      AtomPtr spLigPiAtom1 = ligPiAtoms[0];
      AtomPtr spLigPiAtom2 = ligPiAtoms[1];
      InteractionCenter *pIntnCenter(
          new InteractionCenter(spPseudoAtom, spLigPiAtom1, spLigPiAtom2));
      m_ligAromList.push_back(pIntnCenter);
    }
  }

  // Now the guanidinium groups (positively charged carbons)
  AtomList ligGuanList = GetAtomListWithPredicate(GetLigand()->GetAtomList(),
                                                  isAtomGuanidiniumCarbon());
  for (AtomListConstIter iter = ligGuanList.begin(); iter != ligGuanList.end();
       iter++) {
    AtomList ligPiAtoms = GetBondedAtomList(*iter);
    AtomPtr spLigPiAtom1 = ligPiAtoms[0];
    AtomPtr spLigPiAtom2 = ligPiAtoms[1];
    InteractionCenter *pIntnCenter(
        new InteractionCenter(*iter, spLigPiAtom1, spLigPiAtom2));
    m_ligGuanList.push_back(pIntnCenter); // Store the interaction center
  }
}

void AromIdxSF::SetupScore() {
  // No further setup required
}

double AromIdxSF::RawScore() const {
  double score = 0.0;
  m_nArom = 0;
  m_nGuan = 0;

  // Check grids are defined
  if (m_spAromGrid.Null() || m_spGuanGrid.Null())
    return score;

  f1prms Rprms = GetRprms(); // Distance params
  f1prms Aprms = GetAprms(); // Donor angle params

  // Loop over ligand aromatic interaction centers
  for (InteractionCenterListConstIter ligIter = m_ligAromList.begin();
       ligIter != m_ligAromList.end(); ligIter++) {
    const Coord &cLig1 = (*ligIter)->GetAtom1Ptr()->GetCoords();
    // Get the list of nearby receptor aromatic centers
    const InteractionCenterList &recepAromList =
        m_spAromGrid->GetInteractionList(cLig1);
    // Get the list of nearby receptor guanidinium centers
    const InteractionCenterList &recepGuanList =
        m_spGuanGrid->GetInteractionList(cLig1);

    double s = AromScore(*ligIter, recepAromList, Rprms, Aprms);
    s += AromScore(*ligIter, recepGuanList, Rprms, Aprms);
    // Double s = PiScore(*ligIter,recepAromList);
    if (s > m_threshold) {
      m_nArom++;
    }
    score += s;
  }

  // Loop over ligand guanidinium interaction centers
  for (InteractionCenterListConstIter ligIter = m_ligGuanList.begin();
       ligIter != m_ligGuanList.end(); ligIter++) {
    const Coord &cLig1 = (*ligIter)->GetAtom1Ptr()->GetCoords();
    // Get the list of nearby receptor aromatic centers
    const InteractionCenterList &recepAromList =
        m_spAromGrid->GetInteractionList(cLig1);

    double s = AromScore(*ligIter, recepAromList, Rprms, Aprms);
    // Double s = 0.0;
    if (s > m_threshold) {
      m_nGuan++;
    }
    score += s;
  }
  return score;
}

// Clear the receptor and ligand grids and lists respectively
// As we are not using smart pointers, there is some memory management to do
void AromIdxSF::ClearReceptor() {
  // Wipe the grids
  m_spAromGrid = InteractionGridPtr();
  m_spGuanGrid = InteractionGridPtr();
  // Delete the receptor interaction centers
  for (InteractionCenterListIter iter = m_recepAromList.begin();
       iter != m_recepAromList.end(); iter++) {
    delete *iter;
  }
  m_recepAromList.clear();
  for (InteractionCenterListIter iter = m_recepGuanList.begin();
       iter != m_recepGuanList.end(); iter++) {
    delete *iter;
  }
  m_recepGuanList.clear();
}
void AromIdxSF::ClearLigand() {
  // Delete the ligand interaction centers
  for (InteractionCenterListIter iter = m_ligAromList.begin();
       iter != m_ligAromList.end(); iter++) {
    delete *iter;
  }
  m_ligAromList.clear();
  for (InteractionCenterListIter iter = m_ligGuanList.begin();
       iter != m_ligGuanList.end(); iter++) {
    delete *iter;
  }
  m_ligGuanList.clear();
}

// DM 25 Oct 2000 - track changes to parameter values in local data members
// ParameterUpdated is invoked by ParamHandler::SetParameter
void AromIdxSF::ParameterUpdated(const std::string &strName) {
  if (strName == _R12) {
    m_R12 = GetParameter(_R12);
  } else if (strName == _DR12MIN) {
    m_DR12Min = GetParameter(_DR12MIN);
  } else if (strName == _DR12MAX) {
    m_DR12Max = GetParameter(_DR12MAX);
  } else if (strName == _DAMIN) {
    m_DAMin = GetParameter(_DAMIN);
  } else if (strName == _DAMAX) {
    m_DAMax = GetParameter(_DAMAX);
  } else if (strName == _THRESHOLD) {
    m_threshold = GetParameter(_THRESHOLD);
  } else {
    BaseIdxSF::OwnParameterUpdated(strName);
    BaseSF::ParameterUpdated(strName);
  }
}

// The actual aromatic score, between a given interaction center and a list of
// near neighbour centers
double AromIdxSF::AromScore(const InteractionCenter *pIC1,
                            const InteractionCenterList &IC2List,
                            const f1prms &Rprms, const f1prms &Aprms) const {
  double s(0.0);
  if (IC2List.empty()) {
    return s;
  }

  bool bAnnotate = isAnnotationEnabled();

  Atom *pAtom1_1 = pIC1->GetAtom1Ptr();
  const Coord &cAtom1_1 = pAtom1_1->GetCoords();
  Atom *pAtom1_2 = pIC1->GetAtom2Ptr();
  Atom *pAtom1_3 = pIC1->GetAtom3Ptr();
  Plane pl1 = Plane(cAtom1_1, pAtom1_2->GetCoords(), pAtom1_3->GetCoords());

  for (InteractionCenterListConstIter IC2Iter = IC2List.begin();
       IC2Iter != IC2List.end(); IC2Iter++) {
    Atom *pAtom2_1 = (*IC2Iter)->GetAtom1Ptr();
    const Coord &cAtom2_1 = pAtom2_1->GetCoords();
    Atom *pAtom2_2 = (*IC2Iter)->GetAtom2Ptr();
    Atom *pAtom2_3 = (*IC2Iter)->GetAtom3Ptr();
    Plane pl2 = Plane(cAtom2_1, pAtom2_2->GetCoords(), pAtom2_3->GetCoords());

    // Calculate average perp. distance from each pi center to the other ring
    // plane
    double R = (std::fabs(DistanceFromPointToPlane(cAtom1_1, pl2)) +
                std::fabs(DistanceFromPointToPlane(cAtom2_1, pl1))) /
               2.0;
    double f = f1(R - Rprms.R0, Rprms);
    // Only calculate average slip angle if f  > 0
    if (f > 0.0) {
      Vector v = cAtom2_1 - cAtom1_1;
      double sa = (std::acos(std::fabs(Dot(v.Unit(), pl1.VNorm()))) +
                   std::acos(std::fabs(Dot(v.Unit(), pl2.VNorm())))) *
                  90.0 / M_PI;
      f *= f1(sa, Aprms);
      if (f > 0.0) {
        s += f;
        if (bAnnotate) {
          // Convention is that favourable annotations are negative, so invert
          // the raw score here
          AnnotationPtr spAnnotation(new Annotation(pAtom1_1, pAtom2_1, R, -f));
          AddAnnotation(spAnnotation);
        }
      }
    }
  }
  return s;
}

// The actual aromatic score, between a given interaction center and a list of
// near neighbour centers
double AromIdxSF::PiScore(const InteractionCenter *pIC1,
                          const InteractionCenterList &IC2List) const {
  m_ss = 0.0; // sigma-sigma
  m_pp = 0.0; // pi-pi
  m_sp = 0.0; // sigma-pi

  if (IC2List.empty()) {
    return 0.0;
  }

  double delta = 0.47;

  Atom *pAtom1_1 = pIC1->GetAtom1Ptr();
  Atom *pAtom1_2 = pIC1->GetAtom2Ptr();
  Atom *pAtom1_3 = pIC1->GetAtom3Ptr();
  Plane pl1 = Plane(pAtom1_1->GetCoords(), pAtom1_2->GetCoords(),
                    pAtom1_3->GetCoords());
  Vector v1 = delta * pl1.VNorm();
  // First calculate and store the sigma and pi-system coords (above and below
  // plane)
  PseudoAtom *pseudo1 = dynamic_cast<PseudoAtom *>(pAtom1_1);
  CoordList sigma1;
  CoordList pi1;
  if (pseudo1) {
    AtomList al1 = pseudo1->GetAtomList();
    for (AtomListConstIter iter = al1.begin(); iter != al1.end(); iter++) {
      sigma1.push_back((*iter)->GetCoords());
      pi1.push_back((*iter)->GetCoords() + v1);
      pi1.push_back((*iter)->GetCoords() - v1);
    }
  }

  for (InteractionCenterListConstIter IC2Iter = IC2List.begin();
       IC2Iter != IC2List.end(); IC2Iter++) {
    Atom *pAtom2_1 = (*IC2Iter)->GetAtom1Ptr();
    Atom *pAtom2_2 = (*IC2Iter)->GetAtom2Ptr();
    Atom *pAtom2_3 = (*IC2Iter)->GetAtom3Ptr();
    Plane pl2 = Plane(pAtom2_1->GetCoords(), pAtom2_2->GetCoords(),
                      pAtom2_3->GetCoords());
    Vector v2 = delta * pl2.VNorm();
    // First calculate and store the sigma and pi-system coords (above and below
    // plane)
    PseudoAtom *pseudo2 = dynamic_cast<PseudoAtom *>(pAtom2_1);
    CoordList sigma2;
    CoordList pi2;
    if (pseudo2) {
      AtomList al2 = pseudo2->GetAtomList();
      for (AtomListConstIter iter = al2.begin(); iter != al2.end(); iter++) {
        sigma2.push_back((*iter)->GetCoords());
        pi2.push_back((*iter)->GetCoords() + v2);
        pi2.push_back((*iter)->GetCoords() - v2);
      }
    }

    for (CoordListConstIter cIter1 = sigma1.begin(); cIter1 != sigma1.end();
         cIter1++) {
      for (CoordListConstIter cIter2 = sigma2.begin(); cIter2 != sigma2.end();
           cIter2++) {
        double r = Length(*cIter1, *cIter2) + 0.0001; // prevent div by zero
        m_ss += (1.0 / r);
      }
      for (CoordListConstIter cIter2 = pi2.begin(); cIter2 != pi2.end();
           cIter2++) {
        double r = Length(*cIter1, *cIter2) + 0.0001; // prevent div by zero
        m_sp += (1.0 / r);
      }
    }
    for (CoordListConstIter cIter1 = pi1.begin(); cIter1 != pi1.end();
         cIter1++) {
      for (CoordListConstIter cIter2 = pi2.begin(); cIter2 != pi2.end();
           cIter2++) {
        double r = Length(*cIter1, *cIter2) + 0.0001; // prevent div by zero
        m_pp += (1.0 / r);
      }
      for (CoordListConstIter cIter2 = sigma2.begin(); cIter2 != sigma2.end();
           cIter2++) {
        double r = Length(*cIter1, *cIter2) + 0.0001; // prevent div by zero
        m_sp += (1.0 / r);
      }
    }
  }
  // Correct the totals for the absolute charges (sigma = +1; pi = -0.5) and
  // convert to kcal
  m_ss *= 332.17;
  m_sp *= (-0.5 * 332.17);
  m_pp *= (0.25 * 332.17);
  return m_ss + m_sp + m_pp;
}
