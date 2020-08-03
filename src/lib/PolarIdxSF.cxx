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

#include "PolarIdxSF.h"
#include "Plane.h"
#include "WorkSpace.h"

#include <loguru.hpp>

#include <functional>

using namespace rxdock;

// Static data members
std::string PolarIdxSF::_CT("PolarIdxSF");
std::string PolarIdxSF::_INCR("INCR");
std::string PolarIdxSF::_ATTR("ATTR");
std::string PolarIdxSF::_THRESHOLD_POS("THRESHOLD_POS");
std::string PolarIdxSF::_THRESHOLD_NEG("THRESHOLD_NEG");

std::string &PolarIdxSF::GetIncr() { return _INCR; }

// NB - Virtual base class constructor (BaseSF) gets called first,
// implicit constructor for BaseInterSF is called second
PolarIdxSF::PolarIdxSF(const std::string &strName)
    : BaseSF(_CT, strName), m_bAttr(true), m_bFlexRec(false), m_bSolvent(false),
      m_nPos(0), m_nNeg(0), m_posThreshold(0.25), m_negThreshold(0.25) {
  LOG_F(2, "PolarIdxSF parameterised constructor");
  // Add parameters
  AddParameter(_INCR, 2.4);
  AddParameter(_ATTR, m_bAttr);
  AddParameter(_THRESHOLD_POS, m_posThreshold);
  AddParameter(_THRESHOLD_NEG, m_negThreshold);
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

PolarIdxSF::~PolarIdxSF() {
  LOG_F(2, "PolarIdxSF destructor");
  ClearReceptor();
  ClearLigand();
  ClearSolvent();
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

// Override BaseSF::ScoreMap to provide additional raw descriptors
void PolarIdxSF::ScoreMap(StringVariantMap &scoreMap) const {
  if (isEnabled()) {
    // XB uncommented next line
    //  EnableAnnotations(true);
    // XB uncommented next line
    //  ClearAnnotationList();
    // We can only annotate the ligand-receptor interactions
    // as the rDock Viewer annotation format is hardwired to expect
    // ligand-receptor atom indices
    double rs = InterScore();
    // XB uncommented next line
    //    EnableAnnotations(false);
    rs += LigandSolventScore();

    // First deal with the inter score which is stored in its natural location
    // in the map
    std::string name = GetFullName();
    scoreMap[name] = rs;
    AddToParentMapEntry(scoreMap, rs);

    // Now deal with the system raw scores which need to be stored in
    // SCORE.INTER.POLAR
    double system_rs =
        ReceptorScore() + SolventScore() + ReceptorSolventScore();
    if (system_rs != 0.0) {
      std::string systemName = BaseSF::_SYSTEM_SF + "." + GetName();
      scoreMap[systemName] = system_rs;
      // increment the parent SCORE.SYSTEM total with the weighted score
      double parentScore = scoreMap[BaseSF::_SYSTEM_SF];
      parentScore += system_rs * GetWeight();
      scoreMap[BaseSF::_SYSTEM_SF] = parentScore;
    }

    // XB uncommented next 6 lines
    //    scoreMap[name+".nhbd"] = m_nPos;
    //    scoreMap[name+".nhba"] = m_nNeg;
    //    std::vector<std::string> annList;
    //    RenderAnnotationList(GetName(),annList);
    //    scoreMap[AnnotationHandler::_ANNOTATION_FIELD] += annList;
    //    ClearAnnotationList();
  }
}

void PolarIdxSF::SetupReceptor() {
  ClearReceptor();
  if (GetReceptor().Null())
    return;
  m_bFlexRec = GetReceptor()->isFlexible();
  double idxIncr = GetParameter(_INCR);
  double maxError = GetMaxError();
  idxIncr += maxError;
  double flexDist =
      2.0; // Additional indexing increment for flexibile OH and NH3
  DockingSitePtr spDS = GetWorkSpace()->GetDockingSite();

  int nCoords = GetReceptor()->GetNumSavedCoords() - 1;
  if (nCoords > 0) {
    AtomList atomList = GetReceptor()->GetAtomList();
    m_spPosGrid = CreateInteractionGrid();
    m_recepPosList = CreateDonorInteractionCenters(atomList);
    m_spNegGrid = CreateInteractionGrid();
    m_recepNegList = CreateAcceptorInteractionCenters(atomList);
    for (int i = 1; i <= nCoords; i++) {
      LOG_F(1, "PolarIdxSF::SetupReceptor: Indexing receptor coords # {}", i);
      GetReceptor()->RevertCoords(i);
      for (InteractionCenterListConstIter iter = m_recepPosList.begin();
           iter != m_recepPosList.end(); iter++) {
        double rvdw = (*iter)->GetAtom1Ptr()->GetVdwRadius();
        m_spPosGrid->SetInteractionLists(*iter, rvdw + idxIncr);
      }
      m_spPosGrid->UniqueInteractionLists();
      for (InteractionCenterListConstIter iter = m_recepNegList.begin();
           iter != m_recepNegList.end(); iter++) {
        double rvdw = (*iter)->GetAtom1Ptr()->GetVdwRadius();
        m_spNegGrid->SetInteractionLists(*iter, rvdw + idxIncr);
      }
      m_spNegGrid->UniqueInteractionLists();
    }
  } else {
    AtomList atomList = spDS->GetAtomList(GetReceptor()->GetAtomList(), 0.0,
                                          GetCorrectedRange());
    m_spPosGrid = CreateInteractionGrid();
    m_spNegGrid = CreateInteractionGrid();
    m_recepPosList = CreateDonorInteractionCenters(atomList);
    m_recepNegList = CreateAcceptorInteractionCenters(atomList);

    // For flexible receptors, separate the interaction centers into rigid and
    // flexible
    if (m_bFlexRec) {
      GetReceptor()->SetAtomSelectionFlags(false);
      GetReceptor()
          ->SelectFlexAtoms(); // This leaves all moveable atoms selected
      isInteractionCenterSelected isSel;
      InteractionCenterListIter fIter = std::stable_partition(
          m_recepPosList.begin(), m_recepPosList.end(), std::not1(isSel));
      std::copy(fIter, m_recepPosList.end(),
                std::back_inserter(m_flexRecPosList));
      m_recepPosList.erase(fIter, m_recepPosList.end());

      fIter = std::stable_partition(m_recepNegList.begin(),
                                    m_recepNegList.end(), std::not1(isSel));
      std::copy(fIter, m_recepNegList.end(),
                std::back_inserter(m_flexRecNegList));
      m_recepNegList.erase(fIter, m_recepNegList.end());

      // Build map of intra-protein interactions (similar to intra-ligand)
      // Include flexible-flexible and flexible-rigid
      int nAtoms = GetReceptor()->GetNumAtoms();
      m_flexRecIntns = InteractionListMap(nAtoms, InteractionCenterList());
      m_flexRecPrtIntns = InteractionListMap(nAtoms, InteractionCenterList());
      if (m_bAttr) {
        BuildIntraMap(m_flexRecPosList, m_flexRecNegList,
                      m_flexRecIntns); // flexible donor - flexible acceptor
        BuildIntraMap(m_flexRecPosList, m_recepNegList,
                      m_flexRecIntns); // flexible donor - rigid acceptor
        BuildIntraMap(m_flexRecNegList, m_recepPosList,
                      m_flexRecIntns); // flexible acceptor - rigid donor
      } else {
        BuildIntraMap(m_flexRecPosList,
                      m_flexRecIntns); // flexible donor - flexible donor
        BuildIntraMap(m_flexRecPosList, m_recepPosList,
                      m_flexRecIntns); // flexible donor - rigid donor
        BuildIntraMap(m_flexRecNegList,
                      m_flexRecIntns); // flexible acceptor - flexible acceptor
        BuildIntraMap(m_flexRecNegList, m_recepNegList,
                      m_flexRecIntns); // flexible acceptor - rigid acceptor
      }
      // We can get away with partitioning the variable interactions just at the
      // beginning. For grosser receptor flexibility we would have to partition
      // periodically during docking
      double partitionDist = GetRange() + flexDist;
      Partition(m_flexRecPosList, m_flexRecNegList, m_flexRecIntns,
                m_flexRecPrtIntns, partitionDist);
      // Index the flexible interaction centers over a larger radius
      // NOTE: WE ASSUME ONLY -OH and -NH3 rotation here (protons can't move
      // more than 2.0A at most) Grosser rotations will require a different
      // approach
      for (InteractionCenterListConstIter iter = m_flexRecPosList.begin();
           iter != m_flexRecPosList.end(); iter++) {
        double rvdw = (*iter)->GetAtom1Ptr()->GetVdwRadius();
        m_spPosGrid->SetInteractionLists(*iter, rvdw + idxIncr + flexDist);
      }
      for (InteractionCenterListConstIter iter = m_flexRecNegList.begin();
           iter != m_flexRecNegList.end(); iter++) {
        double rvdw = (*iter)->GetAtom1Ptr()->GetVdwRadius();
        m_spNegGrid->SetInteractionLists(*iter, rvdw + idxIncr + flexDist);
      }
      LOG_F(1, "{} {}: Intra-receptor score = {}", GetWorkSpace()->GetName(),
            GetFullName(), ReceptorScore());
    }

    // Index the rigid interaction centers as usual
    for (InteractionCenterListConstIter iter = m_recepPosList.begin();
         iter != m_recepPosList.end(); iter++) {
      double rvdw = (*iter)->GetAtom1Ptr()->GetVdwRadius();
      m_spPosGrid->SetInteractionLists(*iter, rvdw + idxIncr);
    }
    for (InteractionCenterListConstIter iter = m_recepNegList.begin();
         iter != m_recepNegList.end(); iter++) {
      double rvdw = (*iter)->GetAtom1Ptr()->GetVdwRadius();
      m_spNegGrid->SetInteractionLists(*iter, rvdw + idxIncr);
    }
  }
}

void PolarIdxSF::SetupLigand() {
  ClearLigand();
  if (GetLigand().Null())
    return;

  AtomList atomList(GetLigand()->GetAtomList());
  m_ligPosList = CreateDonorInteractionCenters(atomList);
  m_ligNegList = CreateAcceptorInteractionCenters(atomList);
}

void PolarIdxSF::SetupSolvent() {
  ClearSolvent();
  ModelList solventList = GetSolvent();
  if (solventList.empty())
    return;

  m_bSolvent = true;

  // Concatenate all solvent atoms into a single list
  AtomList solventAtomList;
  for (ModelListConstIter iter = solventList.begin(); iter != solventList.end();
       ++iter) {
    AtomList atomList = (*iter)->GetAtomList();
    std::copy(atomList.begin(), atomList.end(),
              std::back_inserter(solventAtomList));
  }

  // Create the interaction centers
  m_solventPosList = CreateDonorInteractionCenters(solventAtomList);
  m_solventNegList = CreateAcceptorInteractionCenters(solventAtomList);
  m_solventIntns =
      InteractionListMap(solventAtomList.size(), InteractionCenterList());

  // Build the indexed interaction map
  if (m_bAttr) {
    BuildIntraMap(m_solventPosList, m_solventNegList, m_solventIntns);
  } else {
    BuildIntraMap(m_solventPosList, m_solventIntns);
    BuildIntraMap(m_solventNegList, m_solventIntns);
  }
}

void PolarIdxSF::SetupScore() {
  // No further setup required
}

// Overall score is the sum of:
// ligand-receptor
// ligand-solvent
// intra-receptor
// intra-solvent
// receptor-solvent
double PolarIdxSF::RawScore() const {
  return InterScore() + LigandSolventScore() + ReceptorScore() +
         SolventScore() + ReceptorSolventScore();
}

// Clear the receptor and ligand grids and lists respectively
// As we are not using smart pointers, there is some memory management to do
void PolarIdxSF::ClearReceptor() {
  m_spPosGrid = InteractionGridPtr();
  m_spNegGrid = InteractionGridPtr();
  m_flexRecIntns.clear();
  m_flexRecPrtIntns.clear();
  m_bFlexRec = false;
  DeleteList(m_recepPosList);
  DeleteList(m_flexRecPosList);
  DeleteList(m_recepNegList);
  DeleteList(m_flexRecNegList);
}

void PolarIdxSF::ClearLigand() {
  DeleteList(m_ligPosList);
  DeleteList(m_ligNegList);
}

void PolarIdxSF::ClearSolvent() {
  m_solventIntns.clear();
  m_bSolvent = false;
  DeleteList(m_solventPosList);
  DeleteList(m_solventNegList);
}

// Deletes the interaction center objects contained in the list, then clears the
// list
void PolarIdxSF::DeleteList(InteractionCenterList &icList) {
  for (InteractionCenterListIter iter = icList.begin(); iter != icList.end();
       ++iter) {
    delete *iter;
  }
  icList.clear();
}

// DM 25 Oct 2000 - track changes to parameter values in local data members
// ParameterUpdated is invoked by ParamHandler::SetParameter
void PolarIdxSF::ParameterUpdated(const std::string &strName) {
  // DM 25 Oct 2000 - heavily used params
  if (strName == _ATTR) {
    m_bAttr = GetParameter(_ATTR);
  } else if (strName == _THRESHOLD_POS) {
    m_posThreshold = GetParameter(_THRESHOLD_POS);
  } else if (strName == _THRESHOLD_NEG) {
    m_negThreshold = GetParameter(_THRESHOLD_NEG);
  } else {
    PolarSF::OwnParameterUpdated(strName);
    BaseIdxSF::OwnParameterUpdated(strName);
    BaseSF::ParameterUpdated(strName);
  }
}

// Intra-receptor
double PolarIdxSF::ReceptorScore() const {
  return (m_bFlexRec) ? IntraScore(m_flexRecPosList, m_flexRecNegList,
                                   m_flexRecPrtIntns, m_bAttr)
                      : 0.0;
}

// Intra-solvent
double PolarIdxSF::SolventScore() const {
  return (m_bSolvent) ? IntraScore(m_solventPosList, m_solventNegList,
                                   m_solventIntns, m_bAttr)
                      : 0.0;
}

// Ligand-receptor
double PolarIdxSF::InterScore() const {
  return InterScore(m_ligPosList, m_ligNegList, true);
}

// Receptor-solvent
double PolarIdxSF::ReceptorSolventScore() const {
  return (m_bSolvent) ? InterScore(m_solventPosList, m_solventNegList, false)
                      : 0.0;
}

// Ligand-solvent score: very inefficient
double PolarIdxSF::LigandSolventScore() const {
  double score = 0.0;
  if (!m_bSolvent)
    return score;
  PolarSF::f1prms Rprms = GetRprms();   // Distance params
  PolarSF::f1prms A1prms = GetA1prms(); // Donor angle params
  PolarSF::f1prms A2prms = GetA2prms(); // Acceptor angle params
  if (m_bAttr) {
    // Pos-neg
    for (InteractionCenterListConstIter iter = m_solventPosList.begin();
         iter != m_solventPosList.end(); iter++) {
      score += (*iter)->GetAtom1Ptr()->GetUser1Value() *
               PolarScore(*iter, m_ligNegList, Rprms, A1prms, A2prms);
    }
    // Neg-pos
    for (InteractionCenterListConstIter iter = m_solventNegList.begin();
         iter != m_solventNegList.end(); iter++) {
      score += (*iter)->GetAtom1Ptr()->GetUser1Value() *
               PolarScore(*iter, m_ligPosList, Rprms, A2prms, A1prms);
    }
  } else {
    // pos-pos
    for (InteractionCenterListConstIter iter = m_solventPosList.begin();
         iter != m_solventPosList.end(); iter++) {
      score += (*iter)->GetAtom1Ptr()->GetUser1Value() *
               PolarScore(*iter, m_ligPosList, Rprms, A1prms, A1prms);
    }
    // neg-neg
    for (InteractionCenterListConstIter iter = m_solventNegList.begin();
         iter != m_solventNegList.end(); iter++) {
      score += (*iter)->GetAtom1Ptr()->GetUser1Value() *
               PolarScore(*iter, m_ligNegList, Rprms, A2prms, A2prms);
    }
  }
  return score;
}

// Reusable method for receptor-ligand and receptor-solvent scores
// bCount controls whether to count the positive and negative interaction scores
double PolarIdxSF::InterScore(const InteractionCenterList &posList,
                              const InteractionCenterList &negList,
                              bool bCount) const {
  double score = 0.0; // Total score
  if (bCount) {
    m_nPos = 0;
    m_nNeg = 0;
  }

  // Check grid is defined
  if (m_spPosGrid.Null() || m_spNegGrid.Null())
    return score;

  PolarSF::f1prms Rprms = GetRprms();   // Distance params
  PolarSF::f1prms A1prms = GetA1prms(); // Donor angle params
  PolarSF::f1prms A2prms = GetA2prms(); // Acceptor angle params

  double s(0.0); // Partial scores
  // Ligand HBA
  for (InteractionCenterListConstIter lIter = negList.begin();
       lIter != negList.end(); lIter++) {
    Atom *pLig1 = (*lIter)->GetAtom1Ptr();
    const Coord &cLig1 = pLig1->GetCoords();
    // If this is an attractive potential we calculate the score with all
    // adjacent +ve centres (HBD/M+/guan)
    if (m_bAttr) {
      const InteractionCenterList &rList =
          m_spPosGrid->GetInteractionList(cLig1);
      s = PolarScore(*lIter, rList, Rprms, A2prms, A1prms);
    } else {
      // If this is an repulsive potential we calculate the score with all
      // adjacent HBA
      const InteractionCenterList &rList =
          m_spNegGrid->GetInteractionList(cLig1);
      s = PolarScore(*lIter, rList, Rprms, A2prms, A2prms);
    }
    s *= pLig1->GetUser1Value();
    if (bCount && (std::fabs(s) > m_negThreshold)) {
      m_nNeg++;
    }
    score += s;
  }

  // Ligand HBD
  for (InteractionCenterListConstIter lIter = posList.begin();
       lIter != posList.end(); lIter++) {
    Atom *pLig1 = (*lIter)->GetAtom1Ptr();
    const Coord &cLig1 = pLig1->GetCoords();
    // If this is an attractive potential we calculate the score with all
    // adjacent HBA
    if (m_bAttr) {
      const InteractionCenterList &rList =
          m_spNegGrid->GetInteractionList(cLig1);
      s = PolarScore(*lIter, rList, Rprms, A1prms, A2prms);
    } else {
      // If this is an repulsive potential we calculate the score with all
      // adjacent +ve centres (HBD/M+/guan)
      const InteractionCenterList &rList =
          m_spPosGrid->GetInteractionList(cLig1);
      s = PolarScore(*lIter, rList, Rprms, A1prms, A1prms);
    }
    s *= pLig1->GetUser1Value();
    if (bCount && (std::fabs(s) > m_posThreshold)) {
      m_nPos++;
    }
    score += s;
  }
  return score;
}
