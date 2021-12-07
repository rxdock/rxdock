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

#include "rxdock/CavityGridSF.h"
#include "rxdock/ChromPositionRefData.h"
#include "rxdock/LigandFlexData.h"
#include "rxdock/ReceptorFlexData.h"
#include "rxdock/SolventFlexData.h"
#include "rxdock/WorkSpace.h"

#include <functional>
#include <loguru.hpp>

using namespace rxdock;

// Static data members
const std::string CavityGridSF::_CT = "CavityGridSF";
const std::string CavityGridSF::_RMAX = "rmax";
const std::string CavityGridSF::_QUADRATIC = "quadratic";

// NB - Virtual base class constructor (BaseSF) gets called first,
// implicit constructor for BaseInterSF is called second
CavityGridSF::CavityGridSF(const std::string &strName)
    : BaseSF(_CT, strName), m_maxDist(0.0), m_rMax(0.1), m_bQuadratic(false) {
  LOG_F(2, "CavityGridSF parameterised constructor");
  // Add parameters
  AddParameter(_RMAX, m_rMax);
  AddParameter(_QUADRATIC, m_bQuadratic);
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

CavityGridSF::~CavityGridSF() {
  LOG_F(2, "CavityGridSF destructor");
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

void CavityGridSF::SetupReceptor() {
  // DM 09 Apr 2002 - we can reuse the precalculated distance grid directly from
  // the docking site
  m_spGrid = RealGridPtr();
  DockingSitePtr spDS = GetWorkSpace()->GetDockingSite();
  if (spDS.Null())
    return;
  m_spGrid = spDS->GetGrid();
  if (m_spGrid.Null())
    return;
  m_maxDist = m_spGrid->MaxValue();
}

void CavityGridSF::SetupLigand() {}

void CavityGridSF::SetupSolvent() {}

void CavityGridSF::SetupScore() {
  m_atomList.clear();
  WorkSpace *pWorkSpace = GetWorkSpace();
  if (pWorkSpace) {
    HeavyAtomFactory atomFactory(pWorkSpace->GetModels());
    m_atomList = atomFactory.GetAtomList();
  }
}

double CavityGridSF::RawScore() const {
  double score(0.0);

  // Check grid is defined
  if (m_spGrid.Null())
    return score;

  // Quadratic penalty function
  if (m_bQuadratic) {
    for (AtomRListConstIter iter = m_atomList.begin(); iter != m_atomList.end();
         iter++) {
      const Coord &c = (*iter)->GetCoords();
      // Check if atom is off grid, if so default to max grid value
      double r =
          m_spGrid->isValid(c) ? m_spGrid->GetSmoothedValue(c) : m_maxDist;
      double dr = r - m_rMax;
      if (dr > 0.0) {
        score += dr * dr;
      }
    }
  }
  // Linear penalty function
  else {
    for (AtomRListConstIter iter = m_atomList.begin(); iter != m_atomList.end();
         iter++) {
      const Coord &c = (*iter)->GetCoords();
      // Check if atom is off grid, if so default to max grid value
      double r =
          m_spGrid->isValid(c) ? m_spGrid->GetSmoothedValue(c) : m_maxDist;
      double dr = r - m_rMax;
      if (dr > 0.0) {
        score += dr;
      }
    }
  }
  return score;
}

// DM 25 Oct 2000 - track changes to parameter values in local data members
// ParameterUpdated is invoked by ParamHandler::SetParameter
void CavityGridSF::ParameterUpdated(const std::string &strName) {
  if (strName == _RMAX) {
    m_rMax = GetParameter(_RMAX);
  } else if (strName == _QUADRATIC) {
    m_bQuadratic = GetParameter(_QUADRATIC);
  } else {
    BaseSF::ParameterUpdated(strName);
  }
}

CavityGridSF::HeavyAtomFactory::HeavyAtomFactory(ModelList modelList) {
  for (ModelListConstIter iter = modelList.begin(); iter != modelList.end();
       iter++) {
    if ((*iter).Ptr() != nullptr) {
      FlexData *pFlexData = (*iter)->GetFlexData();
      if (pFlexData != nullptr) {
        pFlexData->Accept(*this);
      }
    }
  }
}

void CavityGridSF::HeavyAtomFactory::VisitReceptorFlexData(
    ReceptorFlexData *pFlexData) {
  // Receptor atoms: do nothing (the only movable atoms are hydrogens at
  // present)
}

void CavityGridSF::HeavyAtomFactory::VisitLigandFlexData(
    LigandFlexData *pFlexData) {
  // Extract all non-hydrogen atoms from ligand
  AtomList atomList = pFlexData->GetModel()->GetAtomList();
  std::copy_if(atomList.begin(), atomList.end(), std::back_inserter(m_atomList),
               std::not1(isAtomicNo_eq(1)));
}

void CavityGridSF::HeavyAtomFactory::VisitSolventFlexData(
    SolventFlexData *pFlexData) {
  // Extract all non-hydrogen atoms from solvent but only if the solvent can
  // translate. As there are no dihedrals to worry about, the water oxygens are
  // fixed in position for orientation-only flexibility.
  std::string transModeStr =
      pFlexData->GetParameter(LigandFlexData::_TRANS_MODE);
  double maxTrans = pFlexData->GetParameter(LigandFlexData::_MAX_TRANS);
  ChromElement::eMode eTransMode = ChromElement::StrToMode(transModeStr);
  if ((eTransMode == ChromElement::FREE) ||
      ((eTransMode == ChromElement::TETHERED) && (maxTrans > 0.0))) {
    AtomList atomList = pFlexData->GetModel()->GetAtomList();
    std::copy_if(atomList.begin(), atomList.end(),
                 std::back_inserter(m_atomList), std::not1(isAtomicNo_eq(1)));
  }
}
