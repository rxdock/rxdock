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

#include "FlexAtomFactory.h"
#include "ChromElement.h"
#include "LigandFlexData.h"
#include "ReceptorFlexData.h"
#include "SolventFlexData.h"

using namespace rxdock;

FlexAtomFactory::FlexAtomFactory(Model *pModel) { Visit(pModel); }

FlexAtomFactory::FlexAtomFactory(ModelList modelList) {
  for (ModelListConstIter iter = modelList.begin(); iter != modelList.end();
       iter++) {
    Visit(*iter);
  }
}

void FlexAtomFactory::Visit(Model *pModel) {
  if (pModel != nullptr) {
    FlexData *pFlexData = pModel->GetFlexData();
    if (pFlexData != nullptr) {
      pFlexData->Accept(*this);
    } else {
      // If flex data is null, then model is rigid by default
      // Therefore all atoms are fixed
      AtomList atomList = pModel->GetAtomList();
      pModel->SetAtomUser2Values(0.0);
      std::copy(atomList.begin(), atomList.end(),
                std::back_inserter(m_fixedAtomList));
    }
  }
}

void FlexAtomFactory::Clear() {
  m_fixedAtomList.clear();
  m_tetheredAtomList.clear();
  m_freeAtomList.clear();
}

void FlexAtomFactory::VisitReceptorFlexData(ReceptorFlexData *pFlexData) {
  Model *pModel = pFlexData->GetModel();
  AtomList atomList = pModel->GetAtomList();
  // Divide the receptor atoms into fixed and tethered (pendant OH/NH3+)
  pModel->SetAtomSelectionFlags(false);
  pModel->SetAtomUser2Values(0.0);
  pModel->SelectFlexAtoms();
  for (AtomListIter iter = atomList.begin(); iter != atomList.end(); ++iter) {
    if ((*iter)->GetSelectionFlag()) {
      m_tetheredAtomList.push_back(*iter);
      // This code assumes that the only flexible atoms in the receptor are
      // terminal OH/NH3+ hydrogens. Max rotational displacement from current
      // coords is equal to twice the bond length.
      AtomList bondedAtomList = GetBondedAtomList(*iter);
      if (bondedAtomList.size() == 1) {
        Atom *parent = bondedAtomList.front();
        double bondLength = Length(parent->GetCoords(), (*iter)->GetCoords());
        (*iter)->SetUser2Value(2.0 * bondLength);
      }
    } else {
      m_fixedAtomList.push_back(*iter);
    }
  }
}

void FlexAtomFactory::VisitLigandFlexData(LigandFlexData *pFlexData) {
  Model *pModel = pFlexData->GetModel();
  AtomList atomList = pModel->GetAtomList();
  pModel->SetAtomUser2Values(0.0);
  std::string transModeStr =
      pFlexData->GetParameter(LigandFlexData::_TRANS_MODE);
  std::string rotModeStr = pFlexData->GetParameter(LigandFlexData::_ROT_MODE);
  ChromElement::eMode eTransMode = ChromElement::StrToMode(transModeStr);
  ChromElement::eMode eRotMode = ChromElement::StrToMode(rotModeStr);
  // Internally flexible ligand, or free translation
  if (pModel->isFlexible() || (eTransMode == ChromElement::FREE)) {
    std::copy(atomList.begin(), atomList.end(),
              std::back_inserter(m_freeAtomList));
  } else if (eTransMode == ChromElement::FIXED) {
    // Completely fixed
    if (eRotMode == ChromElement::FIXED) {
      std::copy(atomList.begin(), atomList.end(),
                std::back_inserter(m_fixedAtomList));
    }
    // Fixed translation; tethered or free rotation
    // Therefore max translation of any atom depends on its distance from the
    // COM For simplicity, we don't distinguish between tethered and free
    // rotation so assume full rotation through +/- 180 deg is possible.
    else {
      PrincipalAxes prAxes = GetPrincipalAxesOfAtoms(atomList);
      for (AtomListIter iter = atomList.begin(); iter != atomList.end();
           ++iter) {
        m_tetheredAtomList.push_back(*iter);
        double radius = Length((*iter)->GetCoords(), prAxes.com);
        (*iter)->SetUser2Value(2.0 * radius);
      }
    }
  } else if (eTransMode == ChromElement::TETHERED) {
    double maxTrans = pFlexData->GetParameter(LigandFlexData::_MAX_TRANS);
    // Tethered translation; fixed rotation
    // Therefore max translation of any atom is a constant (+/- maxTrans)
    if (eRotMode == ChromElement::FIXED) {
      std::copy(atomList.begin(), atomList.end(),
                std::back_inserter(m_tetheredAtomList));
      pModel->SetAtomUser2Values(maxTrans);
    }
    // Tethered translation, tethered or free rotation
    // Therefore max translation of any atom depends on its distance from the
    // COM and on maxTrans For simplicity, we don't distinguish between tethered
    // and free rotation so assume full rotation through +/- 180 deg is
    // possible.
    else {
      PrincipalAxes prAxes = GetPrincipalAxesOfAtoms(atomList);
      for (AtomListIter iter = atomList.begin(); iter != atomList.end();
           ++iter) {
        m_tetheredAtomList.push_back(*iter);
        double radius = Length((*iter)->GetCoords(), prAxes.com);
        (*iter)->SetUser2Value((2.0 * radius) + maxTrans);
      }
    }
  }
}

void FlexAtomFactory::VisitSolventFlexData(SolventFlexData *pFlexData) {
  VisitLigandFlexData(pFlexData);
}
