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

#include "ChromFactory.h"
#include "Chrom.h"
#include "ChromDihedralElement.h"
#include "ChromOccupancyElement.h"
#include "ChromPositionElement.h"
#include "LigandFlexData.h"
#include "Model.h"
#include "ReceptorFlexData.h"
#include "SolventFlexData.h"

#include <loguru.hpp>

using namespace rxdock;

ChromFactory::ChromFactory() { m_pChrom = new Chrom(); }

ChromElement *ChromFactory::GetChrom() const { return m_pChrom; }

void ChromFactory::VisitReceptorFlexData(ReceptorFlexData *pFlexData) {
  Model *pModel = pFlexData->GetModel();
  DockingSite *pDockSite = pFlexData->GetDockingSite();
  if (pModel && pDockSite) {
    double flexDistance =
        pFlexData->GetParameter(ReceptorFlexData::_FLEX_DISTANCE);
    double dihedralStepSize =
        pFlexData->GetParameter(ReceptorFlexData::_DIHEDRAL_STEP);
    // Trap the combination of flexible OH/NH3 AND multiple receptor
    // conformations We do not support both of these simultaneously
    if (pModel->GetNumSavedCoords() > 1) {
      std::string message(
          "The combination of flexible OH/NH3 groups AND multiple receptor "
          "conformations is not supported currently");
      throw InvalidRequest(_WHERE_, message);
    }

    // Find all the terminal OH and NH3+ bonds within range of the docking
    // volume
    BondList rotBondList = GetRotatableBondList(pModel->GetBondList());
    DockingSite::isAtomInRange bIsInRange(pDockSite->GetGrid(), 0.0,
                                          flexDistance);
    isBondToOH isOH;
    isBondToNH3 isNH3;
    BondList modelMutatorBondList;
    AtomList noTetheredAtoms;
    for (BondListConstIter iter = rotBondList.begin();
         iter != rotBondList.end(); ++iter) {
      if (isOH(*iter) || isNH3(*iter)) {
        Atom *pAtom1 = (*iter)->GetAtom1Ptr();
        Atom *pAtom2 = (*iter)->GetAtom2Ptr();
        if (bIsInRange(pAtom1) && bIsInRange(pAtom2)) {
          m_pChrom->Add(new ChromDihedralElement(*iter, noTetheredAtoms,
                                                 dihedralStepSize));
          modelMutatorBondList.push_back(*iter);
        }
      }
    }
    // Create the legacy ModelMutator object
    // needed for storing the flexible interaction maps
    if (!modelMutatorBondList.empty()) {
      m_spMutator = ModelMutatorPtr(
          new ModelMutator(pModel, modelMutatorBondList, noTetheredAtoms));
    } else {
      m_spMutator.SetNull();
    }
  }
}

void ChromFactory::VisitLigandFlexData(LigandFlexData *pFlexData) {
  Model *pModel = pFlexData->GetModel();
  DockingSite *pDockSite = pFlexData->GetDockingSite();
  if (pModel && pDockSite) {
    double transStepSize = pFlexData->GetParameter(LigandFlexData::_TRANS_STEP);
    double rotStepSize = pFlexData->GetParameter(LigandFlexData::_ROT_STEP);
    double dihedralStepSize =
        pFlexData->GetParameter(LigandFlexData::_DIHEDRAL_STEP);
    std::string transModeStr =
        pFlexData->GetParameter(LigandFlexData::_TRANS_MODE);
    std::string rotModeStr = pFlexData->GetParameter(LigandFlexData::_ROT_MODE);
    std::string dihedralModeStr =
        pFlexData->GetParameter(LigandFlexData::_DIHEDRAL_MODE);
    double maxTrans = pFlexData->GetParameter(LigandFlexData::_MAX_TRANS);
    double maxRot = pFlexData->GetParameter(LigandFlexData::_MAX_ROT);
    double maxDihedral = pFlexData->GetParameter(LigandFlexData::_MAX_DIHEDRAL);

    // Convert from sampling mode strings to enum values
    ChromElement::eMode transMode = ChromElement::StrToMode(transModeStr);
    ChromElement::eMode rotMode = ChromElement::StrToMode(rotModeStr);
    ChromElement::eMode dihedralMode = ChromElement::StrToMode(dihedralModeStr);

    AtomList tetheredAtoms = pModel->GetTetheredAtomList();
    BondList rotBondList =
        GetBondListWithPredicate(pModel->GetBondList(), isBondRotatable());

    // If we are in tethered mode, ensure the step size is not larger than the
    // tethered range
    if ((transMode == ChromElement::TETHERED) && (transStepSize > maxTrans)) {
      transStepSize = maxTrans;
      pFlexData->SetParameter(LigandFlexData::_TRANS_STEP, transStepSize);
    }
    if ((rotMode == ChromElement::TETHERED) && (rotStepSize > maxRot)) {
      rotStepSize = maxRot;
      pFlexData->SetParameter(LigandFlexData::_ROT_STEP, rotStepSize);
    }
    if ((dihedralMode == ChromElement::TETHERED) &&
        (dihedralStepSize > maxDihedral)) {
      dihedralStepSize = maxDihedral;
      pFlexData->SetParameter(LigandFlexData::_DIHEDRAL_STEP, dihedralStepSize);
    }

    // Dihedrals
    if (dihedralMode != ChromElement::FIXED) {
      for (BondListConstIter iter = rotBondList.begin();
           iter != rotBondList.end(); ++iter) {
        m_pChrom->Add(new ChromDihedralElement(
            *iter, tetheredAtoms, dihedralStepSize, dihedralMode, maxDihedral));
      }
    }

    // Position
    // TODO: Check this code is valid. Does not add the chromosome element
    // if translation and rotation are fixed.
    if ((transMode != ChromElement::FIXED) ||
        (rotMode != ChromElement::FIXED)) {
      // Don't forget that whole body rotation code is in radians (not degrees)
      m_pChrom->Add(new ChromPositionElement(
          pModel, pDockSite, transStepSize, rotStepSize * M_PI / 180.0,
          transMode, rotMode, maxTrans, maxRot * M_PI / 180.0));
    }
    // Create the legacy ModelMutator object
    // needed for storing the flexible interaction maps
    if (!rotBondList.empty()) {
      m_spMutator =
          ModelMutatorPtr(new ModelMutator(pModel, rotBondList, tetheredAtoms));
    } else {
      m_spMutator.SetNull();
    }
  }
}

void ChromFactory::VisitSolventFlexData(SolventFlexData *pFlexData) {
  VisitLigandFlexData(pFlexData);
  Model *pModel = pFlexData->GetModel();
  DockingSite *pDockSite = pFlexData->GetDockingSite();
  if (pModel && pDockSite) {
    double occupancyProb = pFlexData->GetParameter(SolventFlexData::_OCCUPANCY);
    double stepSize = pFlexData->GetParameter(SolventFlexData::_OCCUPANCY_STEP);
    // If occupancy prob is non-zero, and less than one, then add
    // an occupancy element to the chromosome
    if ((occupancyProb > 0.0) && (occupancyProb < 1.0)) {
      double threshold = 1.0 - occupancyProb;
      m_pChrom->Add(new ChromOccupancyElement(pModel, stepSize, threshold));
      LOG_F(INFO, "Solvent model {} has variable occupancy", pModel->GetName());
    }
    // occupancy prob is either zero or one.
    // just need to force the model into either the enabled or disabled state
    // No need for chromosome element
    else if (occupancyProb <= 0.0) {
      pModel->SetOccupancy(0.0);
      LOG_F(WARNING, "Solvent model {} is permanently disabled",
            pModel->GetName());
    } else if (occupancyProb >= 1.0) {
      pModel->SetOccupancy(1.0);
      LOG_F(INFO, "Solvent model {} is permanently enabled", pModel->GetName());
    }
  }
}
