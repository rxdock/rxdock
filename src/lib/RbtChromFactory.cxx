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

#include "RbtChromFactory.h"
#include "RbtChrom.h"
#include "RbtChromDihedralElement.h"
#include "RbtChromOccupancyElement.h"
#include "RbtChromPositionElement.h"
#include "RbtLigandFlexData.h"
#include "RbtModel.h"
#include "RbtReceptorFlexData.h"
#include "RbtSolventFlexData.h"

using namespace rxdock;

RbtChromFactory::RbtChromFactory() { m_pChrom = new RbtChrom(); }

RbtChromElement *RbtChromFactory::GetChrom() const { return m_pChrom; }

void RbtChromFactory::VisitReceptorFlexData(RbtReceptorFlexData *pFlexData) {
  RbtModel *pModel = pFlexData->GetModel();
  RbtDockingSite *pDockSite = pFlexData->GetDockingSite();
  if (pModel && pDockSite) {
    double flexDistance =
        pFlexData->GetParameter(RbtReceptorFlexData::_FLEX_DISTANCE);
    double dihedralStepSize =
        pFlexData->GetParameter(RbtReceptorFlexData::_DIHEDRAL_STEP);
    // Trap the combination of flexible OH/NH3 AND multiple receptor
    // conformations We do not support both of these simultaneously
    if (pModel->GetNumSavedCoords() > 1) {
      std::string message(
          "The combination of flexible OH/NH3 groups AND multiple receptor "
          "conformations is not supported currently");
      throw RbtInvalidRequest(_WHERE_, message);
    }

    // Find all the terminal OH and NH3+ bonds within range of the docking
    // volume
    RbtBondList rotBondList = GetRotatableBondList(pModel->GetBondList());
    RbtDockingSite::isAtomInRange bIsInRange(pDockSite->GetGrid(), 0.0,
                                             flexDistance);
    isBondToOH isOH;
    isBondToNH3 isNH3;
    RbtBondList modelMutatorBondList;
    RbtAtomList noTetheredAtoms;
    for (RbtBondListConstIter iter = rotBondList.begin();
         iter != rotBondList.end(); ++iter) {
      if (isOH(*iter) || isNH3(*iter)) {
        RbtAtom *pAtom1 = (*iter)->GetAtom1Ptr();
        RbtAtom *pAtom2 = (*iter)->GetAtom2Ptr();
        if (bIsInRange(pAtom1) && bIsInRange(pAtom2)) {
          m_pChrom->Add(new RbtChromDihedralElement(*iter, noTetheredAtoms,
                                                    dihedralStepSize));
          modelMutatorBondList.push_back(*iter);
        }
      }
    }
    // Create the legacy ModelMutator object
    // needed for storing the flexible interaction maps
    if (!modelMutatorBondList.empty()) {
      m_spMutator = RbtModelMutatorPtr(
          new RbtModelMutator(pModel, modelMutatorBondList, noTetheredAtoms));
    } else {
      m_spMutator.SetNull();
    }
  }
}

void RbtChromFactory::VisitLigandFlexData(RbtLigandFlexData *pFlexData) {
  RbtModel *pModel = pFlexData->GetModel();
  RbtDockingSite *pDockSite = pFlexData->GetDockingSite();
  if (pModel && pDockSite) {
    double transStepSize =
        pFlexData->GetParameter(RbtLigandFlexData::_TRANS_STEP);
    double rotStepSize = pFlexData->GetParameter(RbtLigandFlexData::_ROT_STEP);
    double dihedralStepSize =
        pFlexData->GetParameter(RbtLigandFlexData::_DIHEDRAL_STEP);
    std::string transModeStr =
        pFlexData->GetParameter(RbtLigandFlexData::_TRANS_MODE);
    std::string rotModeStr =
        pFlexData->GetParameter(RbtLigandFlexData::_ROT_MODE);
    std::string dihedralModeStr =
        pFlexData->GetParameter(RbtLigandFlexData::_DIHEDRAL_MODE);
    double maxTrans = pFlexData->GetParameter(RbtLigandFlexData::_MAX_TRANS);
    double maxRot = pFlexData->GetParameter(RbtLigandFlexData::_MAX_ROT);
    double maxDihedral =
        pFlexData->GetParameter(RbtLigandFlexData::_MAX_DIHEDRAL);

    // Convert from sampling mode strings to enum values
    RbtChromElement::eMode transMode = RbtChromElement::StrToMode(transModeStr);
    RbtChromElement::eMode rotMode = RbtChromElement::StrToMode(rotModeStr);
    RbtChromElement::eMode dihedralMode =
        RbtChromElement::StrToMode(dihedralModeStr);

    RbtAtomList tetheredAtoms = pModel->GetTetheredAtomList();
    RbtBondList rotBondList =
        GetBondListWithPredicate(pModel->GetBondList(), isBondRotatable());

    // If we are in tethered mode, ensure the step size is not larger than the
    // tethered range
    if ((transMode == RbtChromElement::TETHERED) &&
        (transStepSize > maxTrans)) {
      transStepSize = maxTrans;
      pFlexData->SetParameter(RbtLigandFlexData::_TRANS_STEP, transStepSize);
    }
    if ((rotMode == RbtChromElement::TETHERED) && (rotStepSize > maxRot)) {
      rotStepSize = maxRot;
      pFlexData->SetParameter(RbtLigandFlexData::_ROT_STEP, rotStepSize);
    }
    if ((dihedralMode == RbtChromElement::TETHERED) &&
        (dihedralStepSize > maxDihedral)) {
      dihedralStepSize = maxDihedral;
      pFlexData->SetParameter(RbtLigandFlexData::_DIHEDRAL_STEP,
                              dihedralStepSize);
    }

    // Dihedrals
    if (dihedralMode != RbtChromElement::FIXED) {
      for (RbtBondListConstIter iter = rotBondList.begin();
           iter != rotBondList.end(); ++iter) {
        m_pChrom->Add(new RbtChromDihedralElement(
            *iter, tetheredAtoms, dihedralStepSize, dihedralMode, maxDihedral));
      }
    }

    // Position
    // TODO: Check this code is valid. Does not add the chromosome element
    // if translation and rotation are fixed.
    if ((transMode != RbtChromElement::FIXED) ||
        (rotMode != RbtChromElement::FIXED)) {
      // Don't forget that whole body rotation code is in radians (not degrees)
      m_pChrom->Add(new RbtChromPositionElement(
          pModel, pDockSite, transStepSize, rotStepSize * M_PI / 180.0,
          transMode, rotMode, maxTrans, maxRot * M_PI / 180.0));
    }
    // Create the legacy ModelMutator object
    // needed for storing the flexible interaction maps
    if (!rotBondList.empty()) {
      m_spMutator = RbtModelMutatorPtr(
          new RbtModelMutator(pModel, rotBondList, tetheredAtoms));
    } else {
      m_spMutator.SetNull();
    }
  }
}

void RbtChromFactory::VisitSolventFlexData(RbtSolventFlexData *pFlexData) {
  VisitLigandFlexData(pFlexData);
  RbtModel *pModel = pFlexData->GetModel();
  RbtDockingSite *pDockSite = pFlexData->GetDockingSite();
  if (pModel && pDockSite) {
    double occupancyProb =
        pFlexData->GetParameter(RbtSolventFlexData::_OCCUPANCY);
    double stepSize =
        pFlexData->GetParameter(RbtSolventFlexData::_OCCUPANCY_STEP);
    // If occupancy prob is non-zero, and less than one, then add
    // an occupancy element to the chromosome
    if ((occupancyProb > 0.0) && (occupancyProb < 1.0)) {
      double threshold = 1.0 - occupancyProb;
      m_pChrom->Add(new RbtChromOccupancyElement(pModel, stepSize, threshold));
      // std::cout << "INFO Solvent model " << pModel->GetName() << " has
      // variable occupancy" << std::endl;
    }
    // occupancy prob is either zero or one.
    // just need to force the model into either the enabled or disabled state
    // No need for chromosome element
    else if (occupancyProb <= 0.0) {
      pModel->SetOccupancy(0.0);
      std::cout << "WARNING Solvent model " << pModel->GetName()
                << " is permanently disabled" << std::endl;
    } else if (occupancyProb >= 1.0) {
      pModel->SetOccupancy(1.0);
      // std::cout << "INFO Solvent model " << pModel->GetName() << " is
      // permanently enabled" << std::endl;
    }
  }
}
