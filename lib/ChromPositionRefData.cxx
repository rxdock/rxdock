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

#include "rxdock/ChromPositionRefData.h"

using namespace rxdock;

const std::string ChromPositionRefData::_CT = "ChromPositionRefData";
const PrincipalAxes ChromPositionRefData::CARTESIAN_AXES;

ChromPositionRefData::ChromPositionRefData(
    const Model *pModel, const DockingSite *pDockSite, double transStepSize,
    double rotStepSize, ChromElement::eMode transMode,
    ChromElement::eMode rotMode, double maxTrans, double maxRot)
    : m_transStepSize(transStepSize), m_rotStepSize(rotStepSize),
      m_transMode(transMode), m_rotMode(rotMode), m_length(6), m_xOverLength(2),
      m_maxTrans(maxTrans), m_maxRot(maxRot) {
  AtomList atomList = pModel->GetAtomList();
  // Tethered substructure atom list (may be empty)
  AtomList tetheredAtomList = pModel->GetTetheredAtomList();
  // If we have a tethered substructure, use this as the reference atom list
  // to calculate centre of mass and orientation, else use all atoms
  m_refAtoms = (tetheredAtomList.empty()) ? atomList : tetheredAtomList;
  // All atoms are movable, but use std::copy to strip off the smart pointers
  std::copy(atomList.begin(), atomList.end(),
            std::back_inserter(m_movableAtoms));
  GetModelValue(m_initialCom, m_initialOrientation);
  m_initialQuat = m_initialOrientation.ToQuat();
  pDockSite->GetCoordList(m_startCoords);
  // Check for zero ranges in TETHERED mode and convert to FIXED
  if ((m_transMode == ChromElement::TETHERED) && (m_maxTrans <= 0.0)) {
    m_transMode = ChromElement::FIXED;
    m_maxTrans = 0.0;
  }
  if ((m_rotMode == ChromElement::TETHERED) && (m_maxRot <= 0.0)) {
    m_rotMode = ChromElement::FIXED;
    m_maxRot = 0.0;
  }
  // Remove degrees of freedom for fixed modes
  if (IsTransFixed()) {
    m_length -= 3;
    m_xOverLength--;
  }
  if (IsRotFixed()) {
    m_length -= 3;
    m_xOverLength--;
  }
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

ChromPositionRefData::~ChromPositionRefData() { _RBTOBJECTCOUNTER_DESTR_(_CT); }

void ChromPositionRefData::GetModelValue(Coord &com, Euler &orientation) const {
  // Determine the principal axes and centre of mass of the reference atoms
  PrincipalAxes prAxes = GetPrincipalAxesOfAtoms(m_refAtoms);
  // Determine the quaternion needed to align Cartesian axes with actual
  // molecule principal axes. This represents the absolute orientation of
  // the molecule.
  Quat q = GetQuatFromAlignAxes(CARTESIAN_AXES, prAxes);
  orientation.FromQuat(q);
  com = prAxes.com;
}

void ChromPositionRefData::SetModelValue(const Coord &com,
                                         const Euler &orientation) {
  // Determine the principal axes and centre of mass of the reference atoms
  PrincipalAxes prAxes = GetPrincipalAxesOfAtoms(m_refAtoms);
  // Determine the overall rotation required.
  // 1) Go back to realign with Cartesian axes
  Quat qBack = GetQuatFromAlignAxes(prAxes, CARTESIAN_AXES);
  // 2) Go forward to the desired orientation
  Quat qForward = orientation.ToQuat();
  // 3 Combine the two rotations
  Quat q = qForward * qBack;
  for (AtomRListIter iter = m_movableAtoms.begin();
       iter != m_movableAtoms.end(); ++iter) {
    (*iter)->Translate(-prAxes.com); // Move to origin
    (*iter)->RotateUsingQuat(q);     // Rotate
    (*iter)->Translate(com);         // Move to new centre of mass
  }
}
