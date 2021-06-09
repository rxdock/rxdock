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

#include "rxdock/RandLigTransform.h"

#include <loguru.hpp>

using namespace rxdock;

// Static data member for class type
const std::string RandLigTransform::_CT = "RandLigTransform";
// Parameter names
const std::string RandLigTransform::_TORS_STEP = "torsion-step";

////////////////////////////////////////
// Constructors/destructors
RandLigTransform::RandLigTransform(const std::string &strName)
    : BaseUniMolTransform(_CT, strName), m_rand(GetRandInstance()) {
  LOG_F(2, "RandLigTransform parameterised constructor");
  // Add parameters
  AddParameter(_TORS_STEP, 180);
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RandLigTransform::~RandLigTransform() {
  LOG_F(2, "RandLigTransform destructor");
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

////////////////////////////////////////
// Protected methods
///////////////////
void RandLigTransform::SetupTransform() {
  // Clear the rotable bond list from the previous model
  m_rotableBonds.clear();
  // Check for undefined ligand
  if (GetLigand().Null())
    return;
  m_rotableBonds = GetRotatableBondList(GetLigand()->GetBondList());
}

////////////////////////////////////////
// Private methods
///////////////////
// Pure virtual in BaseTransform
// Actually apply the transform
void RandLigTransform::Execute() {
  ModelPtr spLigand = GetLigand();
  if (spLigand.Null())
    return;
  double torsStep = GetParameter(_TORS_STEP);
  for (BondListIter iter = m_rotableBonds.begin(); iter != m_rotableBonds.end();
       iter++) {
    double thetaDeg = 2.0 * torsStep * m_rand.GetRandom01() - torsStep;
    spLigand->RotateBond(*iter, thetaDeg, false);
  }
}
