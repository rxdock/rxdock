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

#include "rxdock/AlignTransform.h"
#include "rxdock/DockingSite.h"
#include "rxdock/WorkSpace.h"

#include <fmt/ostream.h>
#include <loguru.hpp>

using namespace rxdock;

// Static data member for class type
const std::string AlignTransform::_CT = "AlignTransform";
// Parameter names
const std::string AlignTransform::_COM = "center-of-mass";
const std::string AlignTransform::_AXES = "axes";

////////////////////////////////////////
// Constructors/destructors
AlignTransform::AlignTransform(const std::string &strName)
    : BaseBiMolTransform(_CT, strName), m_rand(GetRandInstance()),
      m_totalSize(0) {
  LOG_F(2, "AlignTransform parameterised constructor");
  // Add parameters
  AddParameter(_COM, "ALIGN");
  AddParameter(_AXES, "ALIGN");
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

AlignTransform::~AlignTransform() {
  LOG_F(2, "AlignTransform destructor");
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

////////////////////////////////////////
// Protected methods
///////////////////
void AlignTransform::SetupReceptor() {
  m_cavities.clear();
  m_cumulSize.clear();
  m_totalSize = 0;
  // Return now if docking site not registered with workspace
  DockingSitePtr spDS = GetWorkSpace()->GetDockingSite();
  if (spDS.Null())
    return;

  m_cavities = spDS->GetCavityList();
  if (m_cavities.empty())
    return;

  // Determine cumulative sizes, for selecting cavities with a probability
  // proportional to their size
  for (CavityListConstIter cIter = m_cavities.begin();
       cIter != m_cavities.end(); cIter++) {
    m_totalSize += (*cIter)->GetNumCoords();
    m_cumulSize.push_back(m_totalSize);
  }
}

void AlignTransform::SetupLigand() {}    // No ligand setup required
void AlignTransform::SetupTransform() {} // No further setup required

////////////////////////////////////////
// Private methods
///////////////////
// Pure virtual in BaseTransform
// Actually apply the transform
void AlignTransform::Execute() {
  ModelPtr spLigand(GetLigand());
  if (spLigand.Null() || m_cavities.empty())
    return;

  // Select a cavity at random, weighted by each cavity size
  int iRnd = m_rand.GetRandomInt(m_totalSize);
  int iCavity(0);
  for (iCavity = 0; iRnd >= m_cumulSize[iCavity]; iCavity++)
    ;
  CavityPtr spCavity = m_cavities[iCavity];

  // Get the coord list and principal axes for the cavity
  const CoordList &coordList = spCavity->GetCoordList();
  const PrincipalAxes &prAxes = spCavity->GetPrincipalAxes();

  // Get the alignment parameters
  std::string strPlaceCOM = GetParameter(_COM);
  std::string strPlaceAxes = GetParameter(_AXES);

  // 1. Ligand translation
  // A. Random
  if ((strPlaceCOM == "RANDOM") && !coordList.empty()) {
    // Select a coord at random
    int iRand = m_rand.GetRandomInt(coordList.size());
    Coord asCavityCoord = coordList[iRand];
    LOG_F(1, "Translating ligand COM to active site coord #{}: {}", iRand,
          asCavityCoord);
    // Translate the ligand center of mass to the selected coord
    spLigand->SetCenterOfMass(asCavityCoord);
  }
  // B. Active site center of mass
  else if (strPlaceCOM == "ALIGN") {
    LOG_F(1, "Translating ligand COM to active site COM: {}", prAxes.com);
    spLigand->SetCenterOfMass(prAxes.com);
  }

  // 2. Ligand axes
  // A. Random rotation around random axis
  if (strPlaceAxes == "RANDOM") {
    double thetaDeg = 180.0 * m_rand.GetRandom01();
    Coord axis = m_rand.GetRandomUnitVector();
    LOG_F(1, "Rotating ligand by {} deg around axis={} through COM", thetaDeg,
          axis);
    spLigand->Rotate(axis, thetaDeg);
  }
  // B. Align ligand principal axes with principal axes of active site
  else if (strPlaceAxes == "ALIGN") {
    spLigand->AlignPrincipalAxes(
        prAxes,
        false); // false = don't translate COM as we've already done it above
    LOG_F(1, "Aligning ligand principal axes with active site principal axes");
    // Make random 180 deg rotations around each of the principal axes
    if (m_rand.GetRandom01() < 0.5) {
      spLigand->Rotate(prAxes.axis1, 180.0, prAxes.com);
      LOG_F(1, "180 deg rotation around PA#1");
    }
    if (m_rand.GetRandom01() < 0.5) {
      spLigand->Rotate(prAxes.axis2, 180.0, prAxes.com);
      LOG_F(1, "180 deg rotation around PA#2");
    }
    if (m_rand.GetRandom01() < 0.5) {
      spLigand->Rotate(prAxes.axis3, 180.0, prAxes.com);
      LOG_F(1, "180 deg rotation around PA#3");
    }
  }
}
