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

#include "AlignTransform.h"
#include "DockingSite.h"
#include "WorkSpace.h"

using namespace rxdock;

// Static data member for class type
std::string AlignTransform::_CT("AlignTransform");
// Parameter names
std::string AlignTransform::_COM("COM");
std::string AlignTransform::_AXES("AXES");

////////////////////////////////////////
// Constructors/destructors
AlignTransform::AlignTransform(const std::string &strName)
    : BaseBiMolTransform(_CT, strName), m_rand(GetRandInstance()),
      m_totalSize(0) {
  // Add parameters
  AddParameter(_COM, "ALIGN");
  AddParameter(_AXES, "ALIGN");
#ifdef _DEBUG
  std::cout << _CT << " parameterised constructor" << std::endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

AlignTransform::~AlignTransform() {
#ifdef _DEBUG
  std::cout << _CT << " destructor" << std::endl;
#endif //_DEBUG
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

  int iTrace = GetTrace();

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
    if (iTrace > 1) {
      std::cout << "Translating ligand COM to active site coord #" << iRand
                << ": " << asCavityCoord << std::endl;
    }
    // Translate the ligand center of mass to the selected coord
    spLigand->SetCenterOfMass(asCavityCoord);
  }
  // B. Active site center of mass
  else if (strPlaceCOM == "ALIGN") {
    if (iTrace > 1) {
      std::cout << "Translating ligand COM to active site COM: " << prAxes.com
                << std::endl;
    }
    spLigand->SetCenterOfMass(prAxes.com);
  }

  // 2. Ligand axes
  // A. Random rotation around random axis
  if (strPlaceAxes == "RANDOM") {
    double thetaDeg = 180.0 * m_rand.GetRandom01();
    Coord axis = m_rand.GetRandomUnitVector();
    if (iTrace > 1) {
      std::cout << "Rotating ligand by " << thetaDeg
                << " deg around axis=" << axis << " through COM" << std::endl;
    }
    spLigand->Rotate(axis, thetaDeg);
  }
  // B. Align ligand principal axes with principal axes of active site
  else if (strPlaceAxes == "ALIGN") {
    spLigand->AlignPrincipalAxes(
        prAxes,
        false); // false = don't translate COM as we've already done it above
    if (iTrace > 1) {
      std::cout
          << "Aligning ligand principal axes with active site principal axes"
          << std::endl;
    }
    // Make random 180 deg rotations around each of the principal axes
    if (m_rand.GetRandom01() < 0.5) {
      spLigand->Rotate(prAxes.axis1, 180.0, prAxes.com);
      if (iTrace > 1) {
        std::cout << "180 deg rotation around PA#1" << std::endl;
      }
    }
    if (m_rand.GetRandom01() < 0.5) {
      spLigand->Rotate(prAxes.axis2, 180.0, prAxes.com);
      if (iTrace > 1) {
        std::cout << "180 deg rotation around PA#2" << std::endl;
      }
    }
    if (m_rand.GetRandom01() < 0.5) {
      spLigand->Rotate(prAxes.axis3, 180.0, prAxes.com);
      if (iTrace > 1) {
        std::cout << "180 deg rotation around PA#3" << std::endl;
      }
    }
  }
}
