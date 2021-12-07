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

#include "rxdock/BaseIdxSF.h"
#include "rxdock/WorkSpace.h"

#include <loguru.hpp>

using namespace rxdock;

// Static data members
const std::string BaseIdxSF::_CT = "BaseIdxSF";
const std::string BaseIdxSF::_GRIDSTEP = "grid-step";
const std::string BaseIdxSF::_BORDER = "border";

BaseIdxSF::BaseIdxSF() : m_gridStep(0.5), m_border(1.0) {
  LOG_F(2, "BaseIdxSF default constructor");
  // Add parameters
  AddParameter(_GRIDSTEP, m_gridStep);
  AddParameter(_BORDER, m_border);
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

BaseIdxSF::~BaseIdxSF() { _RBTOBJECTCOUNTER_DESTR_(_CT); }

double BaseIdxSF::GetGridStep() const { return m_gridStep; }

void BaseIdxSF::SetGridStep(double step) { SetParameter(_GRIDSTEP, step); }

double BaseIdxSF::GetBorder() const { return m_border; }

void BaseIdxSF::SetBorder(double border) { SetParameter(_BORDER, border); }

// DM 10 Apr 2002
// I know, I know, grids should be templated to avoid the need for two different
// CreateGrid methods...
InteractionGridPtr BaseIdxSF::CreateInteractionGrid() const {
  // Create a grid covering the docking site
  DockingSitePtr spDS = GetWorkSpace()->GetDockingSite();
  if (spDS.Null())
    return InteractionGridPtr();

  // Extend grid by _BORDER, mainly to allow for the possibility of polar
  // hydrogens falling outside of the docking site
  Coord minCoord = spDS->GetMinCoord() - m_border;
  Coord maxCoord = spDS->GetMaxCoord() + m_border;
  Vector recepExtent = maxCoord - minCoord;
  Vector gridStep(m_gridStep, m_gridStep, m_gridStep);
  Eigen::Vector3d nXYZ = recepExtent.xyz.array() / gridStep.xyz.array();
  unsigned int nX = static_cast<unsigned int>(nXYZ(0)) + 1;
  unsigned int nY = static_cast<unsigned int>(nXYZ(1)) + 1;
  unsigned int nZ = static_cast<unsigned int>(nXYZ(2)) + 1;
  return InteractionGridPtr(
      new InteractionGrid(minCoord, gridStep, nX, nY, nZ));
}

NonBondedGridPtr BaseIdxSF::CreateNonBondedGrid() const {
  // Create a grid covering the docking site
  DockingSitePtr spDS = GetWorkSpace()->GetDockingSite();
  if (spDS.Null())
    return InteractionGridPtr();

  // Extend grid by _BORDER, mainly to allow for the possibility of polar
  // hydrogens falling outside of the docking site
  Coord minCoord = spDS->GetMinCoord() - m_border;
  Coord maxCoord = spDS->GetMaxCoord() + m_border;
  Vector recepExtent = maxCoord - minCoord;
  Vector gridStep(m_gridStep, m_gridStep, m_gridStep);
  Eigen::Vector3d nXYZ = recepExtent.xyz.array() / gridStep.xyz.array();
  unsigned int nX = static_cast<unsigned int>(nXYZ(0)) + 1;
  unsigned int nY = static_cast<unsigned int>(nXYZ(1)) + 1;
  unsigned int nZ = static_cast<unsigned int>(nXYZ(2)) + 1;
  return NonBondedGridPtr(new NonBondedGrid(minCoord, gridStep, nX, nY, nZ));
}

NonBondedHHSGridPtr BaseIdxSF::CreateNonBondedHHSGrid() const {
  // Create a grid covering the docking site
  DockingSitePtr spDS = GetWorkSpace()->GetDockingSite();
  if (spDS.Null())
    return NonBondedHHSGridPtr();

  // Extend grid by _BORDER, mainly to allow for the possibility of polar
  // hydrogens falling outside of the docking site
  Coord minCoord = spDS->GetMinCoord() - m_border;
  Coord maxCoord = spDS->GetMaxCoord() + m_border;
  Vector recepExtent = maxCoord - minCoord;
  Vector gridStep(m_gridStep, m_gridStep, m_gridStep);
  Eigen::Vector3d nXYZ = recepExtent.xyz.array() / gridStep.xyz.array();
  unsigned int nX = static_cast<unsigned int>(nXYZ(0)) + 1;
  unsigned int nY = static_cast<unsigned int>(nXYZ(1)) + 1;
  unsigned int nZ = static_cast<unsigned int>(nXYZ(2)) + 1;
  return NonBondedHHSGridPtr(
      new NonBondedHHSGrid(minCoord, gridStep, nX, nY, nZ));
}

double BaseIdxSF::GetMaxError() const {
  // maxError is half a grid diagonal. This is the tolerance we have to allow
  // when indexing the receptor atoms on the grid When retrieving the nearest
  // neighbours to a ligand atom, we do the lookup on the nearest grid point to
  // the ligand atom. However the actual ligand atom - receptor atom distance
  // may be maxError Angstroms closer than the grid point - receptor atom
  // distance used in the indexing.
  return 0.5 * std::sqrt(3.0) * m_gridStep;
}

// DM 12 Apr 2002
// Returns the maximum range of the scoring function,
// corrected for max grid error, and grid border around docking site
// This should be used by subclasses for selecting the receptor atoms to index
// GetCorrectedRange() = GetRange() + GetMaxError() + GetBorder()
double BaseIdxSF::GetCorrectedRange() const {
  return GetRange() + GetMaxError() + m_border;
}

// As this has a virtual base class we need a separate OwnParameterUpdated
// which can be called by concrete subclass ParameterUpdated methods
// See Stroustrup C++ 3rd edition, p395, on programming virtual base classes
void BaseIdxSF::OwnParameterUpdated(const std::string &strName) {
  if (strName == _GRIDSTEP) {
    m_gridStep = GetParameter(_GRIDSTEP);
  } else if (strName == _BORDER) {
    m_border = GetParameter(_BORDER);
  }
}
