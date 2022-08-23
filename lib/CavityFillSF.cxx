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

#include "rxdock/CavityFillSF.h"
#include "rxdock/DockingSite.h"
#include "rxdock/WorkSpace.h"

#include <fmt/ostream.h>
#include <loguru.hpp>

using namespace rxdock;

// Static data members
const std::string CavityFillSF::_CT = "CavityFillSF";

// NB - Virtual base class constructor (BaseSF) gets called first,
// implicit constructor for BaseInterSF is called second
CavityFillSF::CavityFillSF(const std::string &strName) : BaseSF(_CT, strName) {
  // Add parameters
  LOG_F(2, "CavityFillSF parameterised constructor");
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

CavityFillSF::~CavityFillSF() {
  LOG_F(2, "CavityFillSF destructor");
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

void CavityFillSF::SetupReceptor() {
  m_spGrid = FFTGridPtr();
  DockingSitePtr spDS = GetWorkSpace()->GetDockingSite();
  if (spDS.Null())
    return;

  // Recreate the cavity grid
  CavityList cavList = spDS->GetCavityList();
  if (cavList.empty())
    return;
  Vector gridStep = cavList.front()->GetGridStep();
  double border = 10.0;
  Coord minCoord = spDS->GetMinCoord() - border;
  Coord maxCoord = spDS->GetMaxCoord() + border;
  Vector recepExtent = maxCoord - minCoord;
  Eigen::Vector3d nXYZ = recepExtent.xyz.array() / gridStep.xyz.array();
  unsigned int nX = static_cast<unsigned int>(nXYZ(0)) + 1;
  unsigned int nY = static_cast<unsigned int>(nXYZ(1)) + 1;
  unsigned int nZ = static_cast<unsigned int>(nXYZ(2)) + 1;
  m_spGrid = new FFTGrid(minCoord, gridStep, nX, nY, nZ);

  AtomList atomList = GetReceptor()->GetAtomList();
  for (AtomListConstIter iter = atomList.begin(); iter != atomList.end();
       iter++) {
    double r = (**iter).GetVdwRadius();
    m_spGrid->SetSphere((**iter).GetCoords(), r + 0.3, -1.0, true);
  }

  LOG_F(INFO, "EXCLUDE RECEPTOR VOLUME");
  LOG_F(INFO, "N(excluded)={}", m_spGrid->Count(-1.0));
  LOG_F(INFO, "N(unallocated)={}", m_spGrid->Count(0.0));

  for (CavityListConstIter iter = cavList.begin(); iter != cavList.end();
       ++iter) {
    const CoordList coordList = (*iter)->GetCoordList();
    for (CoordListConstIter iter2 = coordList.begin(); iter2 != coordList.end();
         ++iter2) {
      m_spGrid->SetValue(*iter2, 1.0);
    }
  }

  LOG_F(INFO, "DEFINE CAVITY VOLUME");
  LOG_F(INFO, "N(excluded)={}", m_spGrid->Count(-1.0));
  LOG_F(INFO, "N(cavity)={}", m_spGrid->Count(1.0));
  LOG_F(INFO, "N(unallocated)={}", m_spGrid->Count(0.0));

  // double largeR = 4.0;
  // Map with a large solvent sphere
  // m_spGrid->SetAccessible(largeR, 1.0, -1.0, 0.0, false);
  // LOG_F(INFO, "MAP WITH LARGE SPHERE");
  // LOG_F(INFO, "N(excluded)={}", m_spGrid->Count(-1.0));
  // LOG_F(INFO, "N(cavity)={}", m_spGrid->Count(1.0));
  // LOG_F(INFO, "N(unallocated)={}", m_spGrid->Count(0.0));
}

void CavityFillSF::SetupLigand() {
  m_ligAtomList.clear();
  if (GetLigand().Null())
    return;

  m_ligAtomList = GetLigand()->GetAtomList();
}

void CavityFillSF::SetupScore() {}

double CavityFillSF::RawScore() const {
  // Check grid is defined
  if (m_spGrid.Null())
    return 0.0;

  FFTGridPtr gridCopy = new FFTGrid(*m_spGrid);

  for (AtomListConstIter iter = m_ligAtomList.begin();
       iter != m_ligAtomList.end(); ++iter) {
    double r = (*iter)->GetVdwRadius();
    gridCopy->SetSphere((*iter)->GetCoords(), r + 0.3, -1.0, true);
  }

  LOG_F(INFO, "EXCLUDE LIGAND VOLUME");
  LOG_F(INFO, "N(excluded)={}", gridCopy->Count(-1.0));
  LOG_F(INFO, "N(cavity)={}", gridCopy->Count(1.0));
  LOG_F(INFO, "N(unallocated)={}", gridCopy->Count(0.0));

  // Map with a small solvent sphere
  gridCopy->SetAccessible(2.0, 1.0, -1.0, 0.0, false);
  gridCopy->SetAccessible(1.5, 1.0, -1.0, 2.0, false);

  LOG_F(INFO, "VOID DETECTION");
  LOG_F(INFO, "N(excluded)={}", gridCopy->Count(-1.0));
  LOG_F(INFO, "N(cavity)={}", gridCopy->Count(1.0));
  LOG_F(INFO, "N(voids)={}", gridCopy->Count(2.0));
  LOG_F(INFO, "N(unallocated)={}", gridCopy->Count(0.0));

  gridCopy->ReplaceValue(-1.0, 0.0);
  gridCopy->ReplaceValue(1.0, 0.0);
  std::ofstream dumpFile("void.grd");
  gridCopy->PrintInsightGrid(dumpFile);
  dumpFile.close();

  // Find the contiguous regions of cavity grid points
  unsigned int minSize = 50;
  FFTPeakMap peakMap = gridCopy->FindPeaks(2.0, minSize);

  // Convert peaks to cavity format
  CavityList cavityList;
  Vector gridStep = gridCopy->GetGridStep();
  for (FFTPeakMapConstIter pIter = peakMap.begin(); pIter != peakMap.end();
       pIter++) {
    FFTPeakPtr spPeak((*pIter).second);
    CoordList coordList = gridCopy->GetCoordList(spPeak->points);
    CavityPtr spCavity = CavityPtr(new Cavity(coordList, gridStep));
    cavityList.push_back(spCavity);
  }

  // Sort cavities by volume
  std::sort(cavityList.begin(), cavityList.end(), CavityPtrCmp_Volume());

  LOG_F(1, "Printing {} cavities sorted by volume", cavityList.size());
  for (CavityListConstIter cIter = cavityList.begin();
       cIter != cavityList.end(); cIter++) {
    LOG_F(1, "{}", (**cIter));
  }

  return 0.0;
}

// DM 25 Oct 2000 - track changes to parameter values in local data members
// ParameterUpdated is invoked by ParamHandler::SetParameter
void CavityFillSF::ParameterUpdated(const std::string &strName) {
  BaseSF::ParameterUpdated(strName);
}

void rxdock::to_json(json &j, const CavityFillSF &sf) {
  json atomList;
  for (const auto &cIter : sf.m_ligAtomList) {
    json atom = *cIter;
    atomList.push_back(atom);
  }

  j = json{{"atoms", atomList}};
  if (atomList.size() > 0 && !sf.m_spGrid.Null()) {
    j["fft-grid"] = *sf.m_spGrid;
  }
}

void rxdock::from_json(const json &j, CavityFillSF &sf) {
  sf.m_ligAtomList.clear();
  sf.m_ligAtomList.reserve(j.at("atoms").size());
  for (auto &atom : j.at("atoms")) {
    AtomPtr spAtom = AtomPtr(new Atom(atom));
    sf.m_ligAtomList.push_back(spAtom);
  }
  if (sf.m_ligAtomList.size() > 0) {
    sf.m_spGrid = FFTGridPtr(new FFTGrid(j.at("fft-grid")));
  }
}
