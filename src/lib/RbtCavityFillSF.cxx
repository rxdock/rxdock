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

#include "RbtCavityFillSF.h"
#include "RbtDockingSite.h"
#include "RbtWorkSpace.h"

using namespace rxdock;

// Static data members
std::string RbtCavityFillSF::_CT("RbtCavityFillSF");

// NB - Virtual base class constructor (RbtBaseSF) gets called first,
// implicit constructor for RbtBaseInterSF is called second
RbtCavityFillSF::RbtCavityFillSF(const std::string &strName)
    : RbtBaseSF(_CT, strName) {
  // Add parameters
#ifdef _DEBUG
  std::cout << _CT << " parameterised constructor" << std::endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtCavityFillSF::~RbtCavityFillSF() {
#ifdef _DEBUG
  std::cout << _CT << " destructor" << std::endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

void RbtCavityFillSF::SetupReceptor() {
  m_spGrid = RbtFFTGridPtr();
  RbtDockingSitePtr spDS = GetWorkSpace()->GetDockingSite();
  if (spDS.Null())
    return;

  int iTrace = GetTrace();

  // Recreate the cavity grid
  RbtCavityList cavList = spDS->GetCavityList();
  if (cavList.empty())
    return;
  RbtVector gridStep = cavList.front()->GetGridStep();
  double border = 10.0;
  RbtCoord minCoord = spDS->GetMinCoord() - border;
  RbtCoord maxCoord = spDS->GetMaxCoord() + border;
  RbtVector recepExtent = maxCoord - minCoord;
  Eigen::Vector3d nXYZ = recepExtent.xyz.array() / gridStep.xyz.array();
  unsigned int nX = static_cast<unsigned int>(nXYZ(0)) + 1;
  unsigned int nY = static_cast<unsigned int>(nXYZ(1)) + 1;
  unsigned int nZ = static_cast<unsigned int>(nXYZ(2)) + 1;
  m_spGrid = new RbtFFTGrid(minCoord, gridStep, nX, nY, nZ);

  RbtAtomList atomList = GetReceptor()->GetAtomList();
  for (RbtAtomListConstIter iter = atomList.begin(); iter != atomList.end();
       iter++) {
    double r = (**iter).GetVdwRadius();
    m_spGrid->SetSphere((**iter).GetCoords(), r + 0.3, -1.0, true);
  }

  if (iTrace > 1) {
    std::cout << std::endl << "EXCLUDE RECEPTOR VOLUME" << std::endl;
    std::cout << "N(excluded)=" << m_spGrid->Count(-1.0) << std::endl;
    std::cout << "N(unallocated)=" << m_spGrid->Count(0.0) << std::endl;
  }

  for (RbtCavityListConstIter iter = cavList.begin(); iter != cavList.end();
       ++iter) {
    const RbtCoordList coordList = (*iter)->GetCoordList();
    for (RbtCoordListConstIter iter2 = coordList.begin();
         iter2 != coordList.end(); ++iter2) {
      m_spGrid->SetValue(*iter2, 1.0);
    }
  }

  if (iTrace > 1) {
    std::cout << std::endl << "DEFINE CAVITY VOLUME" << std::endl;
    std::cout << "N(excluded)=" << m_spGrid->Count(-1.0) << std::endl;
    std::cout << "N(cavity)=" << m_spGrid->Count(1.0) << std::endl;
    std::cout << "N(unallocated)=" << m_spGrid->Count(0.0) << std::endl;
  }

  // RbtDouble largeR = 4.0;
  // Map with a large solvent sphere
  // m_spGrid->SetAccessible(largeR,1.0,-1.0,0.0,false);
  // if (iTrace > 1) {
  //  std::cout << std::endl << "MAP WITH LARGE SPHERE" << std::endl;
  //  std::cout << "N(excluded)=" << m_spGrid->Count(-1.0) << std::endl;
  //  std::cout << "N(cavity)=" << m_spGrid->Count(1.0) << std::endl;
  //  std::cout << "N(unallocated)=" << m_spGrid->Count(0.0) << std::endl;
  //}
}

void RbtCavityFillSF::SetupLigand() {
  m_ligAtomList.clear();
  if (GetLigand().Null())
    return;

  m_ligAtomList = GetLigand()->GetAtomList();
}

void RbtCavityFillSF::SetupScore() {}

double RbtCavityFillSF::RawScore() const {
  // Check grid is defined
  if (m_spGrid.Null())
    return 0.0;

  int iTrace = GetTrace();
  RbtFFTGridPtr gridCopy = new RbtFFTGrid(*m_spGrid);

  for (RbtAtomListConstIter iter = m_ligAtomList.begin();
       iter != m_ligAtomList.end(); ++iter) {
    double r = (*iter)->GetVdwRadius();
    gridCopy->SetSphere((*iter)->GetCoords(), r + 0.3, -1.0, true);
  }

  if (iTrace > 1) {
    std::cout << std::endl << "EXCLUDE LIGAND VOLUME" << std::endl;
    std::cout << "N(excluded)=" << gridCopy->Count(-1.0) << std::endl;
    std::cout << "N(cavity)=" << gridCopy->Count(1.0) << std::endl;
    std::cout << "N(unallocated)=" << gridCopy->Count(0.0) << std::endl;
  }

  // Map with a small solvent sphere
  gridCopy->SetAccessible(2.0, 1.0, -1.0, 0.0, false);
  gridCopy->SetAccessible(1.5, 1.0, -1.0, 2.0, false);

  if (iTrace > 1) {
    std::cout << std::endl << "VOID DETECTION" << std::endl;
    std::cout << "N(excluded)=" << gridCopy->Count(-1.0) << std::endl;
    std::cout << "N(cavity)=" << gridCopy->Count(1.0) << std::endl;
    std::cout << "N(voids)=" << gridCopy->Count(2.0) << std::endl;
    std::cout << "N(unallocated)=" << gridCopy->Count(0.0) << std::endl;
  }

  gridCopy->ReplaceValue(-1.0, 0.0);
  gridCopy->ReplaceValue(1.0, 0.0);
  std::ofstream dumpFile("void.grd");
  gridCopy->PrintInsightGrid(dumpFile);
  dumpFile.close();

  // Find the contiguous regions of cavity grid points
  unsigned int minSize = 50;
  RbtFFTPeakMap peakMap = gridCopy->FindPeaks(2.0, minSize);

  // Convert peaks to cavity format
  RbtCavityList cavityList;
  RbtVector gridStep = gridCopy->GetGridStep();
  for (RbtFFTPeakMapConstIter pIter = peakMap.begin(); pIter != peakMap.end();
       pIter++) {
    RbtFFTPeakPtr spPeak((*pIter).second);
    RbtCoordList coordList = gridCopy->GetCoordList(spPeak->points);
    RbtCavityPtr spCavity = RbtCavityPtr(new RbtCavity(coordList, gridStep));
    cavityList.push_back(spCavity);
  }

  // Sort cavities by volume
  std::sort(cavityList.begin(), cavityList.end(), RbtCavityPtrCmp_Volume());

  if (iTrace > 0) {
    for (RbtCavityListConstIter cIter = cavityList.begin();
         cIter != cavityList.end(); cIter++) {
      std::cout << (**cIter) << std::endl;
    }
  }

  return 0.0;
}

// DM 25 Oct 2000 - track changes to parameter values in local data members
// ParameterUpdated is invoked by RbtParamHandler::SetParameter
void RbtCavityFillSF::ParameterUpdated(const std::string &strName) {
  RbtBaseSF::ParameterUpdated(strName);
}
