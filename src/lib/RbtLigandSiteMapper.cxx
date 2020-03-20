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

#include "RbtLigandSiteMapper.h"
#include "RbtFFTGrid.h"
#include "RbtMdlFileSource.h"
#include "RbtWorkSpace.h"

#include <functional>

// Static data member for class type
std::string RbtLigandSiteMapper::_CT("RbtLigandSiteMapper");
std::string RbtLigandSiteMapper::_REF_MOL("REF_MOL");
std::string RbtLigandSiteMapper::_VOL_INCR("VOL_INCR");
std::string RbtLigandSiteMapper::_SMALL_SPHERE("SMALL_SPHERE");
std::string RbtLigandSiteMapper::_GRIDSTEP("GRIDSTEP");
std::string RbtLigandSiteMapper::_RADIUS("RADIUS");
std::string RbtLigandSiteMapper::_MIN_VOLUME("MIN_VOLUME");
std::string RbtLigandSiteMapper::_MAX_CAVITIES("MAX_CAVITIES");

RbtLigandSiteMapper::RbtLigandSiteMapper(const std::string &strName)
    : RbtSiteMapper(_CT, strName) {
  // Add parameters
  AddParameter(_REF_MOL, "ref.sd");
  AddParameter(_VOL_INCR, 0.0);
  AddParameter(_SMALL_SPHERE, 1.5);
  AddParameter(_GRIDSTEP, 0.5);
  AddParameter(_RADIUS, 10.0);
  AddParameter(_MIN_VOLUME, 100.0); // Min cavity volume in A^3
  AddParameter(_MAX_CAVITIES, 99);  // Max number of cavities to return
#ifdef _DEBUG
  std::cout << _CT << " parameterised constructor" << std::endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtLigandSiteMapper::~RbtLigandSiteMapper() {
#ifdef _DEBUG
  std::cout << _CT << " destructor" << std::endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

RbtCavityList RbtLigandSiteMapper::operator()() {
  RbtCavityList cavityList;
  RbtModelPtr spReceptor(GetReceptor());
  if (spReceptor.Null())
    return cavityList;

  std::string strRef = GetParameter(_REF_MOL);
  double rVolIncr = GetParameter(_VOL_INCR);
  double smallR = GetParameter(_SMALL_SPHERE);
  double step = GetParameter(_GRIDSTEP);
  double radius = GetParameter(_RADIUS);
  double minVol = GetParameter(_MIN_VOLUME);
  unsigned int maxCavities = GetParameter(_MAX_CAVITIES);
  int iTrace = GetTrace();

  // Grid values
  const double recVal = -1.0; // Receptor volume
  const double cavVal = 1.0;  // Ligand expanded volume

  // Convert from min volume (in A^3) to min size (number of grid points)
  double minSize = minVol / (step * step * step);

  // Get the reference atom list from the mol/SD file
  // Simple extension would be to read all the records from the SD file
  // to construct a docking site from the superimposition of all the reference
  // ligands
  std::string strMolFile = Rbt::GetRbtFileName("data/ligands", strRef);
  RbtMolecularFileSourcePtr spMdlFileSource(
      new RbtMdlFileSource(strMolFile, false, false, true));
  // Only include non-H reference atoms in the mapping
  RbtAtomList refAtomList = Rbt::GetAtomList(spMdlFileSource->GetAtomList(),
                                             std::not1(Rbt::isAtomicNo_eq(1)));

  RbtFFTGridPtr spGrid;

  // Only include non-H receptor atoms in the mapping
  RbtAtomList atomList = Rbt::GetAtomList(spReceptor->GetAtomList(),
                                          std::not1(Rbt::isAtomicNo_eq(1)));
  RbtVector gridStep(step, step, step);

  // Construct the grid to cover the ligand coords + radius
  RbtCoordList refCoordList;
  Rbt::GetCoordList(refAtomList, refCoordList);
  double border = radius + smallR + step;
  RbtCoord minCoord = Rbt::Min(refCoordList) - border;
  RbtCoord maxCoord = Rbt::Max(refCoordList) + border;
  RbtVector extent = maxCoord - minCoord;
  Eigen::Vector3d nXYZ = extent.xyz.array() / gridStep.xyz.array();
  unsigned int nX = static_cast<unsigned int>(nXYZ(0)) + 1;
  unsigned int nY = static_cast<unsigned int>(nXYZ(1)) + 1;
  unsigned int nZ = static_cast<unsigned int>(nXYZ(2)) + 1;
  spGrid = RbtFFTGridPtr(new RbtFFTGrid(minCoord, gridStep, nX, nY, nZ));
  spGrid->SetAllValues(recVal);

  // Clear the spheres around each ligand atom
  for (RbtCoordListConstIter iter = refCoordList.begin();
       iter != refCoordList.end(); iter++) {
    spGrid->SetSphere((*iter), radius, 0.0, true);
  }

  if (iTrace > 1) {
    std::cout << std::endl << "INITIALISATION" << std::endl;
    std::cout << "Radius=" << radius << std::endl;
    std::cout << "Border=" << border << std::endl;
    std::cout << "N(excluded)=" << spGrid->Count(recVal) << std::endl;
    std::cout << "N(unallocated)=" << spGrid->Count(0.0) << std::endl;
  }

  // Now exclude the receptor volume
  for (RbtAtomListConstIter iter = atomList.begin(); iter != atomList.end();
       iter++) {
    double r = (**iter).GetVdwRadius();
    spGrid->SetSphere((**iter).GetCoords(), r + rVolIncr, recVal, true);
  }

  if (iTrace > 1) {
    std::cout << std::endl << "EXCLUDE RECEPTOR VOLUME" << std::endl;
    std::cout << "N(excluded)=" << spGrid->Count(recVal) << std::endl;
    std::cout << "N(unallocated)=" << spGrid->Count(0.0) << std::endl;
  }

  // Map with a small solvent sphere
  spGrid->SetAccessible(smallR, 0.0, recVal, cavVal, false);

  if (iTrace > 1) {
    std::cout << std::endl << "FINAL CAVITIES" << std::endl;
    std::cout << "N(excluded)=" << spGrid->Count(recVal) << std::endl;
    std::cout << "N(unallocated)=" << spGrid->Count(0.0) << std::endl;
    std::cout << "N(cavities)=" << spGrid->Count(cavVal) << std::endl;
    std::cout << std::endl << "Min cavity size=" << minSize << std::endl;
  }

  // Find the contiguous regions of cavity grid points
  RbtFFTPeakMap peakMap =
      spGrid->FindPeaks(cavVal, static_cast<unsigned int>(minSize));

  // Convert peaks to cavity format
  for (RbtFFTPeakMapConstIter pIter = peakMap.begin(); pIter != peakMap.end();
       pIter++) {
    RbtFFTPeakPtr spPeak((*pIter).second);
    RbtCoordList coordList = spGrid->GetCoordList(spPeak->points);
    RbtCavityPtr spCavity = RbtCavityPtr(new RbtCavity(coordList, gridStep));
    cavityList.push_back(spCavity);
  }

  // Sort cavities by volume
  std::sort(cavityList.begin(), cavityList.end(),
            Rbt::RbtCavityPtrCmp_Volume());

  if (iTrace > 0) {
    for (RbtCavityListConstIter cIter = cavityList.begin();
         cIter != cavityList.end(); cIter++) {
      std::cout << (**cIter) << std::endl;
    }
  }

  // Limit the number of cavities if necessary
  if (cavityList.size() > maxCavities) {
    if (iTrace > 0) {
      std::cout << std::endl
                << cavityList.size() << " cavities identified - limit to "
                << maxCavities << " largest cavities" << std::endl;
    }
    cavityList.erase(cavityList.begin() + maxCavities, cavityList.end());
  }

  return cavityList;
}
