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

#include "LigandSiteMapper.h"
#include "FFTGrid.h"
#include "MdlFileSource.h"
#include "WorkSpace.h"

#include <functional>

using namespace rxdock;

// Static data member for class type
std::string LigandSiteMapper::_CT("LigandSiteMapper");
std::string LigandSiteMapper::_REF_MOL("REF_MOL");
std::string LigandSiteMapper::_VOL_INCR("VOL_INCR");
std::string LigandSiteMapper::_SMALL_SPHERE("SMALL_SPHERE");
std::string LigandSiteMapper::_GRIDSTEP("GRIDSTEP");
std::string LigandSiteMapper::_RADIUS("RADIUS");
std::string LigandSiteMapper::_MIN_VOLUME("MIN_VOLUME");
std::string LigandSiteMapper::_MAX_CAVITIES("MAX_CAVITIES");

LigandSiteMapper::LigandSiteMapper(const std::string &strName)
    : SiteMapper(_CT, strName) {
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

LigandSiteMapper::~LigandSiteMapper() {
#ifdef _DEBUG
  std::cout << _CT << " destructor" << std::endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

CavityList LigandSiteMapper::operator()() {
  CavityList cavityList;
  ModelPtr spReceptor(GetReceptor());
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
  std::string strMolFile = GetDataFileName("data/ligands", strRef);
  MolecularFileSourcePtr spMdlFileSource(
      new MdlFileSource(strMolFile, false, false, true));
  // Only include non-H reference atoms in the mapping
  AtomList refAtomList = GetAtomListWithPredicate(
      spMdlFileSource->GetAtomList(), std::not1(isAtomicNo_eq(1)));

  FFTGridPtr spGrid;

  // Only include non-H receptor atoms in the mapping
  AtomList atomList = GetAtomListWithPredicate(spReceptor->GetAtomList(),
                                               std::not1(isAtomicNo_eq(1)));
  Vector gridStep(step, step, step);

  // Construct the grid to cover the ligand coords + radius
  CoordList refCoordList;
  GetCoordList(refAtomList, refCoordList);
  double border = radius + smallR + step;
  Coord minCoord = Min(refCoordList) - border;
  Coord maxCoord = Max(refCoordList) + border;
  Vector extent = maxCoord - minCoord;
  Eigen::Vector3d nXYZ = extent.xyz.array() / gridStep.xyz.array();
  unsigned int nX = static_cast<unsigned int>(nXYZ(0)) + 1;
  unsigned int nY = static_cast<unsigned int>(nXYZ(1)) + 1;
  unsigned int nZ = static_cast<unsigned int>(nXYZ(2)) + 1;
  spGrid = FFTGridPtr(new FFTGrid(minCoord, gridStep, nX, nY, nZ));
  spGrid->SetAllValues(recVal);

  // Clear the spheres around each ligand atom
  for (CoordListConstIter iter = refCoordList.begin();
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
  for (AtomListConstIter iter = atomList.begin(); iter != atomList.end();
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
  FFTPeakMap peakMap =
      spGrid->FindPeaks(cavVal, static_cast<unsigned int>(minSize));

  // Convert peaks to cavity format
  for (FFTPeakMapConstIter pIter = peakMap.begin(); pIter != peakMap.end();
       pIter++) {
    FFTPeakPtr spPeak((*pIter).second);
    CoordList coordList = spGrid->GetCoordList(spPeak->points);
    CavityPtr spCavity = CavityPtr(new Cavity(coordList, gridStep));
    cavityList.push_back(spCavity);
  }

  // Sort cavities by volume
  std::sort(cavityList.begin(), cavityList.end(), CavityPtrCmp_Volume());

  if (iTrace > 0) {
    for (CavityListConstIter cIter = cavityList.begin();
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
