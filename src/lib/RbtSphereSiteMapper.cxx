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

#include "RbtSphereSiteMapper.h"
#include "RbtFFTGrid.h"

#include <functional>

// Static data member for class type
std::string RbtSphereSiteMapper::_CT("RbtSphereSiteMapper");
std::string RbtSphereSiteMapper::_VOL_INCR("VOL_INCR");
std::string RbtSphereSiteMapper::_SMALL_SPHERE("SMALL_SPHERE");
std::string RbtSphereSiteMapper::_LARGE_SPHERE("LARGE_SPHERE");
std::string RbtSphereSiteMapper::_GRIDSTEP("GRIDSTEP");
std::string RbtSphereSiteMapper::_CENTER("CENTER");
std::string RbtSphereSiteMapper::_RADIUS("RADIUS");
std::string RbtSphereSiteMapper::_MIN_VOLUME("MIN_VOLUME");
std::string RbtSphereSiteMapper::_MAX_CAVITIES("MAX_CAVITIES");

RbtSphereSiteMapper::RbtSphereSiteMapper(const std::string &strName)
    : RbtSiteMapper(_CT, strName) {
  // Add parameters
  AddParameter(_VOL_INCR, 0.0);
  AddParameter(_SMALL_SPHERE, 1.5);
  AddParameter(_LARGE_SPHERE, 4.0);
  AddParameter(_GRIDSTEP, 0.5);
  AddParameter(_CENTER, RbtCoord());
  AddParameter(_RADIUS, 10.0);
  AddParameter(_MIN_VOLUME, 100.0); // Min cavity volume in A^3
  AddParameter(_MAX_CAVITIES, 99);  // Max number of cavities to return
#ifdef _DEBUG
  std::cout << _CT << " parameterised constructor" << std::endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtSphereSiteMapper::~RbtSphereSiteMapper() {
#ifdef _DEBUG
  std::cout << _CT << " destructor" << std::endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

RbtCavityList RbtSphereSiteMapper::operator()() {
  RbtCavityList cavityList;
  RbtModelPtr spReceptor(GetReceptor());
  if (spReceptor.Null())
    return cavityList;

  double rVolIncr = GetParameter(_VOL_INCR);
  double smallR = GetParameter(_SMALL_SPHERE);
  double largeR = GetParameter(_LARGE_SPHERE);
  double step = GetParameter(_GRIDSTEP);
  RbtCoord center = GetParameter(_CENTER);
  double radius = GetParameter(_RADIUS);
  double minVol = GetParameter(_MIN_VOLUME);
  unsigned int maxCavities = GetParameter(_MAX_CAVITIES);
  int iTrace = GetTrace();

  // Grid values
  const double recVal = -1.0;  // Receptor volume
  const double larVal = -0.75; // Accessible to large sphere
  const double excVal = -0.5;  // Excluded from calculation
  const double borVal =
      -0.25; // Border region (used for mapping large sphere only)
  const double cavVal = 1.0; // Cavities

  // Convert from min volume (in A^3) to min size (number of grid points)
  double minSize = minVol / (step * step * step);
  RbtFFTGridPtr spReceptorGrid;
  // Only include non-H receptor atoms in the mapping
  RbtAtomList atomList = Rbt::GetAtomList(spReceptor->GetAtomList(),
                                          std::not1(Rbt::isAtomicNo_eq(1)));
  RbtVector gridStep(step, step, step);

  // We extend the grid by 2*largeR on each side to eliminate edge effects in
  // the cavity mapping
  RbtCoord minCoord = center - radius;
  RbtCoord maxCoord = center + radius;
  double border = 2.0 * (largeR + step);
  minCoord -= border;
  maxCoord += border;
  RbtVector recepExtent = maxCoord - minCoord;
  unsigned int nX = int(recepExtent.x / gridStep.x) + 1;
  unsigned int nY = int(recepExtent.y / gridStep.y) + 1;
  unsigned int nZ = int(recepExtent.z / gridStep.z) + 1;
  spReceptorGrid =
      RbtFFTGridPtr(new RbtFFTGrid(minCoord, gridStep, nX, nY, nZ));
  center = spReceptorGrid->GetGridCenter();

  // Initialise the grid with a zero-value region in the user-specified sphere
  // surrounded by a border region of thickness = large sphere radius
  spReceptorGrid->SetAllValues(excVal);
  spReceptorGrid->SetSphere(center, radius + largeR, borVal, true);
  spReceptorGrid->SetSphere(center, radius, 0.0, true);
  if (iTrace > 1) {
    std::cout << std::endl << "INITIALISATION" << std::endl;
    std::cout << "Center=" << center << std::endl;
    std::cout << "Radius=" << radius << std::endl;
    std::cout << "Border=" << border << std::endl;
    std::cout << "N(excluded)=" << spReceptorGrid->Count(excVal) << std::endl;
    std::cout << "N(border)=" << spReceptorGrid->Count(borVal) << std::endl;
    std::cout << "N(unallocated)=" << spReceptorGrid->Count(0.0) << std::endl;
  }

  // Set all vdW volume grid points to the value -1
  // Add increment to all vdw radii
  // Iterate over all receptor atoms as we want to include those whose centers
  // are outside the active site, but whose vdw volumes overlap the active site
  // region.
  for (RbtAtomListConstIter iter = atomList.begin(); iter != atomList.end();
       iter++) {
    double r = (**iter).GetVdwRadius();
    spReceptorGrid->SetSphere((**iter).GetCoords(), r + rVolIncr, recVal, true);
  }

  if (iTrace > 1) {
    std::cout << std::endl << "EXCLUDE RECEPTOR VOLUME" << std::endl;
    std::cout << "N(receptor)=" << spReceptorGrid->Count(recVal) << std::endl;
    std::cout << "N(excluded)=" << spReceptorGrid->Count(excVal) << std::endl;
    std::cout << "N(border)=" << spReceptorGrid->Count(borVal) << std::endl;
    std::cout << "N(unallocated)=" << spReceptorGrid->Count(0.0) << std::endl;
  }

  // Now map the solvent accessible regions with a large sphere
  // We first map the border region, which will also sweep out and exclude
  // regions of the user-specified inner region This is the first key step for
  // preventing edge effects.
  spReceptorGrid->SetAccessible(largeR, borVal, recVal, larVal, false);
  if (iTrace > 1) {
    std::cout << std::endl
              << "EXCLUDE LARGE SPHERE (Border region)" << std::endl;
    std::cout << "N(receptor)=" << spReceptorGrid->Count(recVal) << std::endl;
    std::cout << "N(large sphere)=" << spReceptorGrid->Count(larVal)
              << std::endl;
    std::cout << "N(excluded)=" << spReceptorGrid->Count(excVal) << std::endl;
    std::cout << "N(border)=" << spReceptorGrid->Count(borVal) << std::endl;
    std::cout << "N(unallocated)=" << spReceptorGrid->Count(0.0) << std::endl;
  }
  spReceptorGrid->SetAccessible(largeR, 0.0, recVal, larVal, false);
  if (iTrace > 1) {
    std::cout << std::endl
              << "EXCLUDE LARGE SPHERE (Unallocated inner region)" << std::endl;
    std::cout << "N(receptor)=" << spReceptorGrid->Count(recVal) << std::endl;
    std::cout << "N(large sphere)=" << spReceptorGrid->Count(larVal)
              << std::endl;
    std::cout << "N(excluded)=" << spReceptorGrid->Count(excVal) << std::endl;
    std::cout << "N(border)=" << spReceptorGrid->Count(borVal) << std::endl;
    std::cout << "N(unallocated)=" << spReceptorGrid->Count(0.0) << std::endl;
  }

  // Finally with a smaller radius. This is the region we want to search for
  // peaks in, so set to a positive value But first we need to replace all
  // non-zero values with the receptor volume value otherwise SetAccessible will
  // not work properly. This is the other key step for preventing edge effects
  spReceptorGrid->ReplaceValue(borVal, recVal);
  spReceptorGrid->ReplaceValue(excVal, recVal);
  spReceptorGrid->ReplaceValue(larVal, recVal);
  spReceptorGrid->SetAccessible(smallR, 0.0, recVal, cavVal, false);

  if (iTrace > 1) {
    std::cout << std::endl << "FINAL CAVITIES" << std::endl;
    std::cout << "N(receptor)=" << spReceptorGrid->Count(recVal) << std::endl;
    std::cout << "N(large sphere)=" << spReceptorGrid->Count(larVal)
              << std::endl;
    std::cout << "N(excluded)=" << spReceptorGrid->Count(excVal) << std::endl;
    std::cout << "N(border)=" << spReceptorGrid->Count(borVal) << std::endl;
    std::cout << "N(unallocated)=" << spReceptorGrid->Count(0.0) << std::endl;
    std::cout << "N(cavities)=" << spReceptorGrid->Count(cavVal) << std::endl;
    std::cout << std::endl << "Min cavity size=" << minSize << std::endl;
  }

  // Find the contiguous regions of cavity grid points
  RbtFFTPeakMap peakMap = spReceptorGrid->FindPeaks(cavVal, minSize);

  // Convert peaks to cavity format
  for (RbtFFTPeakMapConstIter pIter = peakMap.begin(); pIter != peakMap.end();
       pIter++) {
    RbtFFTPeakPtr spPeak((*pIter).second);
    RbtCoordList coordList = spReceptorGrid->GetCoordList(spPeak->points);
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
