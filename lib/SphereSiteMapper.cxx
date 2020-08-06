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

#include "rxdock/SphereSiteMapper.h"
#include "rxdock/FFTGrid.h"

#include <fmt/ostream.h>
#include <loguru.hpp>

#include <functional>

using namespace rxdock;

// Static data member for class type
std::string SphereSiteMapper::_CT("SphereSiteMapper");
std::string SphereSiteMapper::_VOL_INCR("VOL_INCR");
std::string SphereSiteMapper::_SMALL_SPHERE("SMALL_SPHERE");
std::string SphereSiteMapper::_LARGE_SPHERE("LARGE_SPHERE");
std::string SphereSiteMapper::_GRIDSTEP("GRIDSTEP");
std::string SphereSiteMapper::_CENTER("CENTER");
std::string SphereSiteMapper::_RADIUS("RADIUS");
std::string SphereSiteMapper::_MIN_VOLUME("MIN_VOLUME");
std::string SphereSiteMapper::_MAX_CAVITIES("MAX_CAVITIES");

SphereSiteMapper::SphereSiteMapper(const std::string &strName)
    : SiteMapper(_CT, strName) {
  LOG_F(2, "SphereSiteMapper parameterised constructor");
  // Add parameters
  AddParameter(_VOL_INCR, 0.0);
  AddParameter(_SMALL_SPHERE, 1.5);
  AddParameter(_LARGE_SPHERE, 4.0);
  AddParameter(_GRIDSTEP, 0.5);
  AddParameter(_CENTER, Coord());
  AddParameter(_RADIUS, 10.0);
  AddParameter(_MIN_VOLUME, 100.0); // Min cavity volume in A^3
  AddParameter(_MAX_CAVITIES, 99);  // Max number of cavities to return
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

SphereSiteMapper::~SphereSiteMapper() {
  LOG_F(2, "SphereSiteMapper destructor");
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

CavityList SphereSiteMapper::operator()() {
  CavityList cavityList;
  ModelPtr spReceptor(GetReceptor());
  if (spReceptor.Null())
    return cavityList;

  double rVolIncr = GetParameter(_VOL_INCR);
  double smallR = GetParameter(_SMALL_SPHERE);
  double largeR = GetParameter(_LARGE_SPHERE);
  double step = GetParameter(_GRIDSTEP);
  Coord center = GetParameter(_CENTER);
  double radius = GetParameter(_RADIUS);
  double minVol = GetParameter(_MIN_VOLUME);
  unsigned int maxCavities = GetParameter(_MAX_CAVITIES);

  // Grid values
  const double recVal = -1.0;  // Receptor volume
  const double larVal = -0.75; // Accessible to large sphere
  const double excVal = -0.5;  // Excluded from calculation
  const double borVal =
      -0.25; // Border region (used for mapping large sphere only)
  const double cavVal = 1.0; // Cavities

  // Convert from min volume (in A^3) to min size (number of grid points)
  double minSize = minVol / (step * step * step);
  FFTGridPtr spReceptorGrid;
  // Only include non-H receptor atoms in the mapping
  AtomList atomList = GetAtomListWithPredicate(spReceptor->GetAtomList(),
                                               std::not1(isAtomicNo_eq(1)));
  Vector gridStep(step, step, step);

  // We extend the grid by 2*largeR on each side to eliminate edge effects in
  // the cavity mapping
  Coord minCoord = center - radius;
  Coord maxCoord = center + radius;
  double border = 2.0 * (largeR + step);
  minCoord -= border;
  maxCoord += border;
  Vector recepExtent = maxCoord - minCoord;
  Eigen::Vector3d nXYZ = recepExtent.xyz.array() / gridStep.xyz.array();
  unsigned int nX = static_cast<unsigned int>(nXYZ(0)) + 1;
  unsigned int nY = static_cast<unsigned int>(nXYZ(1)) + 1;
  unsigned int nZ = static_cast<unsigned int>(nXYZ(2)) + 1;
  spReceptorGrid = FFTGridPtr(new FFTGrid(minCoord, gridStep, nX, nY, nZ));
  center = spReceptorGrid->GetGridCenter();

  // Initialise the grid with a zero-value region in the user-specified sphere
  // surrounded by a border region of thickness = large sphere radius
  spReceptorGrid->SetAllValues(excVal);
  spReceptorGrid->SetSphere(center, radius + largeR, borVal, true);
  spReceptorGrid->SetSphere(center, radius, 0.0, true);
  LOG_F(INFO, "INITIALISATION");
  LOG_F(INFO, "Center={}", center);
  LOG_F(INFO, "Radius={}", radius);
  LOG_F(INFO, "Border={}", border);
  LOG_F(INFO, "N(excluded)={}", spReceptorGrid->Count(excVal));
  LOG_F(INFO, "N(border)={}", spReceptorGrid->Count(borVal));
  LOG_F(INFO, "N(unallocated)={}", spReceptorGrid->Count(0.0));

  // Set all vdW volume grid points to the value -1
  // Add increment to all vdw radii
  // Iterate over all receptor atoms as we want to include those whose centers
  // are outside the active site, but whose vdw volumes overlap the active site
  // region.
  for (AtomListConstIter iter = atomList.begin(); iter != atomList.end();
       iter++) {
    double r = (**iter).GetVdwRadius();
    spReceptorGrid->SetSphere((**iter).GetCoords(), r + rVolIncr, recVal, true);
  }
  LOG_F(INFO, "EXCLUDE RECEPTOR VOLUME");
  LOG_F(INFO, "N(receptor)={}", spReceptorGrid->Count(recVal));
  LOG_F(INFO, "N(excluded)={}", spReceptorGrid->Count(excVal));
  LOG_F(INFO, "N(border)={}", spReceptorGrid->Count(borVal));
  LOG_F(INFO, "N(unallocated)={}", spReceptorGrid->Count(0.0));

  // Now map the solvent accessible regions with a large sphere
  // We first map the border region, which will also sweep out and exclude
  // regions of the user-specified inner region This is the first key step for
  // preventing edge effects.
  spReceptorGrid->SetAccessible(largeR, borVal, recVal, larVal, false);
  LOG_F(INFO, "EXCLUDE LARGE SPHERE (Border region)");
  LOG_F(INFO, "N(receptor)={}", spReceptorGrid->Count(recVal));
  LOG_F(INFO, "N(large sphere)={}", spReceptorGrid->Count(larVal));
  LOG_F(INFO, "N(excluded)={}", spReceptorGrid->Count(excVal));
  LOG_F(INFO, "N(border)={}", spReceptorGrid->Count(borVal));
  LOG_F(INFO, "N(unallocated)={}", spReceptorGrid->Count(0.0));

  spReceptorGrid->SetAccessible(largeR, 0.0, recVal, larVal, false);
  LOG_F(INFO, "EXCLUDE LARGE SPHERE (Unallocated inner region)");
  LOG_F(INFO, "N(receptor)={}", spReceptorGrid->Count(recVal));
  LOG_F(INFO, "N(large sphere)={}", spReceptorGrid->Count(larVal));
  LOG_F(INFO, "N(excluded)={}", spReceptorGrid->Count(excVal));
  LOG_F(INFO, "N(border)={}", spReceptorGrid->Count(borVal));
  LOG_F(INFO, "N(unallocated)={}", spReceptorGrid->Count(0.0));

  // Finally with a smaller radius. This is the region we want to search for
  // peaks in, so set to a positive value But first we need to replace all
  // non-zero values with the receptor volume value otherwise SetAccessible will
  // not work properly. This is the other key step for preventing edge effects
  spReceptorGrid->ReplaceValue(borVal, recVal);
  spReceptorGrid->ReplaceValue(excVal, recVal);
  spReceptorGrid->ReplaceValue(larVal, recVal);
  spReceptorGrid->SetAccessible(smallR, 0.0, recVal, cavVal, false);
  LOG_F(INFO, "FINAL CAVITIES");
  LOG_F(INFO, "N(receptor)={}", spReceptorGrid->Count(recVal));
  LOG_F(INFO, "N(large sphere)={}", spReceptorGrid->Count(larVal));
  LOG_F(INFO, "N(excluded)={}", spReceptorGrid->Count(excVal));
  LOG_F(INFO, "N(border)={}", spReceptorGrid->Count(borVal));
  LOG_F(INFO, "N(unallocated)={}", spReceptorGrid->Count(0.0));
  LOG_F(INFO, "N(cavities)={}", spReceptorGrid->Count(cavVal));
  LOG_F(INFO, "Min cavity size={}", minSize);

  // Find the contiguous regions of cavity grid points
  FFTPeakMap peakMap =
      spReceptorGrid->FindPeaks(cavVal, static_cast<unsigned int>(minSize));

  // Convert peaks to cavity format
  for (FFTPeakMapConstIter pIter = peakMap.begin(); pIter != peakMap.end();
       pIter++) {
    FFTPeakPtr spPeak((*pIter).second);
    CoordList coordList = spReceptorGrid->GetCoordList(spPeak->points);
    CavityPtr spCavity = CavityPtr(new Cavity(coordList, gridStep));
    cavityList.push_back(spCavity);
  }

  // Sort cavities by volume
  std::sort(cavityList.begin(), cavityList.end(), CavityPtrCmp_Volume());

  for (CavityListConstIter cIter = cavityList.begin();
       cIter != cavityList.end(); cIter++) {
    LOG_F(1, "{}", (**cIter));
  }

  // Limit the number of cavities if necessary
  if (cavityList.size() > maxCavities) {
    LOG_F(INFO, "{} cavities identified - limit to {} largest cavities",
          cavityList.size(), maxCavities);
    cavityList.erase(cavityList.begin() + maxCavities, cavityList.end());
  }

  return cavityList;
}
