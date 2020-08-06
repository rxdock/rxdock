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

#include <algorithm> //for min, max, count
#include <cstring>
#include <iomanip>
#include <queue>

#include "rxdock/FFTGrid.h"
#include "rxdock/FileError.h"

#include <fmt/ostream.h>
#include <loguru.hpp>

using namespace rxdock;

// Static data members
std::string FFTGrid::_CT("FFTGrid");

////////////////////////////////////////
// Constructors/destructors
// Construct a NXxNYxNZ grid running from gridMin at gridStep resolution
FFTGrid::FFTGrid(const Coord &gridMin, const Coord &gridStep, unsigned int NX,
                 unsigned int NY, unsigned int NZ, unsigned int NPad)
    : RealGrid(gridMin, gridStep, NX, NY, NZ, NPad) {
  _RBTOBJECTCOUNTER_CONSTR_("FFTGrid");
}

// Constructor reading params from binary stream
FFTGrid::FFTGrid(json j) : RealGrid(j) {
  j.get_to(*this);
  _RBTOBJECTCOUNTER_CONSTR_("FFTGrid");
}

// Default destructor
FFTGrid::~FFTGrid() { _RBTOBJECTCOUNTER_DESTR_("FFTGrid"); }

// Copy constructor
FFTGrid::FFTGrid(const FFTGrid &grid) : RealGrid(grid) {
  CopyGrid(grid);
  _RBTOBJECTCOUNTER_COPYCONSTR_("FFTGrid");
}

// Copy constructors taking base class arguments
FFTGrid::FFTGrid(const RealGrid &grid) : RealGrid(grid) {
  _RBTOBJECTCOUNTER_COPYCONSTR_("FFTGrid");
}

FFTGrid::FFTGrid(const BaseGrid &grid) : RealGrid(grid) {
  _RBTOBJECTCOUNTER_COPYCONSTR_("FFTGrid");
}

// Copy assignment
FFTGrid &FFTGrid::operator=(const FFTGrid &grid) {
  if (this != &grid) {
    RealGrid::operator=(grid);
    CopyGrid(grid);
  }
  return *this;
}
// Copy assignment taking base class arguments
FFTGrid &FFTGrid::operator=(const RealGrid &grid) {
  if (this != &grid) {
    RealGrid::operator=(grid);
  }
  return *this;
}
// Copy assignment taking base class arguments
FFTGrid &FFTGrid::operator=(const BaseGrid &grid) {
  if (this != &grid) {
    RealGrid::operator=(grid);
  }
  return *this;
}

////////////////////////////////////////
// Virtual functions for reading/writing grid data to streams in
// text and binary format
// Subclasses should provide their own private OwnPrint
// method to handle subclass data members, and override the public
// Print method

// Text output
void FFTGrid::Print(std::ostream &ostr) const {
  RealGrid::Print(ostr);
  OwnPrint(ostr);
}

// Find the coords of all (separate) peaks above the threshold value
// whose volumes are not less than minVol
// Returns a map of FFTPeaks
FFTPeakMap FFTGrid::FindPeaks(double threshold, unsigned int minVol) {
  FFTPeakMap peakMap; // Initialise the return set

  // First compile a list of all points higher than the threshold
  std::set<unsigned int> stillToProcess;
  float *nrData = GetGridData();

  // DM 21 Jan 2000 - take account of tolerance when assessing threshold
  threshold -= GetTolerance();

  for (unsigned int i = 0; i < GetN(); i++) {
    if (nrData[i] > threshold)
      stillToProcess.insert(i);
  }

  LOG_F(1, "{} data points found higher than {}", stillToProcess.size(),
        threshold);

  // Repeat while we still have data points to process
  while (!stillToProcess.empty()) {

    // Start a new peak set
    std::set<unsigned int> currentPeak;
    // Queue of points to be added to current peak
    std::queue<unsigned int> toAddToPeak;

    // Seed the queue with the first unprocessed data point
    std::set<unsigned int>::iterator iter = stillToProcess.begin();
    unsigned int iXYZ0 = *iter;
    unsigned int peakPos = iXYZ0;
    float peakHeight = nrData[peakPos];
    toAddToPeak.push(iXYZ0);
    stillToProcess.erase(iter);

    LOG_F(1, "Seeding new peak at point {}", iXYZ0);

    // Repeat while the peak keeps growing
    while (!toAddToPeak.empty()) {
      // Take a new point from the front of the pending queue
      unsigned int iXYZ0 = toAddToPeak.front();
      unsigned int iX0 = GetIX(iXYZ0);
      unsigned int iY0 = GetIY(iXYZ0);
      unsigned int iZ0 = GetIZ(iXYZ0);
      toAddToPeak.pop();
      // Check if new point is higher than the current maximum
      float f = nrData[iXYZ0];
      if (f > peakHeight) {
        peakHeight = f;
        peakPos = iXYZ0;
      }
      // Add the point to the current peak
      currentPeak.insert(iXYZ0);
      LOG_F(1, "Popping point {}", iXYZ0);

      // Now check if any neighbours of the point can also be added to the queue

      // Is X+1 within range ?
      if (iX0 < GetNX()) {
        unsigned int iXYZ1 =
            iXYZ0 + GetStrideX(); // Add stride(X) to get to X+1
        // Is X+1 in the "still to process" set
        std::set<unsigned int>::iterator iter = stillToProcess.find(iXYZ1);
        if (iter != stillToProcess.end()) {
          toAddToPeak.push(iXYZ1);
          stillToProcess.erase(iter);
          LOG_F(1, "Pushing X+1 point {}", iXYZ1);
        }
      }
      // Is X-1 within range ?
      if (iX0 > 1) {
        unsigned int iXYZ1 = iXYZ0 - GetStrideX();
        // Is X-1 in the "still to process" set
        std::set<unsigned int>::iterator iter = stillToProcess.find(iXYZ1);
        if (iter != stillToProcess.end()) {
          toAddToPeak.push(iXYZ1);
          stillToProcess.erase(iter);
          LOG_F(1, "Pushing X-1 point {}", iXYZ1);
        }
      }
      // Is Y+1 within range ?
      if (iY0 < GetNY()) {
        unsigned int iXYZ1 = iXYZ0 + GetStrideY();
        // Is Y+1 in the "still to process" set
        std::set<unsigned int>::iterator iter = stillToProcess.find(iXYZ1);
        if (iter != stillToProcess.end()) {
          toAddToPeak.push(iXYZ1);
          stillToProcess.erase(iter);
          LOG_F(1, "Pushing Y+1 point {}", iXYZ1);
        }
      }
      // Is Y-1 within range ?
      if (iY0 > 1) {
        unsigned int iXYZ1 = iXYZ0 - GetStrideY();
        // Is Y-1 in the "still to process" set
        std::set<unsigned int>::iterator iter = stillToProcess.find(iXYZ1);
        if (iter != stillToProcess.end()) {
          toAddToPeak.push(iXYZ1);
          stillToProcess.erase(iter);
          LOG_F(1, "Pushing Y-1 point {}", iXYZ1);
        }
      }
      // Is Z+1 within range ?
      if (iZ0 < GetNZ()) {
        unsigned int iXYZ1 = iXYZ0 + GetStrideZ();
        // Is Z+1 in the "still to process" set
        std::set<unsigned int>::iterator iter = stillToProcess.find(iXYZ1);
        if (iter != stillToProcess.end()) {
          toAddToPeak.push(iXYZ1);
          stillToProcess.erase(iter);
          LOG_F(1, "Pushing Z+1 point {}", iXYZ1);
        }
      }
      // Is Z-1 within range ?
      if (iZ0 > 1) {
        unsigned int iXYZ1 = iXYZ0 - GetStrideZ();
        // Is Z-1 in the "still to process" set
        std::set<unsigned int>::iterator iter = stillToProcess.find(iXYZ1);
        if (iter != stillToProcess.end()) {
          toAddToPeak.push(iXYZ1);
          stillToProcess.erase(iter);
          LOG_F(1, "Pushing Z-1 point {}", iXYZ1);
        }
      }
    }

    // Store current peak in peak map if volume is not less than minVol
    if (currentPeak.size() >= minVol) {
      FFTPeakPtr spPeak(new FFTPeak);
      spPeak->index = peakPos;
      spPeak->height = peakHeight;
      spPeak->coord = GetCoord(peakPos);
      spPeak->volume = currentPeak.size();
      spPeak->points = currentPeak;
      peakMap.insert(std::pair<double, FFTPeakPtr>(peakHeight, spPeak));
      LOG_F(1, "Finished growing peak #{}, pos={}, height={:.3f}, volume={}",
            peakMap.size(), spPeak->coord, spPeak->height, spPeak->volume);
    }
  }
  return peakMap;
}

// Returns the grid point with the maximum value in FFTPeak format
// Just a wrapper around FindMaxValue() (see below)
FFTPeak FFTGrid::FindMaxPeak() const {
  FFTPeak maxPeak;
  unsigned int peakPos = FindMaxValue();
  maxPeak.index = peakPos;
  maxPeak.height = GetValue(peakPos);
  maxPeak.coord = GetCoord(peakPos);
  maxPeak.volume = 1;
  maxPeak.points.insert(peakPos);
  return maxPeak;
}

///////////////////////////////////////////////////////////////////////////
// Protected methods

// Protected method for writing data members for this class to text stream
void FFTGrid::OwnPrint(std::ostream &ostr) const {
  ostr << "Class\t" << _CT << std::endl;
}

///////////////////////////////////////////////////////////////////////////
// Private methods

// Helper function called by copy constructor and assignment operator
// No longer required
void FFTGrid::CopyGrid(const FFTGrid &grid) {}
