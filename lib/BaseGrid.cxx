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
#include <cstring>   //for strlen
#include <iomanip>

#include "rxdock/BaseGrid.h"
#include "rxdock/FileError.h"

#include <loguru.hpp>

using namespace rxdock;

// Static data members
const std::string BaseGrid::_CT = "BaseGrid";

////////////////////////////////////////
// Constructors/destructors
// Construct a NXxNYxNZ grid running from gridMin at gridStep resolution
BaseGrid::BaseGrid(const Coord &gridMin, const Coord &gridStep, unsigned int NX,
                   unsigned int NY, unsigned int NZ, unsigned int NPad)
    : m_NX(NX), m_NY(NY), m_NZ(NZ), m_N(NX * NY * NZ), m_SX(NY * NZ), m_SY(NZ),
      m_SZ(1), m_step(gridStep), m_NPad(NPad) {
  // Set the logical (real-world) dimensions of the grid
  SetGridMin(gridMin);
  _RBTOBJECTCOUNTER_CONSTR_("BaseGrid");
}

// Constructor reading params from binary stream
BaseGrid::BaseGrid(json j) {
  j.get_to(*this);
  _RBTOBJECTCOUNTER_CONSTR_("BaseGrid");
}

// Default destructor
BaseGrid::~BaseGrid() { _RBTOBJECTCOUNTER_DESTR_("BaseGrid"); }

// Copy constructor
BaseGrid::BaseGrid(const BaseGrid &grid) {
  CopyGrid(grid);
  _RBTOBJECTCOUNTER_COPYCONSTR_("BaseGrid");
}

// Copy assignment
BaseGrid &BaseGrid::operator=(const BaseGrid &grid) {
  if (this != &grid) {
    CopyGrid(grid);
  }
  return *this;
}

// Insertion operator (primarily for debugging)
std::ostream &operator<<(std::ostream &s, const BaseGrid &grid) {
  grid.Print(s);
  return s;
}

////////////////////////////////////////
// Virtual functions for reading/writing grid data to streams in
// text and binary format
// Subclasses should provide their own private OwnPrint
// method to handle subclass data members, and override the public
// Print method

// These implementations are not strictly necessary for the base class
// as base class objects cannot be instantiated

// Text output
void BaseGrid::Print(std::ostream &ostr) const { OwnPrint(ostr); }

////////////////////////////////////////
// Public methods
////////////////

/////////////////////////
// Set attribute functions
/////////////////////////

// Set grid min without changing the grid step
void BaseGrid::SetGridMin(const Coord &gridMin) {
  // Determine the integral min and max grid point coords (in multiples of grid
  // step from origin)
  Eigen::Vector3d nXYZMin = gridMin.xyz.array() / m_step.xyz.array();
  m_nXMin = static_cast<int>(std::floor(nXYZMin(0) + 0.5));
  m_nYMin = static_cast<int>(std::floor(nXYZMin(1) + 0.5));
  m_nZMin = static_cast<int>(std::floor(nXYZMin(2) + 0.5));

  // Max integral coords are just min+N-1
  m_nXMax = m_nXMin + static_cast<int>(m_NX) - 1;
  m_nYMax = m_nYMin + static_cast<int>(m_NY) - 1;
  m_nZMax = m_nZMin + static_cast<int>(m_NZ) - 1;

  // Now determine the min and max real-world coords which are covered by the
  // grid Include a border of +/-0.5*gridStep so that grid points lie in the
  // center of each grid interval
  m_min = (Coord(m_nXMin, m_nYMin, m_nZMin) - 0.5) * m_step;
  m_max = (Coord(m_nXMax, m_nYMax, m_nZMax) + 0.5) * m_step;

  // Recalc the pad coords as we've moved the grid
  SetPad(m_NPad);
}

// Set grid step
// The integral min and max grid point coords (in multiples of grid step from
// origin) are left unchanged Only the real-world min and max coords are updated
void BaseGrid::SetGridStep(const Vector &gridStep) {
  m_step = gridStep;
  // Now determine the min and max real-world coords which are covered by the
  // grid Include a border of +/-0.5*gridStep so that grid points lie in the
  // center of each grid interval
  m_min = (Coord(m_nXMin, m_nYMin, m_nZMin) - 0.5) * m_step;
  m_max = (Coord(m_nXMax, m_nYMax, m_nZMax) + 0.5) * m_step;

  // Recalc the pad coords as we've moved the grid
  SetPad(m_NPad);
}

// Set a pad region (cuboid inside the grid which defines the valid coords)
// Coords outside the cuboid are deemed invalid
void BaseGrid::SetPad(unsigned int NPad) {
  m_NPad = NPad;
  m_padMin = m_min + m_step * m_NPad;
  m_padMax = m_max - m_step * m_NPad;
}

// Returns list of real-world coordinates for given set of iXYZ indices
CoordList BaseGrid::GetCoordList(const std::set<unsigned int> &iXYZSet) const {
  CoordList coordList;
  coordList.reserve(iXYZSet.size());
  for (std::set<unsigned int>::const_iterator iter = iXYZSet.begin();
       iter != iXYZSet.end(); iter++)
    coordList.push_back(GetCoord(*iter));
  return coordList;
}

// DM 17 May 1999 - returns the set of valid grid points within a sphere of a
// given center and radius DM 17 Jul 2000 - use std::vector<UInt> and return
// by reference, for performance boost
void BaseGrid::GetSphereIndices(const Coord &c, double radius,
                                std::vector<unsigned int> &sIndices) const {
  sIndices.clear(); // Clear the return list

  // Iterate over the cube defined by the coord and the radius
  // without exceeding the limits of the pad region
  Coord cubeMin = Max(c - radius, m_padMin);
  Coord cubeMax = Min(c + radius, m_padMax);

  // Convert to array indices
  // DM 16 Apr 1999 - check again for indices out of range
  // If cubeMax == m_padMax, iMaxX,Y,Z would be one too high
  unsigned int iMinX = std::max(GetIX(cubeMin), m_NPad);
  unsigned int iMinY = std::max(GetIY(cubeMin), m_NPad);
  unsigned int iMinZ = std::max(GetIZ(cubeMin), m_NPad);
  unsigned int iMaxX = std::min(GetIX(cubeMax) + 1, m_NX - m_NPad);
  unsigned int iMaxY = std::min(GetIY(cubeMax) + 1, m_NY - m_NPad);
  unsigned int iMaxZ = std::min(GetIZ(cubeMax) + 1, m_NZ - m_NPad);
  LOG_F(1, "iMin: {} {} {}; iMax: {} {} {}", iMinX, iMinY, iMinZ, iMaxX, iMaxY,
        iMaxZ);

  // Determine if points are inside sphere by checking x^2 + y^2 + z^2 <= rad^2
  double rad2 = radius * radius;

  double rX; // X-distance to sphere center
  unsigned int iX;
  for (iX = iMinX, rX = GetXCoord(iMinX) - c.xyz(0); iX < iMaxX;
       iX++, rX += m_step.xyz(0)) {
    double rX2 = rX * rX;
    // Only bother with the Y loop if we haven't already exceeded rad2
    if (rX2 <= rad2) {
      double rY; // Y-distance to sphere center
      unsigned int iY;
      for (iY = iMinY, rY = GetYCoord(iMinY) - c.xyz(1); iY < iMaxY;
           iY++, rY += m_step.xyz(1)) {
        double rXY2 = rX2 + rY * rY;
        // Only bother with the Z loop if we haven't already exceeded rad2
        if (rXY2 <= rad2) {
          double rZ;
          unsigned int iZ;
          unsigned int iXYZ;
          for (iZ = iMinZ, iXYZ = GetIXYZ(iX, iY, iMinZ),
              rZ = GetZCoord(iMinZ) - c.xyz(2);
               iZ < iMaxZ; iZ++, iXYZ++, rZ += m_step.xyz(2)) {
            double rXYZ2 = rXY2 + rZ * rZ;
            if (rXYZ2 <= rad2) {
              sIndices.push_back(iXYZ);
            }
          }
        }
      }
    }
  }
}

///////////////////////////////////////////////////////////////////////////
// Protected methods

// Protected method for writing data members for this class to text stream
void BaseGrid::OwnPrint(std::ostream &ostr) const {
  ostr << "Class\t" << _CT << std::endl;
  ostr << "m_min\t" << m_min << std::endl;
  ostr << "m_max\t" << m_max << std::endl;
  ostr << "m_step\t" << m_step << std::endl;
  ostr << "m_padMin\t" << m_padMin << std::endl;
  ostr << "m_padMax\t" << m_padMax << std::endl;
  ostr << "m_NX,m_NY,m_NZ\t" << m_NX << "," << m_NY << "," << m_NZ << std::endl;
  ostr << "m_N\t" << m_N << std::endl;
  ostr << "m_SX,m_SY,m_SZ\t" << m_SX << "," << m_SY << "," << m_SZ << std::endl;
  ostr << "m_NPad\t" << m_NPad << std::endl;
  ostr << "m_nXMin,m_nXMax\t" << m_nXMin << "," << m_nXMax << std::endl;
  ostr << "m_nYMin,m_nYMax\t" << m_nYMin << "," << m_nYMax << std::endl;
  ostr << "m_nZMin,m_nZMax\t" << m_nZMin << "," << m_nZMax << std::endl;
}

///////////////////////////////////////////////////////////////////////////
// Private methods

// Helper function called by copy constructor and assignment operator
// Base class doesn't store any pointers so this is a straightforward
// memberwise copy
void BaseGrid::CopyGrid(const BaseGrid &grid) {
  m_min = grid.m_min;
  m_max = grid.m_max;
  m_step = grid.m_step;
  m_padMin = grid.m_padMin;
  m_padMax = grid.m_padMax;
  m_NX = grid.m_NX;
  m_NY = grid.m_NY;
  m_NZ = grid.m_NZ;
  m_N = grid.m_N;
  m_SX = grid.m_SX;
  m_SY = grid.m_SY;
  m_SZ = grid.m_SZ;
  m_NPad = grid.m_NPad;
  m_nXMin = grid.m_nXMin;
  m_nXMax = grid.m_nXMax;
  m_nYMin = grid.m_nYMin;
  m_nYMax = grid.m_nYMax;
  m_nZMin = grid.m_nZMin;
  m_nZMax = grid.m_nZMax;
}

void rxdock::to_json(json &j, const BaseGrid &grid) {
  j = json{{"min", grid.m_min},
           {"max", grid.m_max},
           {"step", grid.m_step},
           {"pad-min", grid.m_padMin},
           {"pad-max", grid.m_padMax},
           {"nxyz", {grid.m_NX, grid.m_NY, grid.m_NZ}},
           {"n", grid.m_N},
           {"sxyz", {grid.m_SX, grid.m_SY, grid.m_SZ}},
           {"npad", grid.m_NPad},
           {"nxyz-min", {grid.m_nXMin, grid.m_nYMin, grid.m_nZMin}},
           {"nxyz-max", {grid.m_nXMax, grid.m_nYMax, grid.m_nZMax}}};
}

void rxdock::from_json(const json &j, BaseGrid &grid) {
  j.at("min").get_to(grid.m_min);
  j.at("max").get_to(grid.m_max);
  j.at("step").get_to(grid.m_step);
  j.at("pad-min").get_to(grid.m_padMin);
  j.at("pad-max").get_to(grid.m_padMax);
  j.at("nxyz").at(0).get_to(grid.m_NX);
  j.at("nxyz").at(1).get_to(grid.m_NY);
  j.at("nxyz").at(2).get_to(grid.m_NZ);
  j.at("n").get_to(grid.m_N);
  j.at("sxyz").at(0).get_to(grid.m_SX);
  j.at("sxyz").at(1).get_to(grid.m_SY);
  j.at("sxyz").at(2).get_to(grid.m_SZ);
  j.at("npad").get_to(grid.m_NPad);
  j.at("nxyz-min").at(0).get_to(grid.m_nXMin);
  j.at("nxyz-min").at(1).get_to(grid.m_nYMin);
  j.at("nxyz-min").at(2).get_to(grid.m_nZMin);
  j.at("nxyz-min").at(0).get_to(grid.m_nXMax);
  j.at("nxyz-min").at(1).get_to(grid.m_nYMax);
  j.at("nxyz-min").at(2).get_to(grid.m_nZMax);
}
