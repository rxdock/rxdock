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

// Interface to all grid classes. Provides conversion methods between different
// addressing methods: (iXYZ),(iX,iY,iZ), and (x,y,z) coords. Underlying data
// values are managed by subclasses.

#ifndef _RBTBASEGRID_H_
#define _RBTBASEGRID_H_

#include "Config.h"
#include "Coord.h"

#include <set>

namespace rxdock {

class BaseGrid {
public:
  // Class type string
  static std::string _CT;

  ////////////////////////////////////////
  // Constructors/destructors
  // Construct a NXxNYxNZ grid running from gridMin at gridStep resolution
  BaseGrid(const Coord &gridMin, const Vector &gridStep, unsigned int NX,
           unsigned int NY, unsigned int NZ, unsigned int NPad = 0);
  // Constructor reading all params from binary stream
  BaseGrid(json j);
  BaseGrid(const BaseGrid &);            // Copy constructor
  BaseGrid &operator=(const BaseGrid &); // Copy assignment
  virtual ~BaseGrid();                   // Default destructor

  // Friend functions
  // Insertion operator (primarily for debugging)
  friend std::ostream &operator<<(std::ostream &s, const BaseGrid &grid);
  friend void to_json(json &j, const BaseGrid &grid);
  friend void from_json(const json &j, BaseGrid &grid);

  ////////////////////////////////////////
  // Virtual functions for reading/writing grid data to streams in
  // text and binary format
  // Subclasses should provide their own private OwnPrint
  // method to handle subclass data members, and override the public
  // Print method
  virtual void Print(std::ostream &ostr) const; // Text output

  ////////////////////////////////////////
  // Public methods
  ////////////////

  /////////////////////////
  // Get attribute functions
  /////////////////////////

  unsigned int GetNX() const { return m_NX; }
  unsigned int GetNY() const { return m_NY; }
  unsigned int GetNZ() const { return m_NZ; }
  unsigned int GetN() const { return m_N; }
  unsigned int GetStrideX() const { return m_SX; }
  unsigned int GetStrideY() const { return m_SY; }
  unsigned int GetStrideZ() const { return m_SZ; }
  int GetnXMin() const { return m_nXMin; }
  int GetnYMin() const { return m_nYMin; }
  int GetnZMin() const { return m_nZMin; }
  int GetnXMax() const { return m_nXMax; }
  int GetnYMax() const { return m_nYMax; }
  int GetnZMax() const { return m_nZMax; }

  const Coord &GetGridMin() const { return m_min; }
  const Coord &GetGridMax() const { return m_max; }
  const Vector &GetGridStep() const { return m_step; }
  Coord GetGridCenter() const { return (GetGridMax() + GetGridMin()) / 2.0; }
  Coord GetGridSize() const { return GetGridMax() - GetGridMin(); }
  unsigned int GetPad() const { return m_NPad; }
  const Coord &GetPadMin() const { return m_padMin; }
  const Coord &GetPadMax() const { return m_padMax; }

  /////////////////////////
  // Set attribute functions
  /////////////////////////

  // Set grid min without changing the grid step
  void SetGridMin(const Coord &gridMin);
  // Translate the grid by the given vector
  void TranslateGrid(const Vector &vec) { SetGridMin(GetGridMin() + vec); }
  // Set the center of the grid to the given coord
  void SetGridCenter(const Coord &gridCenter) {
    TranslateGrid(gridCenter - GetGridCenter());
  }
  // Set grid step
  void SetGridStep(const Vector &gridStep);
  // Set a pad region (cuboid inside the grid which defines the valid coords)
  void SetPad(unsigned int NPad = 0);

  /////////////////////////
  // Bounds checking and index conversion routines
  /////////////////////////

  // When checking coords, take account of pad region
  bool isValid(const Coord &c) const {
    return (c >= m_padMin) && (c < m_padMax);
  }
  // When checking 3-D indices, take account of pad region
  bool isValid(unsigned int iX, unsigned int iY, unsigned int iZ) const {
    return (iX >= m_NPad) && (iX < m_NX - m_NPad) && (iY >= m_NPad) &&
           (iY < m_NY - m_NPad) && (iZ >= m_NPad) && (iZ < m_NZ - m_NPad);
  }
  // When checking 1-D index, don't take account of pad region
  bool isValid(unsigned int iXYZ) const { return (iXYZ < m_N); }

  // Return iX,iY,iZ indices for given coord
  unsigned int GetIX(const Coord &c) const {
    return static_cast<unsigned int>((c.xyz(0) - m_min.xyz(0)) / m_step.xyz(0));
  }
  unsigned int GetIY(const Coord &c) const {
    return static_cast<unsigned int>((c.xyz(1) - m_min.xyz(1)) / m_step.xyz(1));
  }
  unsigned int GetIZ(const Coord &c) const {
    return static_cast<unsigned int>((c.xyz(2) - m_min.xyz(2)) / m_step.xyz(2));
  }

  // Return iX,iY,iZ indices for given coord
  unsigned int GetIX(double x) const {
    return static_cast<unsigned int>((x - m_min.xyz(0)) / m_step.xyz(0));
  }
  unsigned int GetIY(double y) const {
    return static_cast<unsigned int>((y - m_min.xyz(1)) / m_step.xyz(1));
  }
  unsigned int GetIZ(double z) const {
    return static_cast<unsigned int>((z - m_min.xyz(2)) / m_step.xyz(2));
  }

  // Return iX,iY,iZ indices for given iXYZ index
  unsigned int GetIX(unsigned int iXYZ) const { return iXYZ / m_SX; }
  unsigned int GetIY(unsigned int iXYZ) const { return (iXYZ % m_SX) / m_SY; }
  unsigned int GetIZ(unsigned int iXYZ) const { return (iXYZ % m_SY) / m_SZ; }

  // Return iXYZ index for given iX,iY,iZ indices
  unsigned int GetIXYZ(unsigned int iX, unsigned int iY,
                       unsigned int iZ) const {
    return iX * m_SX + iY * m_SY + iZ * m_SZ;
  }
  // Return iXYZ index for given coord
  unsigned int GetIXYZ(const Coord &c) const {
    return GetIXYZ(GetIX(c), GetIY(c), GetIZ(c));
  }

  // Returns real-world coordinate for given iX,iY,iZ indices
  Coord GetCoord(unsigned int iX, unsigned int iY, unsigned int iZ) const {
    return (Coord(m_nXMin, m_nYMin, m_nZMin) + Coord(iX, iY, iZ)) * m_step;
  }
  // Returns real-world X coordinate for given iX index
  double GetXCoord(unsigned int iX) const {
    return (m_nXMin + static_cast<int>(iX)) * m_step.xyz(0);
  }
  // Returns real-world Y coordinate for given iY index
  double GetYCoord(unsigned int iY) const {
    return (m_nYMin + static_cast<int>(iY)) * m_step.xyz(1);
  }
  // Returns real-world Z coordinate for given iZ index
  double GetZCoord(unsigned int iZ) const {
    return (m_nZMin + static_cast<int>(iZ)) * m_step.xyz(2);
  }

  // Returns real-world coordinate for given iXYZ index
  Coord GetCoord(unsigned int iXYZ) const {
    return GetCoord(GetIX(iXYZ), GetIY(iXYZ), GetIZ(iXYZ));
  }
  // Returns list of real-world coordinates for given set of iXYZ indices
  CoordList GetCoordList(const std::set<unsigned int> &iXYZSet) const;

  // DM 17 May 1999 - returns the set of valid grid points within a sphere of a
  // given center and radius DM 17 Jul 2000 - use std::vector<UInt> and
  // return by reference, for performance boost
  void GetSphereIndices(const Coord &c, double radius,
                        std::vector<unsigned int> &sIndices) const;

protected:
  ////////////////////////////////////////
  // Protected methods
  ///////////////////

  // Protected method for writing data members for this class to text stream
  void OwnPrint(std::ostream &ostr) const;

private:
  ////////////////////////////////////////
  // Private methods
  /////////////////
  BaseGrid(); // Disable default constructor

  // Helper function called by copy constructor and assignment operator
  void CopyGrid(const BaseGrid &);

protected:
  ////////////////////////////////////////
  // Protected data
  ////////////////

private:
  ////////////////////////////////////////
  // Private data
  //////////////
  unsigned int
      m_NX; // Number of grid points in X direction (fixed in constructor)
  unsigned int
      m_NY; // Number of grid points in Y direction (fixed in constructor)
  unsigned int
      m_NZ; // Number of grid points in Z direction (fixed in constructor)
  unsigned int m_N;    // Total number of grid points (fixed in constructor)
  unsigned int m_SX;   // Stride of X
  unsigned int m_SY;   // Stride of Y
  unsigned int m_SZ;   // Stride of Z
  Coord m_min;         // Minimum grid coords (real-world units)
  Coord m_max;         // Maximum grid coords (real-world units)
  Vector m_step;       // Grid step in X,Y,Z (real-world units)
  unsigned int m_NPad; // Defines a zero-padded border around the grid
  Coord m_padMin;      // Minimum pad coords (used to define an active region of
                       // the grid)
  Coord m_padMax;      // Maximum pad coords (used to define an active region of
                       // the grid)
  int m_nXMin;         // Minimum X grid point (in units of grid step)
  int m_nXMax;         // Maximum X grid point (in units of grid step)
  int m_nYMin;         // Minimum Y grid point (in units of grid step)
  int m_nYMax;         // Maximum Y grid point (in units of grid step)
  int m_nZMin;         // Minimum Z grid point (in units of grid step)
  int m_nZMax;         // Maximum Z grid point (in units of grid step)
};

// Useful typedefs
typedef SmartPtr<BaseGrid> BaseGridPtr; // Smart pointer

void to_json(json &j, const BaseGrid &grid);
void from_json(const json &j, BaseGrid &grid);

} // namespace rxdock

#endif //_RBTBASEGRID_H_
