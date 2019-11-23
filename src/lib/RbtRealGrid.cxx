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

#include "RbtFileError.h"
#include "RbtRealGrid.h"

// Static data members
std::string RbtRealGrid::_CT("RbtRealGrid");

////////////////////////////////////////
// Constructors/destructors
// Construct a NXxNYxNZ grid running from gridMin at gridStep resolution
RbtRealGrid::RbtRealGrid(const RbtCoord &gridMin, const RbtCoord &gridStep,
                         unsigned int NX, unsigned int NY, unsigned int NZ,
                         unsigned int NPad)
    : RbtBaseGrid(gridMin, gridStep, NX, NY, NZ, NPad), m_grid(nullptr),
      m_data(nullptr), m_tol(0.001) {
  CreateArrays();
  // Initialise the grid to zero
  SetAllValues(0.0);

  _RBTOBJECTCOUNTER_CONSTR_("RbtRealGrid");
}

// Constructor reading params from binary stream
RbtRealGrid::RbtRealGrid(std::istream &istr)
    : RbtBaseGrid(istr), m_grid(nullptr), m_data(nullptr) {
  // Base class constructor has already read the grid dimensions
  // etc, so all we have to do here is created the array
  // and read in the grid values
  CreateArrays();
  OwnRead(istr);
  _RBTOBJECTCOUNTER_CONSTR_("RbtRealGrid");
}

// Default destructor
RbtRealGrid::~RbtRealGrid() {
  ClearArrays();
  _RBTOBJECTCOUNTER_DESTR_("RbtRealGrid");
}

// Copy constructor
RbtRealGrid::RbtRealGrid(const RbtRealGrid &grid)
    : RbtBaseGrid(grid), m_grid(nullptr), m_data(nullptr) {
  // Base class constructor has already been called
  // so we just need to create the array and copy the array values
  CreateArrays();
  CopyGrid(grid);
  _RBTOBJECTCOUNTER_COPYCONSTR_("RbtRealGrid");
}

// Copy constructor taking a base class argument
// Sets up the grid dimensions, then creates an empty data array
RbtRealGrid::RbtRealGrid(const RbtBaseGrid &grid)
    : RbtBaseGrid(grid), m_grid(nullptr), m_data(nullptr) {
  CreateArrays();
  SetAllValues(0.0);
  _RBTOBJECTCOUNTER_COPYCONSTR_("RbtRealGrid");
}

// Copy assignment
RbtRealGrid &RbtRealGrid::operator=(const RbtRealGrid &grid) {
  if (this != &grid) {
    // Clear the current grid
    ClearArrays();
    // In this case we need to explicitly call the base class operator=
    RbtBaseGrid::operator=(grid);
    CreateArrays();
    CopyGrid(grid);
  }
  return *this;
}
// Copy assignment taking a base class argument
// Sets up the grid dimensions, then creates an empty data array
RbtRealGrid &RbtRealGrid::operator=(const RbtBaseGrid &grid) {
  if (this != &grid) {
    // Clear the current grid
    ClearArrays();
    // In this case we need to explicitly call the base class operator=
    RbtBaseGrid::operator=(grid);
    CreateArrays();
    SetAllValues(0.0);
  }
  return *this;
}

////////////////////////////////////////
// Virtual functions for reading/writing grid data to streams in
// text and binary format
// Subclasses should provide their own private OwnPrint,OwnWrite, OwnRead
// methods to handle subclass data members, and override the public
// Print,Write and Read methods

// Text output
void RbtRealGrid::Print(std::ostream &ostr) const {
  RbtBaseGrid::Print(ostr);
  OwnPrint(ostr);
}

// Binary output
void RbtRealGrid::Write(std::ostream &ostr) const {
  RbtBaseGrid::Write(ostr);
  OwnWrite(ostr);
}

// Binary input
void RbtRealGrid::Read(std::istream &istr) {
  // Clear the current grid before reading the new grid dimensions
  ClearArrays();
  RbtBaseGrid::Read(istr); // Base class read
  // Now create the new grid and read in the grid values
  CreateArrays();
  OwnRead(istr);
}

////////////////////////////////////////
// Public methods
////////////////

/////////////////////////
// Get/Set value functions
/////////////////////////

// DM 20 Jul 2000 - get values smoothed by trilinear interpolation
// D. Oberlin and H.A. Scheraga, J. Comp. Chem. (1998) 19, 71.
double RbtRealGrid::GetSmoothedValue(const RbtCoord &c) const {
  const RbtCoord &gridMin = GetGridMin();
  const RbtVector &gridStep = GetGridStep();
  double rx = 1.0 / gridStep.x; // reciprocal of grid step (x)
  double ry = 1.0 / gridStep.y; // reciprocal of grid step (y)
  double rz = 1.0 / gridStep.z; // reciprocal of grid step (z)
  // Get lower left corner grid point
  //(not necessarily the nearest grid point as returned by GetIX() etc)
  // Need to shift the int(..) argument by half a grid step
  unsigned int iX = int(rx * (c.x - gridMin.x) - 0.5) + 1;
  unsigned int iY = int(ry * (c.y - gridMin.y) - 0.5) + 1;
  unsigned int iZ = int(rz * (c.z - gridMin.z) - 0.5) + 1;
#ifdef _DEBUG
  std::cout << "GetSmoothedValue" << c << "\tiX,iY,iZ=" << iX << "\t" << iY
            << "\t" << iZ << std::endl;
#endif //_DEBUG
  // Check this point (iX,iY,iZ) and (iX+1,iY+1,iZ+1) are all in bounds
  // else return the unsmoothed GetValue(c)
  if (!isValid(iX, iY, iZ) || !isValid(iX + 1, iY + 1, iZ + 1)) {
#ifdef _DEBUG
    std::cout << "Out of bounds" << std::endl;
#endif //_DEBUG
    return GetValue(c);
  }
  // p is the vector relative to the lower left corner
  RbtVector p = c - GetCoord(iX, iY, iZ);
#ifdef _DEBUG
  std::cout << "p=" << p << std::endl;
#endif //_DEBUG
  // Set up B0 and B1 for each of x,y,z axes
  // std::vector<double> bx(2);
  // std::vector<double> by(2);
  // std::vector<double> bz(2);
  // bx[1] = rx * p.x;
  // bx[0] = 1.0 - bx[1];
  // by[1] = ry * p.y;
  // by[0] = 1.0 - by[1];
  // bz[1] = rz * p.z;
  // bz[0] = 1.0 - bz[1];

  double val(0.0);
  // DM 3/5/2005 - fully unroll this loop
  double bx1 = rx * p.x;
  double bx0 = 1.0 - bx1;
  double by1 = ry * p.y;
  double by0 = 1.0 - by1;
  double bz1 = rz * p.z;
  double bz0 = 1.0 - bz1;
  double bx0by0 = bx0 * by0;
  double bx0by1 = bx0 * by1;
  double bx1by0 = bx1 * by0;
  double bx1by1 = bx1 * by1;
  val += m_grid[iX][iY][iZ] * bx0by0 * bz0;
  val += m_grid[iX][iY][iZ + 1] * bx0by0 * bz1;
  val += m_grid[iX][iY + 1][iZ] * bx0by1 * bz0;
  val += m_grid[iX][iY + 1][iZ + 1] * bx0by1 * bz1;
  val += m_grid[iX + 1][iY][iZ] * bx1by0 * bz0;
  val += m_grid[iX + 1][iY][iZ + 1] * bx1by0 * bz1;
  val += m_grid[iX + 1][iY + 1][iZ] * bx1by1 * bz0;
  val += m_grid[iX + 1][iY + 1][iZ + 1] * bx1by1 * bz1;
  // for (RbtUInt i = 0; i < 2; i++) {
  //  for (RbtUInt j = 0; j < 2; j++) {
  //    for (RbtUInt k = 0; k < 2; k++) {
  //      val += m_grid[iX+i][iY+j][iZ+k]*bx[i]*by[j]*bz[k];
  //    }
  //  }
  //}
  return val;
}

// Set all grid points to the given value
void RbtRealGrid::SetAllValues(double val) {
  for (unsigned int i = 0; i < GetN(); i++) {
    m_data[i] = val;
  }
}

// Replaces all grid points between oldValMin and oldValMax with newVal
void RbtRealGrid::ReplaceValueRange(double oldValMin, double oldValMax,
                                    double newVal) {
  for (unsigned int i = 0; i < GetN(); i++) {
    float d = m_data[i];
    if ((d >= oldValMin) && (d < oldValMax))
      m_data[i] = newVal;
  }
}

// Set all grid points within radius of coord to the given value
// If bOverwrite is false, does not replace non-zero values
// If bOverwrite is true, all grid points are set to the new value
void RbtRealGrid::SetSphere(const RbtCoord &c, double radius, double val,
                            bool bOverwrite) {
  std::vector<unsigned int> sphereIndices;
  GetSphereIndices(c, radius, sphereIndices);
  SetValues(sphereIndices, val, bOverwrite);
}

// Set all grid points with radii between rad1 and rad2 from coord to the given
// value If bOverwrite is false, does not replace non-zero values If bOverwrite
// is true, all grid points are set the new value
void RbtRealGrid::SetSurface(const RbtCoord &c, double innerRad,
                             double outerRad, double val, bool bOverwrite) {
  std::vector<unsigned int> oIndices; // List of sphere points for outer sphere
  std::vector<unsigned int> iIndices; // List of sphere points for outer sphere
  std::vector<unsigned int>
      sIndices; // List of points for surface region (points in outer
                // sphere, but not inner sphere)
  GetSphereIndices(c, outerRad, oIndices);
  GetSphereIndices(c, innerRad, iIndices);
  // std::sort(oIndices.begin(),oIndices.end());
  // std::sort(iIndices.begin(),iIndices.end());
  std::set_difference(oIndices.begin(), oIndices.end(), iIndices.begin(),
                      iIndices.end(), std::back_inserter(sIndices));
  SetValues(sIndices, val, bOverwrite);
}

// Sets all grid points with value=oldValue, which are adjacent to those with
// value=adjacentValue, to value=newValue
//+/- tolerance is applied to oldValue and adjacentValue
void RbtRealGrid::CreateSurface(double oldVal, double adjVal, double newVal) {
  // Iterate over the cuboid defined by the pad coords
  unsigned int iMinX = GetPad() + 1;
  unsigned int iMinY = GetPad() + 1;
  unsigned int iMinZ = GetPad() + 1;
  unsigned int iMaxX = GetNX() - GetPad();
  unsigned int iMaxY = GetNY() - GetPad();
  unsigned int iMaxZ = GetNZ() - GetPad();

  for (unsigned int iX = iMinX; iX <= iMaxX; iX++) {
    for (unsigned int iY = iMinY; iY <= iMaxY; iY++) {
      for (unsigned int iZ = iMinZ; iZ <= iMaxZ; iZ++) {
        // We have a match with oldVal
        if (std::fabs(m_grid[iX][iY][iZ] - oldVal) < m_tol) {
          // Check the six adjacent points for a match with adjVal
          if (((iX > iMinX) &&
               (std::fabs(m_grid[iX - 1][iY][iZ] - adjVal) < m_tol)) ||
              ((iX < iMaxX) &&
               (std::fabs(m_grid[iX + 1][iY][iZ] - adjVal) < m_tol)) ||
              ((iY > iMinY) &&
               (std::fabs(m_grid[iX][iY - 1][iZ] - adjVal) < m_tol)) ||
              ((iY < iMaxY) &&
               (std::fabs(m_grid[iX][iY + 1][iZ] - adjVal) < m_tol)) ||
              ((iZ > iMinZ) &&
               (std::fabs(m_grid[iX][iY][iZ - 1] - adjVal) < m_tol)) ||
              ((iZ < iMaxZ) &&
               (std::fabs(m_grid[iX][iY][iZ + 1] - adjVal) < m_tol)))
            m_grid[iX][iY][iZ] = newVal;
        }
      }
    }
  }
}

// DM 16 Apr 1999 - helper function for determining solvent accessible regions
// Returns true if any of the grid points within a sphere around the central
// coord have the specified value
bool RbtRealGrid::isValueWithinSphere(const RbtCoord &c, double radius,
                                      double val) {
  std::vector<unsigned int> sphereIndices;
  GetSphereIndices(c, radius, sphereIndices);
  return isValueWithinList(sphereIndices, val);
}

// Sets all grid points with value=oldValue, which have no grid points with
// value=adjacentValue within a sphere of given radius, to value=newValue
//+/- tolerance is applied to oldValue and adjacentValue
void RbtRealGrid::SetAccessible(double radius, double oldVal, double adjVal,
                                double newVal, bool bCenterOnly) {
  // Iterate over the cuboid defined by the pad coords
  unsigned int iMinX = GetPad() + 1;
  unsigned int iMinY = GetPad() + 1;
  unsigned int iMinZ = GetPad() + 1;
  unsigned int iMaxX = GetNX() - GetPad();
  unsigned int iMaxY = GetNY() - GetPad();
  unsigned int iMaxZ = GetNZ() - GetPad();

  // Work out the maximum no. of grid points in the sphere and reserve enough
  // space in the indices vector. Actually, this is a considerable overestimate
  // (no. of points in the enclosing cube)
  std::vector<unsigned int> sphereIndices;
  unsigned int nMax = (int(radius / GetGridStep().x) + 1) *
                      (int(radius / GetGridStep().y) + 1) *
                      (int(radius / GetGridStep().z) + 1);
  sphereIndices.reserve(nMax);
  for (unsigned int iX = iMinX; iX <= iMaxX; iX++) {
    for (unsigned int iY = iMinY; iY <= iMaxY; iY++) {
      for (unsigned int iZ = iMinZ; iZ <= iMaxZ; iZ++) {
        // We have a match with oldVal
        if (std::fabs(m_grid[iX][iY][iZ] - oldVal) < m_tol) {
          RbtCoord c = GetCoord(iX, iY, iZ);
          // Check the sphere around this grid point
          GetSphereIndices(c, radius, sphereIndices);
          if (!isValueWithinList(sphereIndices, adjVal)) {
            if (bCenterOnly)
              m_grid[iX][iY][iZ] = newVal; // Set just the center grid point
            else
              // SetValues(sphereIndices,newVal,false);//Set all grid points in
              // the sphere
              SetValues(sphereIndices, newVal,
                        true); // DM 3 April 2002 - now overwrite ALL gridpoints
                               // with newVal
          }
        }
      }
    }
  }
}

/////////////////////////
// Statistical functions
/////////////////////////

// Returns number of occurrences of a given value range
unsigned int RbtRealGrid::CountRange(double valMin, double valMax) const {
  unsigned int n(0);
  for (unsigned int i = 0; i < GetN(); i++) {
    float d = m_data[i];
    if ((d >= valMin) && (d < valMax))
      n++;
  }
  return n;
}

// Min/max values
double RbtRealGrid::MinValue() const {
  float fMin = m_data[0];
  for (unsigned int i = 0; i < GetN(); i++) {
    if (m_data[i] < fMin)
      fMin = m_data[i];
  }
  return fMin;
}

double RbtRealGrid::MaxValue() const {
  float fMax = m_data[0];
  for (unsigned int i = 0; i < GetN(); i++) {
    if (m_data[i] > fMax)
      fMax = m_data[i];
  }
  return fMax;
}

// iXYZ index of grid point with minimum value
unsigned int RbtRealGrid::FindMinValue() const {
  unsigned int iMin = 0;
  for (unsigned int i = 0; i < GetN(); i++) {
    if (m_data[i] < m_data[iMin])
      iMin = i;
  }
  return iMin;
}

// iXYZ index of grid point with maximum value
unsigned int RbtRealGrid::FindMaxValue() const {
  unsigned int iMax = 0;
  for (unsigned int i = 0; i < GetN(); i++) {
    if (m_data[i] > m_data[iMax])
      iMax = i;
  }
  return iMax;
}

/////////////////////////
// Output functions
/////////////////////////

// Dump grid in a format readable by Insight
void RbtRealGrid::PrintInsightGrid(std::ostream &s) const {
  s << "RBT FFT GRID" << std::endl;
  s << "(1F15.10)" << std::endl;
  s.precision(3);
  s.setf(std::ios_base::fixed, std::ios_base::floatfield);
  s.setf(std::ios_base::right, std::ios_base::adjustfield);
  s << std::setw(8) << GetGridStep().x * (GetNX() - 1) << std::setw(8)
    << GetGridStep().y * (GetNY() - 1) << std::setw(8)
    << GetGridStep().z * (GetNZ() - 1) << std::setw(8) << 90.0 << std::setw(8)
    << 90.0 << std::setw(8) << 90.0 << std::endl;
  // s.setf(0,ios_base::floatfield);
  s << std::setw(5) << GetNX() - 1 << std::setw(5) << GetNY() - 1
    << std::setw(5) << GetNZ() - 1 << std::endl;
  s << std::setw(5) << 1 << std::setw(5) << GetnXMin() << std::setw(5)
    << GetnXMax() << std::setw(5) << GetnYMin() << std::setw(5) << GetnYMax()
    << std::setw(5) << GetnZMin() << std::setw(5) << GetnZMax() << std::endl;
  s.precision(10);
  // s.setf(std::ios_base::fixed,ios_base::floatfield);
  // s.setf(std::ios_base::right,ios_base::adjustfield);
  // Insight expects data in [1][1][1],[2][1][1]..[NX][1][1] order
  for (unsigned int iZ = 1; iZ <= GetNZ(); iZ++) {
    for (unsigned int iY = 1; iY <= GetNY(); iY++) {
      for (unsigned int iX = 1; iX <= GetNX(); iX++) {
        s << std::setw(15) << m_grid[iX][iY][iZ] << std::endl;
      }
    }
  }
}

///////////////////////////////////////////////////////////////////////////
// Protected methods

// Protected method for writing data members for this class to text stream
void RbtRealGrid::OwnPrint(std::ostream &ostr) const {
  ostr << "Class\t" << _CT << std::endl;
  // Iterate over all grid points
  float tol = GetTolerance();
  for (unsigned int iX = 1; iX <= GetNX(); iX++) {
    ostr << std::endl << std::endl << "Plane iX=" << iX << std::endl;
    for (unsigned int iY = 1; iY <= GetNY(); iY++) {
      for (unsigned int iZ = 1; iZ <= GetNZ(); iZ++) {
        float f = m_grid[iX][iY][iZ];
        ostr << ((f < -tol) ? '-' : (f > tol) ? '+' : '.');
      }
      ostr << std::endl;
    }
  }
}

// Protected method for writing data members for this class to binary stream
//(Serialisation)
void RbtRealGrid::OwnWrite(std::ostream &ostr) const {
  // Write the class name as a title so we can check the authenticity of streams
  // on read
  const char *const gridTitle = _CT.c_str();
  int length = strlen(gridTitle);
  Rbt::WriteWithThrow(ostr, (const char *)&length, sizeof(length));
  Rbt::WriteWithThrow(ostr, gridTitle, length);

  // Write all the data members
  Rbt::WriteWithThrow(ostr, (const char *)&m_tol, sizeof(m_tol));
  for (unsigned int i = 0; i < GetN(); i++) {
    Rbt::WriteWithThrow(ostr, (const char *)&m_data[i], sizeof(m_data[i]));
  }
}

// Protected method for reading data members for this class from binary stream
// WARNING: Assumes grid data array has already been created
// and is of the correct size
void RbtRealGrid::OwnRead(std::istream &istr) {
  // Read title
  int length;
  Rbt::ReadWithThrow(istr, (char *)&length, sizeof(length));
  char *gridTitle = new char[length + 1];
  Rbt::ReadWithThrow(istr, gridTitle, length);
  // Add null character to end of string
  gridTitle[length] = '\0';
  // Compare title with class name
  bool match = (_CT == gridTitle);
  delete[] gridTitle;
  if (!match) {
    throw RbtFileParseError(_WHERE_,
                            "Invalid title string in " + _CT + "::Read()");
  }

  // Read all the data members
  Rbt::ReadWithThrow(istr, (char *)&m_tol, sizeof(m_tol));
  for (unsigned int i = 0; i < GetN(); i++) {
    Rbt::ReadWithThrow(istr, (char *)&m_data[i], sizeof(m_data[i]));
  }
}

///////////////////////////////////////////////////////////////////////////
// Private methods

// DM 17 Jul 2000 - analogous to isValueWithinSphere but iterates over arbitrary
// set of IXYZ indices. Private method as there is no error checking on iXYZ
// values out of bounds
bool RbtRealGrid::isValueWithinList(const std::vector<unsigned int> &iXYZList,
                                    double val) {
  for (std::vector<unsigned int>::const_iterator iter = iXYZList.begin();
       iter != iXYZList.end(); iter++) {
    if (std::fabs(m_data[*iter] - val) < m_tol) {
      return true;
    }
  }
  return false; // Got all the way without finding a match
}

// Analogous to SetSphere but iterates over arbitrary set
// of IXYZ indices. Private method as there is no error checking on iXYZ values
// out of bounds If bOverwrite is false, does not replace non-zero values If
// bOverwrite is true, all grid points are set the new value
void RbtRealGrid::SetValues(const std::vector<unsigned int> &iXYZList,
                            double val, bool bOverwrite) {
  for (std::vector<unsigned int>::const_iterator iter = iXYZList.begin();
       iter != iXYZList.end(); iter++) {
    if (bOverwrite || (std::fabs(m_data[*iter]) < m_tol)) {
      m_data[*iter] = val;
    }
  }
}

void RbtRealGrid::CreateArrays() {
  if (m_grid != nullptr) { // Clear existing grid
    ClearArrays();
  }
  int nX = GetNX();
  int nY = GetNY();
  int nZ = GetNZ();
  m_grid = new float **[nX + 1];
  m_grid[1] = new float *[nX * nY + 1];
  m_grid[1][1] = new float[nX * nY * nZ + 1];

  for (int iY = 2; iY <= nY; iY++) {
    m_grid[1][iY] = m_grid[1][iY - 1] + nZ;
  }
  for (int iX = 2; iX <= nX; iX++) {
    m_grid[iX] = m_grid[iX - 1] + nY;
    m_grid[iX][1] = m_grid[iX - 1][1] + nY * nZ;
    for (int iY = 2; iY <= nY; iY++) {
      m_grid[iX][iY] = m_grid[iX][iY - 1] + nZ;
    }
  }
  m_data = &m_grid[1][1][1];
}

void RbtRealGrid::ClearArrays() {
  delete[] m_grid[1][1];
  delete[] m_grid[1];
  delete[] m_grid;
  m_grid = nullptr;
  m_data = nullptr;
}

// Helper function called by copy constructor and assignment operator
// Need to create a new deep copy of the data array
// Gets called after array has been created, and base class copy has been done
void RbtRealGrid::CopyGrid(const RbtRealGrid &grid) {
  m_tol = grid.m_tol;
  float *gridData = grid.GetGridData();
  for (unsigned int i = 0; i < GetN(); i++) {
    m_data[i] = gridData[i];
  }
}
