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

using namespace rxdock;

// Static data members
std::string RbtRealGrid::_CT("RbtRealGrid");

////////////////////////////////////////
// Constructors/destructors
// Construct a NXxNYxNZ grid running from gridMin at gridStep resolution
RbtRealGrid::RbtRealGrid(const RbtCoord &gridMin, const RbtCoord &gridStep,
                         unsigned int NX, unsigned int NY, unsigned int NZ,
                         unsigned int NPad)
    : RbtBaseGrid(gridMin, gridStep, NX, NY, NZ, NPad), m_tol(0.001) {
  CreateArrays();
  // Initialise the grid to zero
  SetAllValues(0.0);

  _RBTOBJECTCOUNTER_CONSTR_("RbtRealGrid");
}

// Constructor reading params from binary stream
RbtRealGrid::RbtRealGrid(std::istream &istr) : RbtBaseGrid(istr) {
  // Base class constructor has already read the grid dimensions
  // etc, so all we have to do here is created the array
  // and read in the grid values
  CreateArrays();
  OwnRead(istr);
  _RBTOBJECTCOUNTER_CONSTR_("RbtRealGrid");
}

// Default destructor
RbtRealGrid::~RbtRealGrid() { _RBTOBJECTCOUNTER_DESTR_("RbtRealGrid"); }

// Copy constructor
RbtRealGrid::RbtRealGrid(const RbtRealGrid &grid) : RbtBaseGrid(grid) {
  // Base class constructor has already been called
  // so we just need to create the array and copy the array values
  CreateArrays();
  CopyGrid(grid);
  _RBTOBJECTCOUNTER_COPYCONSTR_("RbtRealGrid");
}

// Copy constructor taking a base class argument
// Sets up the grid dimensions, then creates an empty data array
RbtRealGrid::RbtRealGrid(const RbtBaseGrid &grid) : RbtBaseGrid(grid) {
  CreateArrays();
  SetAllValues(0.0);
  _RBTOBJECTCOUNTER_COPYCONSTR_("RbtRealGrid");
}

// Copy assignment
RbtRealGrid &RbtRealGrid::operator=(const RbtRealGrid &grid) {
  if (this != &grid) {
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
  Eigen::Vector3d rXYZ =
      gridStep.xyz.Constant(1.0).array() / gridStep.xyz.array();
  double rx = rXYZ(0); // reciprocal of grid step (x)
  double ry = rXYZ(1); // reciprocal of grid step (y)
  double rz = rXYZ(2); // reciprocal of grid step (z)
  // Get lower left corner grid point
  //(not necessarily the nearest grid point as returned by GetIX() etc)
  // Need to shift the int(..) argument by half a grid step
  Eigen::Vector3d diff = c.xyz - gridMin.xyz;
  unsigned int iX = static_cast<unsigned int>(rx * diff(0) - 0.5);
  unsigned int iY = static_cast<unsigned int>(ry * diff(1) - 0.5);
  unsigned int iZ = static_cast<unsigned int>(rz * diff(2) - 0.5);
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
  // bx[1] = rx * p.xyz(0);
  // bx[0] = 1.0 - bx[1];
  // by[1] = ry * p.xyz(1);
  // by[0] = 1.0 - by[1];
  // bz[1] = rz * p.xyz(2);
  // bz[0] = 1.0 - bz[1];

  double val(0.0);
  // DM 3/5/2005 - fully unroll this loop
  double bx1 = rx * p.xyz(0);
  double bx0 = 1.0 - bx1;
  double by1 = ry * p.xyz(1);
  double by0 = 1.0 - by1;
  double bz1 = rz * p.xyz(2);
  double bz0 = 1.0 - bz1;
  double bx0by0 = bx0 * by0;
  double bx0by1 = bx0 * by1;
  double bx1by0 = bx1 * by0;
  double bx1by1 = bx1 * by1;
  val += m_grid(iX, iY, iZ) * bx0by0 * bz0;
  val += m_grid(iX, iY, iZ + 1) * bx0by0 * bz1;
  val += m_grid(iX, iY + 1, iZ) * bx0by1 * bz0;
  val += m_grid(iX, iY + 1, iZ + 1) * bx0by1 * bz1;
  val += m_grid(iX + 1, iY, iZ) * bx1by0 * bz0;
  val += m_grid(iX + 1, iY, iZ + 1) * bx1by0 * bz1;
  val += m_grid(iX + 1, iY + 1, iZ) * bx1by1 * bz0;
  val += m_grid(iX + 1, iY + 1, iZ + 1) * bx1by1 * bz1;
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
void RbtRealGrid::SetAllValues(double val) { m_grid.setConstant(val); }

// Replaces all grid points between oldValMin and oldValMax with newVal
void RbtRealGrid::ReplaceValueRange(double oldValMin, double oldValMax,
                                    double newVal) {
  Eigen::Tensor<bool, 3, Eigen::RowMajor> bGrid =
      m_grid >= m_grid.constant(oldValMin) &&
      m_grid < m_grid.constant(oldValMax);
  m_grid = bGrid.select(m_grid.constant(newVal), m_grid);
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
  unsigned int iMinX = GetPad();
  unsigned int iMinY = GetPad();
  unsigned int iMinZ = GetPad();
  unsigned int iMaxX = GetNX() - GetPad();
  unsigned int iMaxY = GetNY() - GetPad();
  unsigned int iMaxZ = GetNZ() - GetPad();

  // probably could be better written
  for (unsigned int iX = iMinX; iX < iMaxX; iX++) {
    for (unsigned int iY = iMinY; iY < iMaxY; iY++) {
      for (unsigned int iZ = iMinZ; iZ < iMaxZ; iZ++) {
        // We have a match with oldVal
        if (std::fabs(m_grid(iX, iY, iZ) - oldVal) < m_tol) {
          // Check the six adjacent points for a match with adjVal
          if (((iX > iMinX) &&
               (std::fabs(m_grid(iX - 1, iY, iZ) - adjVal) < m_tol)) ||
              ((iX < iMaxX) &&
               (std::fabs(m_grid(iX + 1, iY, iZ) - adjVal) < m_tol)) ||
              ((iY > iMinY) &&
               (std::fabs(m_grid(iX, iY - 1, iZ) - adjVal) < m_tol)) ||
              ((iY < iMaxY) &&
               (std::fabs(m_grid(iX, iY + 1, iZ) - adjVal) < m_tol)) ||
              ((iZ > iMinZ) &&
               (std::fabs(m_grid(iX, iY, iZ - 1) - adjVal) < m_tol)) ||
              ((iZ < iMaxZ) &&
               (std::fabs(m_grid(iX, iY, iZ + 1) - adjVal) < m_tol)))
            m_grid(iX, iY, iZ) = newVal;
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
  unsigned int iMinX = GetPad();
  unsigned int iMinY = GetPad();
  unsigned int iMinZ = GetPad();
  unsigned int iMaxX = GetNX() - GetPad();
  unsigned int iMaxY = GetNY() - GetPad();
  unsigned int iMaxZ = GetNZ() - GetPad();

  // Work out the maximum no. of grid points in the sphere and reserve enough
  // space in the indices vector. Actually, this is a considerable overestimate
  // (no. of points in the enclosing cube)
  std::vector<unsigned int> sphereIndices;
  Eigen::Vector3d gridStepXYZ = GetGridStep().xyz;
  Eigen::Vector3d nXYZ =
      gridStepXYZ.Constant(radius).array() / gridStepXYZ.array();
  unsigned int nMax = (static_cast<unsigned int>(nXYZ(0)) + 1) *
                      (static_cast<unsigned int>(nXYZ(1)) + 1) *
                      (static_cast<unsigned int>(nXYZ(2)) + 1);
  sphereIndices.reserve(nMax);
  // possibly can be better written
  for (unsigned int iX = iMinX; iX < iMaxX; iX++) {
    for (unsigned int iY = iMinY; iY < iMaxY; iY++) {
      for (unsigned int iZ = iMinZ; iZ < iMaxZ; iZ++) {
        // We have a match with oldVal
        if (std::fabs(m_grid(iX, iY, iZ) - oldVal) < m_tol) {
          RbtCoord c = GetCoord(iX, iY, iZ);
          // Check the sphere around this grid point
          GetSphereIndices(c, radius, sphereIndices);
          if (!isValueWithinList(sphereIndices, adjVal)) {
            if (bCenterOnly)
              m_grid(iX, iY, iZ) = newVal; // Set just the center grid point
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
  Eigen::Tensor<bool, 3, Eigen::RowMajor> bGrid =
      m_grid >= m_grid.constant(valMin) && m_grid < m_grid.constant(valMax);
  Eigen::Tensor<unsigned int, 0, Eigen::RowMajor> tN =
      bGrid.cast<unsigned int>().sum();
  return tN(0);
}

// Min/max values
double RbtRealGrid::MinValue() const {
  Eigen::Tensor<float, 0, Eigen::RowMajor> tMinimum = m_grid.minimum();
  return tMinimum(0);
}

double RbtRealGrid::MaxValue() const {
  Eigen::Tensor<float, 0, Eigen::RowMajor> tMaximum = m_grid.maximum();
  return tMaximum(0);
}

// iXYZ index of grid point with minimum value
unsigned int RbtRealGrid::FindMinValue() const {
  unsigned int iMin = 0;
  const float *data = m_grid.data();
  for (unsigned int i = 0; i < GetN(); i++) {
    if (data[i] < data[iMin])
      iMin = i;
  }
  return iMin;
}

// iXYZ index of grid point with maximum value
unsigned int RbtRealGrid::FindMaxValue() const {
  unsigned int iMax = 0;
  const float *data = m_grid.data();
  for (unsigned int i = 0; i < GetN(); i++) {
    if (data[i] > data[iMax])
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
  s << std::setw(8) << GetGridStep().xyz(0) * (GetNX() - 1) << std::setw(8)
    << GetGridStep().xyz(1) * (GetNY() - 1) << std::setw(8)
    << GetGridStep().xyz(2) * (GetNZ() - 1) << std::setw(8) << 90.0
    << std::setw(8) << 90.0 << std::setw(8) << 90.0 << std::endl;
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
  for (unsigned int iZ = 0; iZ < GetNZ(); iZ++) {
    for (unsigned int iY = 0; iY < GetNY(); iY++) {
      for (unsigned int iX = 0; iX < GetNX(); iX++) {
        s << std::setw(15) << m_grid(iX, iY, iZ) << std::endl;
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
  // TODO use chip()
  for (unsigned int iX = 0; iX < GetNX(); iX++) {
    ostr << std::endl << std::endl << "Plane iX=" << iX << std::endl;
    for (unsigned int iY = 0; iY < GetNY(); iY++) {
      for (unsigned int iZ = 0; iZ < GetNZ(); iZ++) {
        float f = m_grid(iX, iY, iZ);
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
  WriteWithThrow(ostr, (const char *)&length, sizeof(length));
  WriteWithThrow(ostr, gridTitle, length);

  // Write all the data members
  WriteWithThrow(ostr, (const char *)&m_tol, sizeof(m_tol));
  const float *data = m_grid.data();
  for (unsigned int i = 0; i < GetN(); i++) {
    WriteWithThrow(ostr, (const char *)&data[i], sizeof(data[i]));
  }
}

// Protected method for reading data members for this class from binary stream
// WARNING: Assumes grid data array has already been created
// and is of the correct size
void RbtRealGrid::OwnRead(std::istream &istr) {
  // Read title
  int length;
  ReadWithThrow(istr, (char *)&length, sizeof(length));
  char *gridTitle = new char[length + 1];
  ReadWithThrow(istr, gridTitle, length);
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
  ReadWithThrow(istr, (char *)&m_tol, sizeof(m_tol));
  float *data = m_grid.data();
  for (unsigned int i = 0; i < GetN(); i++) {
    ReadWithThrow(istr, (char *)&data[i], sizeof(data[i]));
  }
}

///////////////////////////////////////////////////////////////////////////
// Private methods

// DM 17 Jul 2000 - analogous to isValueWithinSphere but iterates over arbitrary
// set of IXYZ indices. Private method as there is no error checking on iXYZ
// values out of bounds
bool RbtRealGrid::isValueWithinList(const std::vector<unsigned int> &iXYZList,
                                    double val) {
  float *data = m_grid.data();
  for (std::vector<unsigned int>::const_iterator iter = iXYZList.begin();
       iter != iXYZList.end(); iter++) {
    if (std::fabs(data[*iter] - val) < m_tol) {
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
  float *data = m_grid.data();
  for (std::vector<unsigned int>::const_iterator iter = iXYZList.begin();
       iter != iXYZList.end(); iter++) {
    if (bOverwrite || (std::fabs(data[*iter]) < m_tol)) {
      data[*iter] = val;
    }
  }
}

void RbtRealGrid::CreateArrays() {
  int nX = GetNX();
  int nY = GetNY();
  int nZ = GetNZ();
  m_grid = Eigen::Tensor<float, 3, Eigen::RowMajor>(nX, nY, nZ);
}

// Helper function called by copy constructor and assignment operator
// Need to create a new deep copy of the data array
// Gets called after array has been created, and base class copy has been done
void RbtRealGrid::CopyGrid(const RbtRealGrid &grid) {
  m_tol = grid.m_tol;
  m_grid = grid.m_grid;
}
