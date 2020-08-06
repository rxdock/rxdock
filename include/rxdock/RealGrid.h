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

// Manages an array of floats representing a 3-D lattice of grid points

#ifndef _RBTREALGRID_H_
#define _RBTREALGRID_H_

#ifdef __PGI
#define EIGEN_DONT_VECTORIZE
#endif
#ifdef _WIN32
#define NOMINMAX
#endif
#include <unsupported/Eigen/CXX11/Tensor>

#include "rxdock/BaseGrid.h"

namespace rxdock {

class RealGrid : public BaseGrid {
public:
  // Class type string
  static std::string _CT;

  ////////////////////////////////////////
  // Constructors/destructors
  // Construct a NXxNYxNZ grid running from gridMin at gridStep resolution
  RBTDLL_EXPORT RealGrid(const Coord &gridMin, const Coord &gridStep,
                         unsigned int NX, unsigned int NY, unsigned int NZ,
                         unsigned int NPad = 0);

  // Constructor reading all params from binary stream
  RBTDLL_EXPORT RealGrid(json j);

  ~RealGrid(); // Default destructor

  // Copy constructor
  RealGrid(const RealGrid &);
  // Copy constructor taking a base class argument
  // Sets up the grid dimensions, then creates an empty data array
  RealGrid(const BaseGrid &);
  // Copy assignment
  RealGrid &operator=(const RealGrid &);
  // Copy assignment taking a base class argument
  // Sets up the grid dimensions, then creates an empty data array
  RealGrid &operator=(const BaseGrid &);

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
  float *GetGridData() {
    float *data = m_grid.data();
    return data;
  }

  /////////////////////////
  // Get/Set value functions
  /////////////////////////

  // Get/Set tolerance used for assessing value equality
  double GetTolerance() const { return m_tol; }
  void SetTolerance(double tol) { m_tol = tol; }

  // Get/Set single grid point value with bounds checking
  double GetValue(const Coord &c) const {
    return isValid(c) ? m_grid(GetIX(c), GetIY(c), GetIZ(c)) : 0.0;
  }

  double GetValue(unsigned int iX, unsigned int iY, unsigned int iZ) const {
    return isValid(iX, iY, iZ) ? m_grid(iX, iY, iZ) : 0.0;
  }

  double GetValue(unsigned int iXYZ) const {
    const float *data = m_grid.data();
    return isValid(iXYZ) ? data[iXYZ] : 0.0;
  }

  // DM 20 Jul 2000 - get values smoothed by trilinear interpolation
  // D. Oberlin and H.A. Scheraga, J. Comp. Chem. (1998) 19, 71.
  RBTDLL_EXPORT double GetSmoothedValue(const Coord &c) const;

  void SetValue(const Coord &c, double val) {
    if (isValid(c))
      m_grid(GetIX(c), GetIY(c), GetIZ(c)) = val;
  }

  void SetValue(unsigned int iX, unsigned int iY, unsigned int iZ, double val) {
    if (isValid(iX, iY, iZ))
      m_grid(iX, iY, iZ) = val;
  }

  void SetValue(unsigned int iXYZ, double val) {
    float *data = m_grid.data();
    if (isValid(iXYZ))
      data[iXYZ] = val;
  }

  // Set all grid points to the given value
  RBTDLL_EXPORT void SetAllValues(double val);
  // Replaces all grid points with a given value (+/- tolerance) with newVal
  void ReplaceValue(double oldVal, double newVal) {
    ReplaceValueRange(oldVal - m_tol, oldVal + m_tol, newVal);
  }
  // Replaces all grid points between oldValMin and oldValMax with newVal
  void ReplaceValueRange(double oldValMin, double oldValMax, double newVal);

  // Set all grid points within radius of coord to the given value
  // If bOverwrite is false, does not replace non-zero values
  // If bOverwrite is true, all grid points are set the new value
  void SetSphere(const Coord &c, double radius, double val,
                 bool bOverwrite = true);

  // Set all grid points with radii between rad1 and rad2 from coord to the
  // given value If bOverwrite is false, does not replace non-zero values If
  // bOverwrite is true, all grid points are set the new value
  void SetSurface(const Coord &c, double innerRad, double outerRad, double val,
                  bool bOverwrite = true);

  // Sets all grid points with value=oldValue, which are adjacent to those with
  // value=adjacentValue, to value=newValue
  //+/- tolerance is applied to oldValue and adjacentValue
  void CreateSurface(double oldVal, double adjVal, double newVal);

  // DM 16 Apr 1999 - helper function for determining solvent accessible regions
  // Returns true if any of the grid points within a sphere around the central
  // coord have the specified value
  bool isValueWithinSphere(const Coord &c, double radius, double val);
  // Sets all grid points with value=oldValue, which have no grid points with
  // value=adjacentValue within a sphere of given radius, to value=newValue
  //+/- tolerance is applied to oldValue and adjacentValue
  // If bCenterOnly is true, just the center of the sphere is set the newValue
  // If bCenterOnly is false, all grid points in the sphere are set to the
  // newValue
  void SetAccessible(double radius, double oldVal, double adjVal, double newVal,
                     bool bCenterOnly = true);

  /////////////////////////
  // Statistical functions
  /////////////////////////

  // Returns number of occurrences of a given value (+/- tolerance)
  unsigned int Count(double val) const {
    return CountRange(val - m_tol, val + m_tol);
  }
  // Returns number of occurrences of a given value range
  unsigned int CountRange(double valMin, double valMax) const;

  // Min/max values
  double MinValue() const;
  double MaxValue() const;

  unsigned int
  FindMinValue() const; // iXYZ index of grid point with minimum value
  unsigned int
  FindMaxValue() const; // iXYZ index of grid point with maximum value

  /////////////////////////
  // I/O functions
  /////////////////////////

  // Dump grid in a format readable by Insight
  RBTDLL_EXPORT void PrintInsightGrid(std::ostream &s) const;

  RBTDLL_EXPORT friend void to_json(json &j, const RealGrid &grid);
  friend void from_json(const json &j, RealGrid &grid);

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
  RealGrid(); // Disable default constructor

  // DM 17 Jul 2000 - analogous to isValueWithinSphere but iterates over
  // arbitrary set of IXYZ indices. Private method as there is no error checking
  // on iXYZ values out of bounds
  bool isValueWithinList(const std::vector<unsigned int> &iXYZList, double val);
  // Analogous to SetSphere but iterates over arbitrary set
  // of IXYZ indices. Private method as there is no error checking on iXYZ
  // values out of bounds If bOverwrite is false, does not replace non-zero
  // values If bOverwrite is true, all grid points are set the new value
  void SetValues(const std::vector<unsigned int> &iXYZList, double val,
                 bool bOverwrite = true);

  void CreateArrays();

  // Helper function called by copy constructor and assignment operator
  void CopyGrid(const RealGrid &);

protected:
  ////////////////////////////////////////
  // Protected data
  ////////////////

private:
  ////////////////////////////////////////
  // Private data
  //////////////
  Eigen::Tensor<float, 3, Eigen::RowMajor>
      m_grid;   // 3-D array, accessed as m_grid(i, j, k), indicies from 0
  double m_tol; // Tolerance for comparing grid values;
};

// Useful typedefs
typedef SmartPtr<RealGrid> RealGridPtr; // Smart pointer
typedef std::vector<RealGridPtr> RealGridList;
typedef RealGridList::iterator RealGridListIter;
typedef RealGridList::const_iterator RealGridListConstIter;

RBTDLL_EXPORT void to_json(json &j, const RealGrid &grid);
void from_json(const json &j, RealGrid &grid);

} // namespace rxdock

#endif //_RBTREALGRID_H_
