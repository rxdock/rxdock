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

// Subclass of RealGrid with extra methods for finding peaks
// The original FFT functionality has long since been abandoned

#ifndef _RBTFFTGRID_H_
#define _RBTFFTGRID_H_

#include "rxdock/RealGrid.h"

namespace rxdock {

// struct for holding FFT Peak info
class FFTPeak {
public:
  // Constructor
  FFTPeak() : index(0), coord(0.0, 0.0, 0.0), height(0.0), volume(0) {}

  // Data attributes - public to avoid need for accessor functions
  unsigned int index;  // Location of peak maximum - iXYZ index into grid
  Coord coord;         // Location of peak maximum - real world vector/coord
                       // corresponding to displacement
  double height;       // Maximum value of peak
  unsigned int volume; // Number of grid points in peak (=points.size())
  std::set<unsigned int> points; // Set of iXYZ peak indices
};

typedef SmartPtr<FFTPeak> FFTPeakPtr; // Smart pointer
// Multi-map of key=peak height,value=FFTPeakPtr
typedef std::multimap<double, FFTPeakPtr> FFTPeakMap;
typedef FFTPeakMap::iterator FFTPeakMapIter;
typedef FFTPeakMap::reverse_iterator FFTPeakMapRIter;
typedef FFTPeakMap::const_iterator FFTPeakMapConstIter;
typedef FFTPeakMap::const_reverse_iterator FFTPeakMapConstRIter;

class FFTGrid : public RealGrid {
public:
  // Class type string
  static std::string _CT;

  ////////////////////////////////////////
  // Constructors/destructors
  // Construct a NXxNYxNZ grid running from gridMin at gridStep resolution
  FFTGrid(const Coord &gridMin, const Coord &gridStep, unsigned int NX,
          unsigned int NY, unsigned int NZ, unsigned int NPad = 0);

  // Constructor reading all params from binary stream
  FFTGrid(json j);

  ~FFTGrid(); // Default destructor

  FFTGrid(const FFTGrid &); // Copy constructor
  // Copy constructors taking base class arguments
  FFTGrid(const RealGrid &);
  FFTGrid(const BaseGrid &);

  FFTGrid &operator=(const FFTGrid &); // Copy assignment
  // Copy assignment taking base class arguments
  FFTGrid &operator=(const RealGrid &);
  FFTGrid &operator=(const BaseGrid &);

  ////////////////////////////////////////
  // Public methods
  ////////////////

  ////////////////////////////////////////
  // Virtual functions for reading/writing grid data to streams in
  // text and binary format
  // Subclasses should provide their own private OwnPrint
  // method to handle subclass data members, and override the public
  // Print method
  virtual void Print(std::ostream &ostr) const; // Text output

  // Find the coords of all (separate) peaks above the threshold value
  // whose volumes are not less than minVol
  // Returns a map of FFTPeaks
  FFTPeakMap FindPeaks(double threshold, unsigned int minVol = 1);
  // Returns the grid point with the maximum value in FFTPeak format
  // Just a wrapper around FindMaxValue() (see below)
  FFTPeak FindMaxPeak() const;

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
  FFTGrid(); // Disable default constructor

  // Helper function called by copy constructor and assignment operator
  void CopyGrid(const FFTGrid &);

protected:
  ////////////////////////////////////////
  // Protected data
  ////////////////

private:
  ////////////////////////////////////////
  // Private data
  //////////////
};

// Useful typedefs
typedef SmartPtr<FFTGrid> FFTGridPtr; // Smart pointer

} // namespace rxdock

#endif //_RBTFFTGRID_H_
