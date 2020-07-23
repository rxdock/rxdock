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

// Simple class for handling cartesian coords (points) and 3-d vectors
// The same class is used to represent both points and vectors, although the
// RbtVector typedef can be used in argument lists to distinguish the two cases

#ifndef _RBTCOORD_H_
#define _RBTCOORD_H_

#include <algorithm> //for min,max
#include <cmath>     //for sqrt
#include <iterator>  //for back_inserter
#include <numeric>   //for accumulate

#ifdef __PGI
#define EIGEN_DONT_VECTORIZE
#endif
#include <Eigen/Core>
#include <Eigen/Geometry>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "RbtConfig.h"

extern std::istream &eatSeps(std::istream &is);

namespace rxdock {

class RbtCoord {
  ///////////////////////////////////////////////
  // Data members
  ///////////////
  //(leave public so we don't need
  // to write accessor functions)
public:
  Eigen::Vector3d xyz;

  ///////////////////////////////////////////////
  // Constructors/destructors:
  ///////////////////////////
public:
  // Default constructor (initialise to zero)
  inline RbtCoord() : xyz(0.0, 0.0, 0.0) {}

  // Constructor with initial values
  inline RbtCoord(double x1, double y1, double z1) : xyz(x1, y1, z1) {}

  // Constructor using Eigen
  inline RbtCoord(Eigen::Vector3d xyz1) : xyz(xyz1) {}

  // Destructor
  virtual ~RbtCoord() {}

  // Copy constructor
  inline RbtCoord(const RbtCoord &coord) { xyz = coord.xyz; }

  // DM 19 Jul 2000 - Read,Write methods to read/write coords to binary streams
  inline std::ostream &Write(std::ostream &ostr) const {
    ostr.write((const char *)&(xyz(0)), sizeof(xyz(0)));
    ostr.write((const char *)&(xyz(1)), sizeof(xyz(1)));
    ostr.write((const char *)&(xyz(2)), sizeof(xyz(2)));
    return ostr;
  }
  inline std::istream &Read(std::istream &istr) {
    istr.read((char *)&(xyz(0)), sizeof(xyz(0)));
    istr.read((char *)&(xyz(1)), sizeof(xyz(1)));
    istr.read((char *)&(xyz(2)), sizeof(xyz(2)));
    return istr;
  }

  ///////////////////////////////////////////////
  // Operator functions:
  /////////////////////

  // Copy assignment (*this = coord)
  inline RbtCoord &operator=(const RbtCoord &coord) {
    if (this != &coord) { // beware of self-assignment
      xyz = coord.xyz;
    }
    return *this;
  }

  // Copy assignment (*this = const)
  inline RbtCoord &operator=(const double &d) {
    xyz.setConstant(d);
    return *this;
  }

  //*this += coord
  inline void operator+=(const RbtCoord &coord) { xyz += coord.xyz; }

  //*this += const
  inline void operator+=(const double &d) { xyz += xyz.Constant(d); }

  //*this -= coord
  inline void operator-=(const RbtCoord &coord) { xyz -= coord.xyz; }

  //*this -= const
  inline void operator-=(const double &d) { xyz -= xyz.Constant(d); }

  //*this *= const
  inline void operator*=(const double &d) { xyz *= d; }

  //*this /= const
  inline void operator/=(const double &d) { xyz /= d; }

  ///////////////////////////////////////////////
  // Friend functions:
  ///////////////////

  // Insertion (output) operator
  friend std::ostream &operator<<(std::ostream &s, const RbtCoord &coord) {
    return s << "(" << coord.xyz(0) << "," << coord.xyz(1) << ","
             << coord.xyz(2) << ")";
  }

  // Input operator
  // DM 15 Apr 1999
  // Valid formats (actually the comma can be any single non-whitespace
  // character):
  // x,y,z
  // (x,y,z)
  // BGD 06 Nov 2002
  // Trying to make it more versatil. The separator now is any
  // amount of white space or commas.
  // x   y,z
  // (x   y, z  )
  // uses the auxiliary function eatSeps

  friend std::istream &operator>>(std::istream &s, RbtCoord &coord) {
    double x = 0, y = 0, z = 0;
    char c = 0;
    s >> c;
    if (c == '(') {
      s >> x;
      eatSeps(s);
      s >> y;
      eatSeps(s);
      s >> z;
      eatSeps(s);
      s >> c;
      if (c != ')')
        s.clear(std::ios_base::badbit); // Missing closing bracket
    } else {
      s.putback(c);
      s >> x;
      eatSeps(s);
      s >> y;
      eatSeps(s);
      s >> z;
    }
    if (s)
      coord = RbtCoord(x, y, z);
    return s;
  }

  // Equality
  inline friend bool operator==(const RbtCoord &coord1,
                                const RbtCoord &coord2) {
    return coord1.xyz == coord2.xyz;
  }

  // Non-equality
  inline friend bool operator!=(const RbtCoord &coord1,
                                const RbtCoord &coord2) {
    return coord1.xyz != coord2.xyz;
  }

  // Greater than
  inline friend bool operator>(const RbtCoord &coord1, const RbtCoord &coord2) {
    return (coord1.xyz.array() > coord2.xyz.array()).all();
  }

  // Greater than or equal
  inline friend bool operator>=(const RbtCoord &coord1,
                                const RbtCoord &coord2) {
    return (coord1.xyz.array() >= coord2.xyz.array()).all();
  }

  // Less than
  inline friend bool operator<(const RbtCoord &coord1, const RbtCoord &coord2) {
    return (coord1.xyz.array() < coord2.xyz.array()).all();
  }

  // Less than or equal
  inline friend bool operator<=(const RbtCoord &coord1,
                                const RbtCoord &coord2) {
    return (coord1.xyz.array() <= coord2.xyz.array()).all();
  }

  // Addition (coord + coord)
  inline friend RbtCoord operator+(const RbtCoord &coord1,
                                   const RbtCoord &coord2) {
    return RbtCoord(coord1.xyz + coord2.xyz);
  }

  // Addition (coord + const)
  inline friend RbtCoord operator+(const RbtCoord &coord1, double d) {
    return RbtCoord(coord1.xyz + coord1.xyz.Constant(d));
  }

  // Addition (const + coord)
  inline friend RbtCoord operator+(double d, const RbtCoord &coord1) {
    return RbtCoord(coord1.xyz.Constant(d) + coord1.xyz);
  }

  // Subtraction (coord - coord)
  inline friend RbtCoord operator-(const RbtCoord &coord1,
                                   const RbtCoord &coord2) {
    return RbtCoord(coord1.xyz - coord2.xyz);
  }

  // Subtraction (coord - const)
  inline friend RbtCoord operator-(const RbtCoord &coord1, double d) {
    return RbtCoord(coord1.xyz - coord1.xyz.Constant(d));
  }

  // Subtraction (const - coord)
  inline friend RbtCoord operator-(double d, const RbtCoord &coord1) {
    return RbtCoord(coord1.xyz.Constant(d) - coord1.xyz);
  }

  // Negation
  inline friend RbtCoord operator-(const RbtCoord &coord) {
    return RbtCoord(-coord.xyz);
  }

  // Dot product
  // DM 7 Jan 1999 - operator* replaced by Dot(coord) member function
  // and by Dot(coord1,coord2) non-member function
  // inline friend RbtDouble operator*(const RbtCoord& coord1, const RbtCoord&
  // coord2) {
  //  return coord1.xyz.dot(coord.xyz);
  //}

  // Scalar product (coord * const)
  inline friend RbtCoord operator*(const RbtCoord &coord, const double &d) {
    return RbtCoord(coord.xyz * d);
  }

  // Scalar product (const * coord)
  inline friend RbtCoord operator*(const double &d, const RbtCoord &coord) {
    return RbtCoord(d * coord.xyz);
  }

  // Scalar product (coord * coord : component-wise multiplication)
  inline friend RbtCoord operator*(const RbtCoord &coord1,
                                   const RbtCoord &coord2) {
    return RbtCoord(coord1.xyz.array() * coord2.xyz.array());
  }

  // Scalar division (coord / const)
  inline friend RbtCoord operator/(const RbtCoord &coord, const double &d) {
    return RbtCoord(coord.xyz / d);
  }

  ///////////////////////////////////////////////
  // Public methods
  ////////////////

  // Mostly applicable to vectors (rather than coords)

  // Returns square of magnitude of vector (or distance from origin)
  inline double Length2() const { return xyz.squaredNorm(); }

  // Returns magnitude of vector (or distance from origin)
  inline double Length() const { return xyz.norm(); }

  // Returns unit vector in same direction
  // Member function (V2 = V1.Unit())
  // Checks for zero length
  inline RbtCoord Unit() const {
    double l = Length();
    return (l > 0) ? *this / l : *this;
  }

  // Cross product member function (V3 = V1.Cross(V2))
  inline RbtCoord Cross(const RbtCoord &v2) const {
    return RbtCoord(xyz.cross(v2.xyz));
  }

  // Dot product member function (D = V1.Dot(V2)
  inline double Dot(const RbtCoord &v2) const { return xyz.dot(v2.xyz); }
};

// Useful typedefs

// DM 4 Mar 1999 - to make it easier to distinguish arguments that are
// treated as -3D vectors from those are treated as 3-D points, use this typedef
typedef RbtCoord RbtVector;

typedef std::vector<RbtCoord> RbtCoordList; // Vector of coords
typedef RbtCoordList::iterator RbtCoordListIter;
typedef RbtCoordList::const_iterator RbtCoordListConstIter;

// DM 19/1/99 typedefs for a map of key=RbtInt, value=RbtCoord
// typedef std::map<RbtString,RbtCoord> RbtStringCoordMap;
// typedef RbtStringCoordMap::iterator RbtStringCoordMapIter;
// typedef RbtStringCoordMap::const_iterator RbtStringCoordMapConstIter;
// DM 08 Feb 1999 - Changed from key=string to key=int
typedef std::map<unsigned int, RbtCoord> RbtUIntCoordMap;
typedef RbtUIntCoordMap::iterator RbtUIntCoordMapIter;
typedef RbtUIntCoordMap::const_iterator RbtUIntCoordMapConstIter;

///////////////////////////////////////////////
// Non-member functions (in rxdock namespace)
//////////////////////////////////////////

//////////////////
// VECTOR FUNCTIONS
//////////////////

// Returns square of magnitude of vector
inline double Length2(const RbtVector &v1) { return v1.Length2(); }

// Returns square of distance between two coords
inline double Length2(const RbtCoord &c1, const RbtCoord &c2) {
  return Length2(c2 - c1);
}

// Returns magnitude of vector
inline double Length(const RbtVector &v1) { return v1.Length(); }

// Returns distance between two coords
inline double Length(const RbtCoord &c1, const RbtCoord &c2) {
  return Length(c2 - c1);
}

// Returns unit vector in same direction
// V2 = Unit(V1)
inline RbtVector Unit(const RbtVector &v1) { return v1.Unit(); }

// Cross product (V3 = Cross(V1,V2))
inline RbtVector Cross(const RbtVector &v11, const RbtVector &v2) {
  return v11.Cross(v2);
}

// Dot product (D = Dot(V1,V2))
inline double Dot(const RbtVector &v1, const RbtVector &v2) {
  return v1.Dot(v2);
}

// Returns angle formed between 2 vectors
inline double Angle(const RbtVector &v1, const RbtVector &v2) {
  // Calculate the product of the lengths
  double d1d2 = v1.Length() * v2.Length();

  // Calculate the sin and cos
  double cos_theta = v1.Dot(v2) / d1d2;
  double sin_theta = (v1.Cross(v2)).Length() / d1d2;

  // Get theta and convert to degrees
  double theta = std::atan2(sin_theta, cos_theta);
  return theta * 180.0 / M_PI;
}

// Returns angle formed between 3 coords coord1..coord2..coord3
inline double Angle(const RbtCoord &coord1, const RbtCoord &coord2,
                    const RbtCoord &coord3) {
  // Determine vectors between coords and call vector version of Angle()
  return Angle(coord1 - coord2, coord3 - coord2);
}

// DM 7 June 1999
// Returns dihedral formed between 3 vectors
inline double Dihedral(const RbtVector &v1, const RbtVector &v2,
                       const RbtVector &v3) {
  // Calculate the cross products
  RbtVector A = v1.Cross(v2);
  RbtVector B = v2.Cross(v3);
  RbtVector C = v2.Cross(A);
  // Calculate the distances
  double rA = A.Length();
  double rB = B.Length();
  double rC = C.Length();
  // Calculate the sin and cos
  double cos_phi = A.Dot(B) / (rA * rB);
  double sin_phi = C.Dot(B) / (rC * rB);
  // Get phi and convert to degrees
  double phi = -std::atan2(sin_phi, cos_phi);
  return phi * 180.0 / M_PI;
}

// Returns dihedral formed between 4 coords c1..c2..c3..c4
inline double Dihedral(const RbtCoord &c1, const RbtCoord &c2,
                       const RbtCoord &c3, const RbtCoord &c4) {
  // Determine vectors between coords and call vector version of Dihedral()
  return Dihedral(c1 - c2, c2 - c3, c3 - c4);
}

//////////////////////
// COORDINATE FUNCTIONS
//////////////////////

// DM 28 Jul 1999 - functions for use by STL transform
inline double ExtractXCoord(const RbtCoord &c) { return c.xyz(0); }
inline double ExtractYCoord(const RbtCoord &c) { return c.xyz(1); }
inline double ExtractZCoord(const RbtCoord &c) { return c.xyz(2); }

// Returns minimum of two coords (component-wise minimum)
inline RbtCoord Min(const RbtCoord &coord1, const RbtCoord &coord2) {
  return RbtCoord(coord1.xyz.array().min(coord2.xyz.array()));
}

// Returns maximum of two coords (component-wise maximum)
inline RbtCoord Max(const RbtCoord &coord1, const RbtCoord &coord2) {
  return RbtCoord(coord1.xyz.array().max(coord2.xyz.array()));
}

// DM 28 Jul 1999 - returns component-wise min coord for coord list (i.e.min x,
// min y, min z)
inline RbtCoord Min(const RbtCoordList &cl) {
  std::vector<double> cx, cy, cz;
  // Extract arrays of x,y,z components separately
  std::transform(cl.begin(), cl.end(), std::back_inserter(cx), ExtractXCoord);
  std::transform(cl.begin(), cl.end(), std::back_inserter(cy), ExtractYCoord);
  std::transform(cl.begin(), cl.end(), std::back_inserter(cz), ExtractZCoord);
  // Compose the min coord from the min of each of the components
  return RbtCoord(*(std::min_element(cx.begin(), cx.end())),
                  *(std::min_element(cy.begin(), cy.end())),
                  *(std::min_element(cz.begin(), cz.end())));
}

// DM 28 Jul 1999 - returns component-wise max coord for coord list (i.e. max x,
// max y, max z)
inline RbtCoord Max(const RbtCoordList &cl) {
  std::vector<double> cx, cy, cz;
  // Extract arrays of x,y,z components separately
  std::transform(cl.begin(), cl.end(), std::back_inserter(cx), ExtractXCoord);
  std::transform(cl.begin(), cl.end(), std::back_inserter(cy), ExtractYCoord);
  std::transform(cl.begin(), cl.end(), std::back_inserter(cz), ExtractZCoord);
  // Compose the max coord from the max of each of the components
  return RbtCoord(*(std::max_element(cx.begin(), cx.end())),
                  *(std::max_element(cy.begin(), cy.end())),
                  *(std::max_element(cz.begin(), cz.end())));
}

// Returns center of "mass" of coords in coord list
inline RbtCoord GetCenterOfAtomicMass(const RbtCoordList &coordList) {
  RbtCoord com;
  com = std::accumulate(coordList.begin(), coordList.end(), com);
  // Check for divide by zero
  return coordList.empty() ? com : com / static_cast<double>(coordList.size());
}

} // namespace rxdock

#endif //_RBTCOORD_H_
