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
// Vector typedef can be used in argument lists to distinguish the two cases

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

#include "Config.h"

extern std::istream &eatSeps(std::istream &is);

namespace rxdock {

class Coord {
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
  inline Coord() : xyz(0.0, 0.0, 0.0) {}

  // Constructor with initial values
  inline Coord(double x1, double y1, double z1) : xyz(x1, y1, z1) {}

  // Constructor using Eigen
  inline Coord(Eigen::Vector3d xyz1) : xyz(xyz1) {}

  // Destructor
  virtual ~Coord() {}

  // Copy constructor
  inline Coord(const Coord &coord) { xyz = coord.xyz; }

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
  inline Coord &operator=(const Coord &coord) {
    if (this != &coord) { // beware of self-assignment
      xyz = coord.xyz;
    }
    return *this;
  }

  // Copy assignment (*this = const)
  inline Coord &operator=(const double &d) {
    xyz.setConstant(d);
    return *this;
  }

  //*this += coord
  inline void operator+=(const Coord &coord) { xyz += coord.xyz; }

  //*this += const
  inline void operator+=(const double &d) { xyz += xyz.Constant(d); }

  //*this -= coord
  inline void operator-=(const Coord &coord) { xyz -= coord.xyz; }

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
  friend std::ostream &operator<<(std::ostream &s, const Coord &coord) {
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

  friend std::istream &operator>>(std::istream &s, Coord &coord) {
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
      coord = Coord(x, y, z);
    return s;
  }

  // Equality
  inline friend bool operator==(const Coord &coord1, const Coord &coord2) {
    return coord1.xyz == coord2.xyz;
  }

  // Non-equality
  inline friend bool operator!=(const Coord &coord1, const Coord &coord2) {
    return coord1.xyz != coord2.xyz;
  }

  // Greater than
  inline friend bool operator>(const Coord &coord1, const Coord &coord2) {
    return (coord1.xyz.array() > coord2.xyz.array()).all();
  }

  // Greater than or equal
  inline friend bool operator>=(const Coord &coord1, const Coord &coord2) {
    return (coord1.xyz.array() >= coord2.xyz.array()).all();
  }

  // Less than
  inline friend bool operator<(const Coord &coord1, const Coord &coord2) {
    return (coord1.xyz.array() < coord2.xyz.array()).all();
  }

  // Less than or equal
  inline friend bool operator<=(const Coord &coord1, const Coord &coord2) {
    return (coord1.xyz.array() <= coord2.xyz.array()).all();
  }

  // Addition (coord + coord)
  inline friend Coord operator+(const Coord &coord1, const Coord &coord2) {
    return Coord(coord1.xyz + coord2.xyz);
  }

  // Addition (coord + const)
  inline friend Coord operator+(const Coord &coord1, double d) {
    return Coord(coord1.xyz + coord1.xyz.Constant(d));
  }

  // Addition (const + coord)
  inline friend Coord operator+(double d, const Coord &coord1) {
    return Coord(coord1.xyz.Constant(d) + coord1.xyz);
  }

  // Subtraction (coord - coord)
  inline friend Coord operator-(const Coord &coord1, const Coord &coord2) {
    return Coord(coord1.xyz - coord2.xyz);
  }

  // Subtraction (coord - const)
  inline friend Coord operator-(const Coord &coord1, double d) {
    return Coord(coord1.xyz - coord1.xyz.Constant(d));
  }

  // Subtraction (const - coord)
  inline friend Coord operator-(double d, const Coord &coord1) {
    return Coord(coord1.xyz.Constant(d) - coord1.xyz);
  }

  // Negation
  inline friend Coord operator-(const Coord &coord) {
    return Coord(-coord.xyz);
  }

  // Dot product
  // DM 7 Jan 1999 - operator* replaced by Dot(coord) member function
  // and by Dot(coord1,coord2) non-member function
  // inline friend Double operator*(const Coord& coord1, const Coord&
  // coord2) {
  //  return coord1.xyz.dot(coord.xyz);
  //}

  // Scalar product (coord * const)
  inline friend Coord operator*(const Coord &coord, const double &d) {
    return Coord(coord.xyz * d);
  }

  // Scalar product (const * coord)
  inline friend Coord operator*(const double &d, const Coord &coord) {
    return Coord(d * coord.xyz);
  }

  // Scalar product (coord * coord : component-wise multiplication)
  inline friend Coord operator*(const Coord &coord1, const Coord &coord2) {
    return Coord(coord1.xyz.array() * coord2.xyz.array());
  }

  // Scalar division (coord / const)
  inline friend Coord operator/(const Coord &coord, const double &d) {
    return Coord(coord.xyz / d);
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
  inline Coord Unit() const {
    double l = Length();
    return (l > 0) ? *this / l : *this;
  }

  // Cross product member function (V3 = V1.Cross(V2))
  inline Coord Cross(const Coord &v2) const { return Coord(xyz.cross(v2.xyz)); }

  // Dot product member function (D = V1.Dot(V2)
  inline double Dot(const Coord &v2) const { return xyz.dot(v2.xyz); }
};

// Useful typedefs

// DM 4 Mar 1999 - to make it easier to distinguish arguments that are
// treated as -3D vectors from those are treated as 3-D points, use this typedef
typedef Coord Vector;

typedef std::vector<Coord> CoordList; // Vector of coords
typedef CoordList::iterator CoordListIter;
typedef CoordList::const_iterator CoordListConstIter;

// DM 19/1/99 typedefs for a map of key=Int, value=Coord
// typedef std::map<String,Coord> StringCoordMap;
// typedef StringCoordMap::iterator StringCoordMapIter;
// typedef StringCoordMap::const_iterator StringCoordMapConstIter;
// DM 08 Feb 1999 - Changed from key=string to key=int
typedef std::map<unsigned int, Coord> UIntCoordMap;
typedef UIntCoordMap::iterator UIntCoordMapIter;
typedef UIntCoordMap::const_iterator UIntCoordMapConstIter;

///////////////////////////////////////////////
// Non-member functions (in rxdock namespace)
//////////////////////////////////////////

//////////////////
// VECTOR FUNCTIONS
//////////////////

// Returns square of magnitude of vector
inline double Length2(const Vector &v1) { return v1.Length2(); }

// Returns square of distance between two coords
inline double Length2(const Coord &c1, const Coord &c2) {
  return Length2(c2 - c1);
}

// Returns magnitude of vector
inline double Length(const Vector &v1) { return v1.Length(); }

// Returns distance between two coords
inline double Length(const Coord &c1, const Coord &c2) {
  return Length(c2 - c1);
}

// Returns unit vector in same direction
// V2 = Unit(V1)
inline Vector Unit(const Vector &v1) { return v1.Unit(); }

// Cross product (V3 = Cross(V1,V2))
inline Vector Cross(const Vector &v11, const Vector &v2) {
  return v11.Cross(v2);
}

// Dot product (D = Dot(V1,V2))
inline double Dot(const Vector &v1, const Vector &v2) { return v1.Dot(v2); }

// Returns angle formed between 2 vectors
inline double Angle(const Vector &v1, const Vector &v2) {
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
inline double Angle(const Coord &coord1, const Coord &coord2,
                    const Coord &coord3) {
  // Determine vectors between coords and call vector version of Angle()
  return Angle(coord1 - coord2, coord3 - coord2);
}

// DM 7 June 1999
// Returns dihedral formed between 3 vectors
inline double Dihedral(const Vector &v1, const Vector &v2, const Vector &v3) {
  // Calculate the cross products
  Vector A = v1.Cross(v2);
  Vector B = v2.Cross(v3);
  Vector C = v2.Cross(A);
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
inline double Dihedral(const Coord &c1, const Coord &c2, const Coord &c3,
                       const Coord &c4) {
  // Determine vectors between coords and call vector version of Dihedral()
  return Dihedral(c1 - c2, c2 - c3, c3 - c4);
}

//////////////////////
// COORDINATE FUNCTIONS
//////////////////////

// DM 28 Jul 1999 - functions for use by STL transform
inline double ExtractXCoord(const Coord &c) { return c.xyz(0); }
inline double ExtractYCoord(const Coord &c) { return c.xyz(1); }
inline double ExtractZCoord(const Coord &c) { return c.xyz(2); }

// Returns minimum of two coords (component-wise minimum)
inline Coord Min(const Coord &coord1, const Coord &coord2) {
  return Coord(coord1.xyz.array().min(coord2.xyz.array()));
}

// Returns maximum of two coords (component-wise maximum)
inline Coord Max(const Coord &coord1, const Coord &coord2) {
  return Coord(coord1.xyz.array().max(coord2.xyz.array()));
}

// DM 28 Jul 1999 - returns component-wise min coord for coord list (i.e.min x,
// min y, min z)
inline Coord Min(const CoordList &cl) {
  std::vector<double> cx, cy, cz;
  // Extract arrays of x,y,z components separately
  std::transform(cl.begin(), cl.end(), std::back_inserter(cx), ExtractXCoord);
  std::transform(cl.begin(), cl.end(), std::back_inserter(cy), ExtractYCoord);
  std::transform(cl.begin(), cl.end(), std::back_inserter(cz), ExtractZCoord);
  // Compose the min coord from the min of each of the components
  return Coord(*(std::min_element(cx.begin(), cx.end())),
               *(std::min_element(cy.begin(), cy.end())),
               *(std::min_element(cz.begin(), cz.end())));
}

// DM 28 Jul 1999 - returns component-wise max coord for coord list (i.e. max x,
// max y, max z)
inline Coord Max(const CoordList &cl) {
  std::vector<double> cx, cy, cz;
  // Extract arrays of x,y,z components separately
  std::transform(cl.begin(), cl.end(), std::back_inserter(cx), ExtractXCoord);
  std::transform(cl.begin(), cl.end(), std::back_inserter(cy), ExtractYCoord);
  std::transform(cl.begin(), cl.end(), std::back_inserter(cz), ExtractZCoord);
  // Compose the max coord from the max of each of the components
  return Coord(*(std::max_element(cx.begin(), cx.end())),
               *(std::max_element(cy.begin(), cy.end())),
               *(std::max_element(cz.begin(), cz.end())));
}

// Returns center of "mass" of coords in coord list
inline Coord GetCenterOfAtomicMass(const CoordList &coordList) {
  Coord com;
  com = std::accumulate(coordList.begin(), coordList.end(), com);
  // Check for divide by zero
  return coordList.empty() ? com : com / static_cast<double>(coordList.size());
}

} // namespace rxdock

#endif //_RBTCOORD_H_
