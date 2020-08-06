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

// Simple class for handling quaternions

#ifndef _RBTQUAT_H_
#define _RBTQUAT_H_

#include <cmath> //for sqrt

#include "Coord.h"

namespace rxdock {

class Quat {
  ///////////////////////////////////////////////
  // Data members
  ///////////////
  //(leave public so we don't need
  // to write accessor functions)
public:
  double s; // Scalar component
  Vector v; // Vector component

  ///////////////////////////////////////////////
  // Constructors/destructors:
  ///////////////////////////
public:
  // Default constructor (corresponds to no rotation (null transformation))
  inline Quat() : s(1.0), v(0.0, 0.0, 0.0) {}

  // Constructor with initial values (vector component passed as 3 Double's)
  inline Quat(double s1, double vx, double vy, double vz)
      : s(s1), v(vx, vy, vz) {}

  // Constructor with initial values (vector component passed as Vector)
  // CAUTION: Argument types are the reverse of below
  inline Quat(double s1, const Vector &v1) : s(s1), v(v1) {}

  // Constructor accepting a rotation axis and rotation angle (radians).
  // CAUTION: Argument types are the reverse of above
  inline Quat(const Vector &axis, double phi) {
    double halfPhi(0.5 * phi);
    s = std::cos(halfPhi);
    v = std::sin(halfPhi) * axis.Unit();
  }

  // Destructor
  virtual ~Quat() {}

  // Copy constructor
  inline Quat(const Quat &quat) {
    s = quat.s;
    v = quat.v;
  }

  ///////////////////////////////////////////////
  // Operator functions:
  /////////////////////

  // Copy assignment (*this = coord)
  inline Quat &operator=(const Quat &quat) {
    if (this != &quat) { // beware of self-assignment
      s = quat.s;
      v = quat.v;
    }
    return *this;
  }

  //*this += quat
  inline void operator+=(const Quat &quat) {
    s += quat.s;
    v += quat.v;
  }

  //*this -= quat
  inline void operator-=(const Quat &quat) {
    s -= quat.s;
    v -= quat.v;
  }

  //*this *= const
  inline void operator*=(const double &d) {
    s *= d;
    v *= d;
  }

  //*this /= const
  inline void operator/=(const double &d) {
    s /= d;
    v /= d;
  }

  ///////////////////////////////////////////////
  // Friend functions:
  ///////////////////

  // Insertion operator
  friend std::ostream &operator<<(std::ostream &s, const Quat &quat) {
    return s << quat.s << "," << quat.v;
  }

  // Equality
  inline friend bool operator==(const Quat &quat1, const Quat &quat2) {
    return quat1.s == quat2.s && quat1.v == quat2.v;
  }

  // Non-equality
  inline friend bool operator!=(const Quat &quat1, const Quat &quat2) {
    return quat1.s != quat2.s || quat1.v != quat2.v;
  }

  // Addition
  inline friend Quat operator+(const Quat &quat1, const Quat &quat2) {
    return Quat(quat1.s + quat2.s, quat1.v + quat2.v);
  }

  // Subtraction
  inline friend Quat operator-(const Quat &quat1, const Quat &quat2) {
    return Quat(quat1.s - quat2.s, quat1.v - quat2.v);
  }

  // Negation
  inline friend Quat operator-(const Quat &quat) {
    return Quat(-quat.s, -quat.v);
  }

  // Multiplication (non-commutative)
  inline friend Quat operator*(const Quat &quat1, const Quat &quat2) {
    return Quat(quat1.s * quat2.s - (quat1.v).Dot(quat2.v),
                quat1.s * quat2.v + quat1.v * quat2.s +
                    (quat1.v).Cross(quat2.v));
  }

  // Scalar product
  inline friend Quat operator*(const Quat &quat, const double &d) {
    return Quat(quat.s * d, quat.v * d);
  }

  // Scalar product
  inline friend Quat operator*(const double &d, const Quat &quat) {
    return Quat(quat.s * d, quat.v * d);
  }

  // Scalar division
  inline friend Quat operator/(const Quat &quat, const double &d) {
    return Quat(quat.s / d, quat.v / d);
  }

  ///////////////////////////////////////////////
  // Public methods
  ////////////////

  // Returns magnitude of quat
  // Member function (L = Q1.Length())
  inline double Length() const { return std::sqrt(s * s + v.Length2()); }

  // Returns unit quat
  // Member function (U1 = Q1.Unit())
  inline Quat Unit() const {
    double l = Length();
    return (l > 0) ? *this / l : *this;
  }

  // Dot product member function (D = Q1.Dot(Q2))
  inline double Dot(const Quat &quat) const {
    return s * quat.s + v.Dot(quat.v);
  }

  // Returns conjugate
  // Member function (Q1* = Q1.Conj())
  inline Quat Conj() const { return Quat(s, -v); }

  // Returns Q W Q* where Q=S+V is a quaternion and W is a 3-D vector
  // Q W Q* = (S^2 - V.V)W + 2S(VxW) + 2V(V.W)
  // where . = Dot product and x = Cross product
  // The end result is a 3-D rotation, assuming q has been obtained from
  // RotationQuat
  inline Coord Rotate(const Coord &w) const {
    return (s * s - v.Dot(v)) * w + 2 * s * v.Cross(w) + 2 * v * v.Dot(w);
  }
};

// Useful typedefs
typedef std::vector<Quat> QuatList; // Vector of quats
typedef QuatList::iterator QuatListIter;
typedef QuatList::const_iterator QuatListConstIter;

///////////////////////////////////////////////
// Non-member functions (in rxdock namespace)
//////////////////////////////////////////

// Returns magnitude of quat (L = Length(Q1))
inline double Length(const Quat &quat) { return quat.Length(); }

// Returns unit quat (U1 = Unit(Q1))
inline Quat Unit(const Quat &quat) { return quat.Unit(); }

// Dot product (D = Dot(Q1,Q2))
inline double Dot(const Quat &quat1, const Quat &quat2) {
  return quat1.Dot(quat2);
}

// Returns conjugate (Q1* = Conj(Q1))
inline Quat Conj(const Quat &quat) { return quat.Conj(); }

} // namespace rxdock

#endif //_RBTQUAT_H_
