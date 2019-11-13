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

// Simple class to manage Euler angles
// Main purpose is to provides conversion methods
// to/from quaternion representation
//
// We adopt a modified "NASA Standard Aeroplane" convention
// in which the y and z axes are swapped relative to the definition in
// http://www.euclideanspace.com/maths/geometry/rotations/euler/index.htm
//(and related links).
// This should be equivalent to the VRML convention.
//
// An object's orientation is defined in terms of:
//
// 1) heading (-PI to +PI) (rotation around z) (y in NASA Aeroplane)
// 2) attitude (-PI/2 to +PI/2) (rotation around y) (z in NASA Aeroplane)
// 3) bank (-PI to +PI) (rotation around x)
//
// See http://www.euclideanspace.com/maths/geometry/rotations/euler/index.htm
// and related links for definitions and conversion equations.
//
// NOTE: This class is designed to provide fast conversion between Euler angles
// and quaternions. It is the callers responsibility to ensure that:
// a) heading, attitude and bank angles are within range
// b) quaternions passed to the constructor and SetQuat are of unit length
// Run-time exceptions may occur if these conditions are not met.
#ifndef RBTEULER_H_
#define RBTEULER_H_

#include "RbtQuat.h"

class RbtEuler {
public:
  // Constructor accepting values for heading, attitude and bank (radians).
  // NOTE: For performance reasons, no checks are made that the values
  // are within range
  RbtEuler(double heading = 0.0, double attitude = 0.0, double bank = 0.0)
      : m_heading(heading), m_attitude(attitude), m_bank(bank) {}
  // Constructor accepting a quaternion.
  // Quaternion is converted to Euler angle representation
  // NOTE: For performance reasons, q is assumed to be of unit length
  // and no further checks are made.
  RbtEuler(const RbtQuat &q) { FromQuat(q); }

  // Returns true if all three Euler angles are within their
  // standard ranges
  bool isStandardised() const {
    return (m_heading >= -M_PI) && (m_heading <= M_PI) &&
           (m_attitude >= -M_PI / 2.0) && (m_attitude <= M_PI / 2.0) &&
           (m_bank >= -M_PI) && (m_bank <= M_PI);
  }
  // Standardises the Euler angles by conversion to/from quaternion
  void Standardise() {
    if (!isStandardised())
      FromQuat(ToQuat());
  }

  // Gets the current heading angle (radians)
  double GetHeading() const { return m_heading; }
  // Gets the current attitude angle (radians)
  double GetAttitude() const { return m_attitude; }
  // Gets the current bank angle (radians)
  double GetBank() const { return m_bank; }
  // Gets the equivalent quaternion representation
  // NOTE: For performance reasons, no checks are made that the current
  // Euler angles are within range prior to the conversion.
  RBTDLL_EXPORT RbtQuat ToQuat() const;
  // Updates the Euler angles from a quaternion representation
  // NOTE: For performance reasons, q is assumed to be of unit length
  // and no further checks are made.
  RBTDLL_EXPORT void FromQuat(const RbtQuat &q);
  // Convenience method to rotate orientation by a quaternion
  void Rotate(const RbtQuat &q) { FromQuat(q * ToQuat()); }
  void Rotate(const RbtVector &axis, double theta) {
    Rotate(RbtQuat(axis, theta));
  }

  // Insertion operator
  friend std::ostream &operator<<(std::ostream &s, const RbtEuler &euler) {
    return s << "Heading = " << euler.m_heading
             << ", Attitude = " << euler.m_attitude
             << ", Bank = " << euler.m_bank;
  }

private:
  double m_heading;  // heading angle (radians)
  double m_attitude; // attitude angle (radians)
  double m_bank;     // bank angle (radians)
};
#endif /*RBTEULER_H_*/
