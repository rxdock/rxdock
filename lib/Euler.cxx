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

#include "rxdock/Euler.h"

using namespace rxdock;

Quat Euler::ToQuat() const {
  double c1 = std::cos(m_heading / 2.0);
  double s1 = std::sin(m_heading / 2.0);
  double c2 = std::cos(m_attitude / 2.0);
  double s2 = std::sin(m_attitude / 2.0);
  double c3 = std::cos(m_bank / 2.0);
  double s3 = std::sin(m_bank / 2.0);
  double c1c2 = c1 * c2;
  double s1s2 = s1 * s2;
  Quat q(c1c2 * c3 - s1s2 * s3, c1c2 * s3 + s1s2 * c3,
         c1 * s2 * c3 - s1 * c2 * s3, s1 * c2 * c3 + c1 * s2 * s3);
  return q;
}

void Euler::FromQuat(const Quat &q) {
  double test = (q.v.xyz(0) * q.v.xyz(2)) + (q.v.xyz(1) * q.s);
  if (test > 0.499999) { // singularity at north pole
    m_heading = 2.0 * std::atan2(q.v.xyz(0), q.s);
    m_attitude = M_PI / 2.0;
    m_bank = 0.0;
  } else if (test < -0.499999) { // singularity at south pole
    m_heading = -2.0 * std::atan2(q.v.xyz(0), q.s);
    m_attitude = -M_PI / 2.0;
    m_bank = 0.0;
  } else {
    Vector sq = q.v * q.v;
    m_heading = std::atan2(2.0 * (q.v.xyz(2) * q.s - q.v.xyz(0) * q.v.xyz(1)),
                           1.0 - 2.0 * (sq.xyz(2) + sq.xyz(1)));
    m_attitude = std::asin(2.0 * test);
    m_bank = std::atan2(2.0 * (q.v.xyz(0) * q.s - q.v.xyz(2) * q.v.xyz(1)),
                        1.0 - 2.0 * (sq.xyz(0) + sq.xyz(1)));
  }
}

void rxdock::to_json(json &j, const Euler &eul) {
  j = json{{"heading", eul.m_heading},
           {"attitude", eul.m_attitude},
           {"bank", eul.m_bank}};
}

void rxdock::from_json(const json &j, Euler &eul) {
  j.at("heading").get_to(eul.m_heading);
  j.at("attitude").get_to(eul.m_attitude);
  j.at("bank").get_to(eul.m_bank);
}