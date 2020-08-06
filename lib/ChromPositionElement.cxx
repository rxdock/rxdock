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

#include "ChromPositionElement.h"

using namespace rxdock;

std::string ChromPositionElement::_CT = "ChromPositionElement";

ChromPositionElement::ChromPositionElement(
    const Model *pModel, const DockingSite *pDockSite, double transStepSize,
    double rotStepSize, ChromElement::eMode transMode,
    ChromElement::eMode rotMode, double maxTrans, double maxRot) {
  m_spRefData =
      new ChromPositionRefData(pModel, pDockSite, transStepSize, rotStepSize,
                               transMode, rotMode, maxTrans, maxRot);
  SyncFromModel();
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

ChromPositionElement::ChromPositionElement(ChromPositionRefDataPtr spRefData,
                                           const Coord &com,
                                           const Euler &orientation)
    : m_spRefData(spRefData), m_com(com), m_orientation(orientation) {
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

ChromPositionElement::~ChromPositionElement() { _RBTOBJECTCOUNTER_DESTR_(_CT); }

void ChromPositionElement::Reset() {
  m_com = m_spRefData->GetInitialCOM();
  m_orientation = m_spRefData->GetInitialOrientation();
}

void ChromPositionElement::Randomise() {
  RandomiseCOM();
  RandomiseOrientation();
}

void ChromPositionElement::RandomiseCOM() {
  double dist;
  Vector axis;
  int numStartCoords;
  int iCoord;
  switch (m_spRefData->GetTransMode()) {
    // TETHERED: Perform a single mutation from the initial COM
    // up to the maximum permitted
  case ChromElement::TETHERED:
    dist = m_spRefData->GetMaxTrans() * GetRand().GetRandom01();
    axis = GetRand().GetRandomUnitVector();
    m_com = m_spRefData->GetInitialCOM() + dist * axis;
    break;
  // FREE: Randomise across docking volume coords
  case ChromElement::FREE:
    numStartCoords = m_spRefData->GetNumStartCoords();
    iCoord = GetRand().GetRandomInt(numStartCoords);
    m_com = m_spRefData->GetStartCoord(iCoord);
    break;
  // FIXED: Revert to initial COM
  default:
    m_com = m_spRefData->GetInitialCOM();
    break;
  }
}

void ChromPositionElement::RandomiseOrientation() {
  double theta;
  Vector axis;
  double heading, attitude, bank;
  switch (m_spRefData->GetRotMode()) {
    // TETHERED: Perform a single mutation from the initial orientation
    // up to the maximum permitted
  case ChromElement::TETHERED:
    m_orientation = m_spRefData->GetInitialOrientation();
    theta = m_spRefData->GetMaxRot() * GetRand().GetRandom01();
    axis = GetRand().GetRandomUnitVector();
    m_orientation.Rotate(axis, theta);
    break;
  // FREE: completely scramble the initial orientation
  case ChromElement::FREE:
    heading = 2.0 * M_PI * GetRand().GetRandom01() - M_PI;
    attitude = M_PI * GetRand().GetRandom01() - 0.5 * M_PI;
    bank = 2.0 * M_PI * GetRand().GetRandom01() - M_PI;
    m_orientation = Euler(heading, attitude, bank);
    break;
  // FIXED: Revert to initial orientation
  default:
    m_orientation = m_spRefData->GetInitialOrientation();
    break;
  }
}

void ChromPositionElement::Mutate(double relStepSize) {
  MutateCOM(relStepSize);
  MutateOrientation(relStepSize);
}

void ChromPositionElement::MutateCOM(double relStepSize) {
  double absTransStepSize;
  double dist;
  Vector axis;
  switch (m_spRefData->GetTransMode()) {
  case ChromElement::TETHERED:
    absTransStepSize = relStepSize * m_spRefData->GetTransStepSize();
    if (absTransStepSize > 0) {
      dist = absTransStepSize * GetRand().GetRandom01();
      axis = GetRand().GetRandomUnitVector();
      m_com += dist * axis;
      CorrectTetheredCOM();
    }
    break;
  case ChromElement::FREE:
    absTransStepSize = relStepSize * m_spRefData->GetTransStepSize();
    if (absTransStepSize > 0) {
      dist = absTransStepSize * GetRand().GetRandom01();
      axis = GetRand().GetRandomUnitVector();
      m_com += dist * axis;
    }
    break;
  // FIXED: Do nothing
  default:
    break;
  }
}

void ChromPositionElement::MutateOrientation(double relStepSize) {
  double absRotStepSize;
  double theta;
  Vector axis;
  switch (m_spRefData->GetRotMode()) {
  case ChromElement::TETHERED:
    absRotStepSize = relStepSize * m_spRefData->GetRotStepSize();
    if (absRotStepSize > 0) {
      theta = absRotStepSize * GetRand().GetRandom01();
      axis = GetRand().GetRandomUnitVector();
      m_orientation.Rotate(axis, theta);
      CorrectTetheredOrientation();
    }
    break;
  case ChromElement::FREE:
    absRotStepSize = relStepSize * m_spRefData->GetRotStepSize();
    if (absRotStepSize > 0) {
      theta = absRotStepSize * GetRand().GetRandom01();
      axis = GetRand().GetRandomUnitVector();
      m_orientation.Rotate(axis, theta);
    }
    break;
  // FIXED: Do nothing
  default:
    break;
  }
}

void ChromPositionElement::SyncFromModel() {
  m_spRefData->GetModelValue(m_com, m_orientation);
}

void ChromPositionElement::SyncToModel() {
  m_spRefData->SetModelValue(m_com, m_orientation);
}

ChromElement *ChromPositionElement::clone() const {
  return new ChromPositionElement(m_spRefData, m_com, m_orientation);
}

void ChromPositionElement::GetVector(std::vector<double> &v) const {
  if (!m_spRefData->IsTransFixed()) {
    v.insert(v.end(), m_com.xyz.data(), m_com.xyz.data() + m_com.xyz.size());
  }
  if (!m_spRefData->IsRotFixed()) {
    v.push_back(m_orientation.GetHeading());
    v.push_back(m_orientation.GetAttitude());
    v.push_back(m_orientation.GetBank());
  }
}

void ChromPositionElement::GetVector(XOverList &v) const {
  if (!m_spRefData->IsTransFixed()) {
    XOverElement comElement(m_com.xyz.data(),
                            m_com.xyz.data() + m_com.xyz.size());
    v.push_back(comElement);
  }
  if (!m_spRefData->IsRotFixed()) {
    XOverElement orientationElement;
    orientationElement.push_back(m_orientation.GetHeading());
    orientationElement.push_back(m_orientation.GetAttitude());
    orientationElement.push_back(m_orientation.GetBank());
    v.push_back(orientationElement);
  }
}

void ChromPositionElement::SetVector(const std::vector<double> &v, int &i) {
  if (VectorOK(v, i)) {
    if (!m_spRefData->IsTransFixed()) {
      // 2013Nov26 (DM) Bug fix to unsafe code.
      // We cannot assume that the multiple increments will be processed
      // left->right. This assumption is broken in g++ 4. Safer to break into
      // separate statements. m_com = Coord(v[i++], v[i++], v[i++]);
      double x(v[i++]);
      double y(v[i++]);
      double z(v[i++]);
      m_com = Coord(x, y, z);
    }
    if (!m_spRefData->IsRotFixed()) {
      // 2013Nov26 (DM) Bug fix to unsafe code.
      // We cannot assume that the multiple increments will be processed
      // left->right. This assumption is broken in g++ 4. Safer to break into
      // separate statements. m_orientation = Euler(v[i++], v[i++], v[i++]);
      double heading(v[i++]);
      double attitude(v[i++]);
      double bank(v[i++]);
      m_orientation = Euler(heading, attitude, bank);
      m_orientation.Standardise();
    }
  } else {
    throw BadArgument(_WHERE_,
                      "Index out of range or insufficient elements remaining");
  }
}

void ChromPositionElement::SetVector(const XOverList &v, int &i) {
  if (VectorOK(v, i)) {
    if (!m_spRefData->IsTransFixed()) {
      XOverElement comElement(v[i++]);
      if (comElement.size() == 3) {
        // As we crossover an intact COM vector there should be no need to check
        // for tethered bounds
        m_com = Coord(comElement[0], comElement[1], comElement[2]);
      } else {
        throw BadArgument(_WHERE_, "comElement vector is of incorrect length");
      }
    }
    if (!m_spRefData->IsRotFixed()) {
      XOverElement orientationElement(v[i++]);
      if (orientationElement.size() == 3) {
        // As we crossover an intact orientation vector there should be no need
        // to check for tethered bounds
        m_orientation = Euler(orientationElement[0], orientationElement[1],
                              orientationElement[2]);
      } else {
        throw BadArgument(_WHERE_,
                          "orientationElement vector is of incorrect length");
      }
    }
  } else {
    throw BadArgument(_WHERE_,
                      "Index out of range or insufficient elements remaining");
  }
}

void ChromPositionElement::GetStepVector(std::vector<double> &v) const {
  if (!m_spRefData->IsTransFixed()) {
    double transStepSize = m_spRefData->GetTransStepSize();
    for (int i = 0; i < 3; ++i) {
      v.push_back(transStepSize);
    }
  }
  if (!m_spRefData->IsRotFixed()) {
    double rotStepSize = m_spRefData->GetRotStepSize();
    for (int i = 0; i < 3; ++i) {
      v.push_back(rotStepSize);
    }
  }
}

double ChromPositionElement::CompareVector(const std::vector<double> &v,
                                           int &i) const {
  double retVal(0.0);
  if (VectorOK(v, i)) {
    if (!m_spRefData->IsTransFixed()) {
      // 2013Nov26 (DM) Bug fix to unsafe code.
      // We cannot assume that the multiple increments will be processed
      // left->right. This assumption is broken in g++ 4. Safer to break into
      // separate statements. Coord otherCom = Coord(v[i++], v[i++],
      // v[i++]);
      double x(v[i++]);
      double y(v[i++]);
      double z(v[i++]);
      Coord otherCom = Coord(x, y, z);
      double transStepSize = m_spRefData->GetTransStepSize();
      // Compare distance between centres of mass
      if (transStepSize > 0.0) {
        double absDiff = Length(m_com, otherCom);
        double relDiff = absDiff / transStepSize;
        retVal = std::max(retVal, relDiff);
      }
    }
    if (!m_spRefData->IsRotFixed()) {
      // 2013Nov26 (DM) Bug fix to unsafe code.
      // We cannot assume that the multiple increments will be processed
      // left->right. This assumption is broken in g++ 4. Safer to break into
      // separate statements. Euler otherOrientation = Euler(v[i++],
      // v[i++], v[i++]);
      double heading(v[i++]);
      double attitude(v[i++]);
      double bank(v[i++]);
      Euler otherOrientation = Euler(heading, attitude, bank);
      double rotStepSize = m_spRefData->GetRotStepSize();
      // Compare orientations
      if (rotStepSize > 0.0) {
        // Determine the difference between the two orientations
        // in terms of the axis/angle needed to align them
        // q.s = std::cos(phi / 2)
        Quat qAlign = otherOrientation.ToQuat() * m_orientation.ToQuat().Conj();
        double cosHalfTheta = qAlign.s;
        if (cosHalfTheta < -1.0) {
          cosHalfTheta = -1.0;
        } else if (cosHalfTheta > 1.0) {
          cosHalfTheta = 1.0;
        }
        double absDiff =
            std::fabs(StandardisedValue(2.0 * std::acos(cosHalfTheta)));
        double relDiff = absDiff / rotStepSize;
        retVal = std::max(retVal, relDiff);
      }
    }
  } else {
    retVal = -1.0;
  }
  return retVal;
}

void ChromPositionElement::Print(std::ostream &s) const {
  s << "COM " << m_com << std::endl;
  s << "EULER " << m_orientation << std::endl;
}

double ChromPositionElement::StandardisedValue(double rotationAngle) {
  while (rotationAngle >= M_PI) {
    rotationAngle -= 2.0 * M_PI;
  }
  while (rotationAngle < -M_PI) {
    rotationAngle += 2.0 * M_PI;
  }
  return rotationAngle;
}

void ChromPositionElement::CorrectTetheredCOM() {
  // If we are out of bounds
  // revert to the initial COM and translate as far as we can
  // towards the new point until we hit the boundary
  double maxTrans = m_spRefData->GetMaxTrans();
  // vector from initial COM to new point
  Vector axis = m_com - m_spRefData->GetInitialCOM();
  // avoid square roots until necessary
  if (axis.Length2() > (maxTrans * maxTrans)) {
    // Stay just inside the boundary sphere
    m_com = m_spRefData->GetInitialCOM() + (0.999 * maxTrans * axis.Unit());
  }
}

void ChromPositionElement::CorrectTetheredOrientation() {
  // Check for orientation out of bounds
  double maxRot = m_spRefData->GetMaxRot();
  Quat qAlign = m_spRefData->GetInitialQuat() * m_orientation.ToQuat().Conj();
  double cosHalfTheta = qAlign.s;
  if (cosHalfTheta < -1.0) {
    cosHalfTheta = -1.0;
  } else if (cosHalfTheta > 1.0) {
    cosHalfTheta = 1.0;
  }
  // Theta is the rotation angle required to realign with reference
  double theta = StandardisedValue(2.0 * std::acos(cosHalfTheta));
  // Deal with pos and neg angles independently as we have to
  // invert the rotation axis for negative angles
  if (theta < -maxRot) {
    Vector axis = -qAlign.v / std::sin(theta / 2.0);
    // Adjust theta to bring the orientation just inside the tethered bound
    theta += 0.999 * maxRot;
    m_orientation.Rotate(axis, theta);
  } else if (theta > maxRot) {
    Vector axis = qAlign.v / std::sin(theta / 2.0);
    // Adjust theta to bring the orientation just inside the tethered bound
    theta -= 0.999 * maxRot;
    m_orientation.Rotate(axis, theta);
  }
}
