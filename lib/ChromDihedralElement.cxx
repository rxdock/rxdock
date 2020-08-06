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

#include "rxdock/ChromDihedralElement.h"

using namespace rxdock;

std::string ChromDihedralElement::_CT = "ChromDihedralElement";

ChromDihedralElement::ChromDihedralElement(BondPtr spBond,
                                           AtomList tetheredAtoms,
                                           double stepSize,
                                           ChromElement::eMode mode,
                                           double maxDihedral)
    : m_value(0.0) {
  m_spRefData = new ChromDihedralRefData(spBond, tetheredAtoms, stepSize, mode,
                                         maxDihedral);
  // Set the initial genotype to match the current phenotype
  SyncFromModel();
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

ChromDihedralElement::ChromDihedralElement(ChromDihedralRefDataPtr spRefData,
                                           double value)
    : m_spRefData(spRefData), m_value(value) {
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

ChromDihedralElement::~ChromDihedralElement() { _RBTOBJECTCOUNTER_DESTR_(_CT); }

void ChromDihedralElement::Reset() { m_value = m_spRefData->GetInitialValue(); }

void ChromDihedralElement::Randomise() {
  double maxDelta = m_spRefData->GetMaxDihedral();
  double delta;
  switch (m_spRefData->GetMode()) {
  case ChromElement::TETHERED:
    delta = 2.0 * maxDelta * GetRand().GetRandom01() - maxDelta;
    m_value = StandardisedValue(m_spRefData->GetInitialValue() + delta);
    break;
  case ChromElement::FREE:
    m_value = 360.0 * GetRand().GetRandom01() - 180.0;
    break;
  default:
    m_value = m_spRefData->GetInitialValue();
    break;
  }
}

void ChromDihedralElement::Mutate(double relStepSize) {
  double absStepSize = relStepSize * m_spRefData->GetStepSize();
  double delta;
  if (absStepSize > 0) {
    switch (m_spRefData->GetMode()) {
    case ChromElement::TETHERED:
      delta = 2.0 * absStepSize * GetRand().GetRandom01() - absStepSize;
      m_value = StandardisedValue(m_value + delta);
      CorrectTetheredDihedral();
      break;
    case ChromElement::FREE:
      delta = 2.0 * absStepSize * GetRand().GetRandom01() - absStepSize;
      m_value = StandardisedValue(m_value + delta);
      break;
    default:
      break;
    }
  }
}

void ChromDihedralElement::SyncFromModel() {
  m_value = m_spRefData->GetModelValue();
}

void ChromDihedralElement::SyncToModel() {
  m_spRefData->SetModelValue(m_value);
}

ChromElement *ChromDihedralElement::clone() const {
  return new ChromDihedralElement(m_spRefData, m_value);
}

void ChromDihedralElement::GetVector(std::vector<double> &v) const {
  v.push_back(m_value);
}

void ChromDihedralElement::GetVector(XOverList &v) const {
  XOverElement dihedralElement;
  dihedralElement.push_back(m_value);
  v.push_back(dihedralElement);
}

void ChromDihedralElement::SetVector(const std::vector<double> &v, int &i) {
  if (VectorOK(v, i)) {
    m_value = StandardisedValue(v[i++]);
  } else {
    throw BadArgument(_WHERE_,
                      "Index out of range or insufficient elements remaining");
  }
}

void ChromDihedralElement::SetVector(const XOverList &v, int &i) {
  if (VectorOK(v, i)) {
    XOverElement dihedralElement(v[i++]);
    if (dihedralElement.size() == 1) {
      m_value = dihedralElement[0];
    } else {
      throw BadArgument(_WHERE_,
                        "dihedralElement vector is of incorrect length");
    }
  } else {
    throw BadArgument(_WHERE_,
                      "Index out of range or insufficient elements remaining");
  }
}

void ChromDihedralElement::GetStepVector(std::vector<double> &v) const {
  v.push_back(m_spRefData->GetStepSize());
}

double ChromDihedralElement::CompareVector(const std::vector<double> &v,
                                           int &i) const {
  double retVal(0.0);
  if (!VectorOK(v, i)) {
    retVal = -1.0;
  } else {
    double otherAngle = v[i++];
    double stepSize = m_spRefData->GetStepSize();
    if (stepSize > 0.0) {
      double absDiff = std::fabs(StandardisedValue(m_value - otherAngle));
      retVal = absDiff / stepSize;
    }
  }
  return retVal;
}

void ChromDihedralElement::Print(std::ostream &s) const {
  s << "DIHEDRAL " << m_value << std::endl;
}

double ChromDihedralElement::StandardisedValue(double dihedralAngle) {
  while (dihedralAngle >= 180) {
    dihedralAngle -= 360.0;
  }
  while (dihedralAngle < -180.0) {
    dihedralAngle += 360.0;
  }
  return dihedralAngle;
}

void ChromDihedralElement::CorrectTetheredDihedral() {
  double maxDelta = m_spRefData->GetMaxDihedral();
  double initialValue = m_spRefData->GetInitialValue();
  double delta = StandardisedValue(m_value - initialValue);
  if (delta > maxDelta) {
    m_value = StandardisedValue(initialValue + maxDelta);
  } else if (delta < -maxDelta) {
    m_value = StandardisedValue(initialValue - maxDelta);
  }
}
