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

#include "RbtChromDihedralElement.h"

std::string RbtChromDihedralElement::_CT = "RbtChromDihedralElement";

RbtChromDihedralElement::RbtChromDihedralElement(RbtBondPtr spBond,
                                                 RbtAtomList tetheredAtoms,
                                                 double stepSize,
                                                 RbtChromElement::eMode mode,
                                                 double maxDihedral)
    : m_value(0.0) {
  m_spRefData = new RbtChromDihedralRefData(spBond, tetheredAtoms, stepSize,
                                            mode, maxDihedral);
  // Set the initial genotype to match the current phenotype
  SyncFromModel();
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtChromDihedralElement::RbtChromDihedralElement(
    RbtChromDihedralRefDataPtr spRefData, double value)
    : m_spRefData(spRefData), m_value(value) {
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtChromDihedralElement::~RbtChromDihedralElement() {
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

void RbtChromDihedralElement::Reset() {
  m_value = m_spRefData->GetInitialValue();
}

void RbtChromDihedralElement::Randomise() {
  double maxDelta = m_spRefData->GetMaxDihedral();
  double delta;
  switch (m_spRefData->GetMode()) {
  case RbtChromElement::TETHERED:
    delta = 2.0 * maxDelta * GetRand().GetRandom01() - maxDelta;
    m_value = StandardisedValue(m_spRefData->GetInitialValue() + delta);
    break;
  case RbtChromElement::FREE:
    m_value = 360.0 * GetRand().GetRandom01() - 180.0;
    break;
  default:
    m_value = m_spRefData->GetInitialValue();
    break;
  }
}

void RbtChromDihedralElement::Mutate(double relStepSize) {
  double absStepSize = relStepSize * m_spRefData->GetStepSize();
  double delta;
  if (absStepSize > 0) {
    switch (m_spRefData->GetMode()) {
    case RbtChromElement::TETHERED:
      delta = 2.0 * absStepSize * GetRand().GetRandom01() - absStepSize;
      m_value = StandardisedValue(m_value + delta);
      CorrectTetheredDihedral();
      break;
    case RbtChromElement::FREE:
      delta = 2.0 * absStepSize * GetRand().GetRandom01() - absStepSize;
      m_value = StandardisedValue(m_value + delta);
      break;
    default:
      break;
    }
  }
}

void RbtChromDihedralElement::SyncFromModel() {
  m_value = m_spRefData->GetModelValue();
}

void RbtChromDihedralElement::SyncToModel() {
  m_spRefData->SetModelValue(m_value);
}

RbtChromElement *RbtChromDihedralElement::clone() const {
  return new RbtChromDihedralElement(m_spRefData, m_value);
}

void RbtChromDihedralElement::GetVector(RbtDoubleList &v) const {
  v.push_back(m_value);
}

void RbtChromDihedralElement::GetVector(RbtXOverList &v) const {
  RbtXOverElement dihedralElement;
  dihedralElement.push_back(m_value);
  v.push_back(dihedralElement);
}

void RbtChromDihedralElement::SetVector(const RbtDoubleList &v,
                                        int &i) throw(RbtError) {
  if (VectorOK(v, i)) {
    m_value = StandardisedValue(v[i++]);
  } else {
    throw RbtBadArgument(
        _WHERE_, "Index out of range or insufficient elements remaining");
  }
}

void RbtChromDihedralElement::SetVector(const RbtXOverList &v,
                                        int &i) throw(RbtError) {
  if (VectorOK(v, i)) {
    RbtXOverElement dihedralElement(v[i++]);
    if (dihedralElement.size() == 1) {
      m_value = dihedralElement[0];
    } else {
      throw RbtBadArgument(_WHERE_,
                           "dihedralElement vector is of incorrect length");
    }
  } else {
    throw RbtBadArgument(
        _WHERE_, "Index out of range or insufficient elements remaining");
  }
}

void RbtChromDihedralElement::GetStepVector(RbtDoubleList &v) const {
  v.push_back(m_spRefData->GetStepSize());
}

double RbtChromDihedralElement::CompareVector(const RbtDoubleList &v,
                                              int &i) const {
  double retVal(0.0);
  if (!VectorOK(v, i)) {
    retVal = -1.0;
  } else {
    double otherAngle = v[i++];
    double stepSize = m_spRefData->GetStepSize();
    if (stepSize > 0.0) {
      double absDiff = fabs(StandardisedValue(m_value - otherAngle));
      retVal = absDiff / stepSize;
    }
  }
  return retVal;
}

void RbtChromDihedralElement::Print(ostream &s) const {
  s << "DIHEDRAL " << m_value << endl;
}

double RbtChromDihedralElement::StandardisedValue(double dihedralAngle) {
  while (dihedralAngle >= 180) {
    dihedralAngle -= 360.0;
  }
  while (dihedralAngle < -180.0) {
    dihedralAngle += 360.0;
  }
  return dihedralAngle;
}

void RbtChromDihedralElement::CorrectTetheredDihedral() {
  double maxDelta = m_spRefData->GetMaxDihedral();
  double initialValue = m_spRefData->GetInitialValue();
  double delta = StandardisedValue(m_value - initialValue);
  if (delta > maxDelta) {
    m_value = StandardisedValue(initialValue + maxDelta);
  } else if (delta < -maxDelta) {
    m_value = StandardisedValue(initialValue - maxDelta);
  }
}
