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

#include "RbtChromOccupancyElement.h"

std::string RbtChromOccupancyElement::_CT = "RbtChromOccupancyElement";

RbtChromOccupancyElement::RbtChromOccupancyElement(RbtModel *pModel,
                                                   double stepSize,
                                                   double threshold)
    : m_value(1.0) {
  m_spRefData = new RbtChromOccupancyRefData(pModel, stepSize, threshold);
  // Set the initial genotype to match the current phenotype
  SyncFromModel();
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtChromOccupancyElement::RbtChromOccupancyElement(
    RbtChromOccupancyRefDataPtr spRefData, double value)
    : m_spRefData(spRefData), m_value(value) {
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtChromOccupancyElement::~RbtChromOccupancyElement() {
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

void RbtChromOccupancyElement::Reset() {
  m_value = m_spRefData->GetInitialValue();
}

void RbtChromOccupancyElement::Randomise() {
  m_value = GetRand().GetRandom01();
}

void RbtChromOccupancyElement::Mutate(double relStepSize) {
  double absStepSize = relStepSize * m_spRefData->GetStepSize();
  double delta;
  if (absStepSize > 0) {
    delta = 2.0 * absStepSize * GetRand().GetRandom01() - absStepSize;
    m_value = StandardisedValue(m_value + delta);
  }
}

void RbtChromOccupancyElement::SyncFromModel() {
  m_value = m_spRefData->GetModelValue();
}

void RbtChromOccupancyElement::SyncToModel() {
  m_spRefData->SetModelValue(m_value);
}

RbtChromElement *RbtChromOccupancyElement::clone() const {
  return new RbtChromOccupancyElement(m_spRefData, m_value);
}

void RbtChromOccupancyElement::GetVector(RbtDoubleList &v) const {
  v.push_back(m_value);
}

void RbtChromOccupancyElement::GetVector(RbtXOverList &v) const {
  RbtXOverElement occupancyElement;
  occupancyElement.push_back(m_value);
  v.push_back(occupancyElement);
}

void RbtChromOccupancyElement::SetVector(const RbtDoubleList &v,
                                         int &i) throw(RbtError) {
  if (VectorOK(v, i)) {
    m_value = StandardisedValue(v[i++]);
  } else {
    throw RbtBadArgument(
        _WHERE_, "Index out of range or insufficient elements remaining");
  }
}

void RbtChromOccupancyElement::SetVector(const RbtXOverList &v,
                                         int &i) throw(RbtError) {
  if (VectorOK(v, i)) {
    RbtXOverElement occupancyElement(v[i++]);
    if (occupancyElement.size() == 1) {
      m_value = occupancyElement[0];
    } else {
      throw RbtBadArgument(_WHERE_,
                           "occupancyElement vector is of incorrect length");
    }
  } else {
    throw RbtBadArgument(
        _WHERE_, "Index out of range or insufficient elements remaining");
  }
}

void RbtChromOccupancyElement::GetStepVector(RbtDoubleList &v) const {
  v.push_back(m_spRefData->GetStepSize());
}

double RbtChromOccupancyElement::CompareVector(const RbtDoubleList &v,
                                               int &i) const {
  double retVal(0.0);
  if (!VectorOK(v, i)) {
    retVal = -1.0;
  } else {
    double otherOccupancy = v[i++];
    double stepSize = m_spRefData->GetStepSize();
    if (stepSize > 0.0) {
      double absDiff = fabs(m_value - otherOccupancy);
      retVal = absDiff / stepSize;
    }
  }
  return retVal;
}

void RbtChromOccupancyElement::Print(ostream &s) const {
  s << "OCCUPANCY " << m_value << endl;
}

double RbtChromOccupancyElement::StandardisedValue(double occupancy) {
  if (occupancy < 0.0) {
    occupancy = 0.0;
  } else if (occupancy > 1.0) {
    occupancy = 1.0;
  }
  return occupancy;
}
