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

#include "rxdock/ChromOccupancyElement.h"

using namespace rxdock;

std::string ChromOccupancyElement::_CT = "ChromOccupancyElement";

ChromOccupancyElement::ChromOccupancyElement(Model *pModel, double stepSize,
                                             double threshold)
    : m_value(1.0) {
  m_spRefData = new ChromOccupancyRefData(pModel, stepSize, threshold);
  // Set the initial genotype to match the current phenotype
  SyncFromModel();
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

ChromOccupancyElement::ChromOccupancyElement(ChromOccupancyRefDataPtr spRefData,
                                             double value)
    : m_spRefData(spRefData), m_value(value) {
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

ChromOccupancyElement::~ChromOccupancyElement() {
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

void ChromOccupancyElement::Reset() {
  m_value = m_spRefData->GetInitialValue();
}

void ChromOccupancyElement::Randomise() { m_value = GetRand().GetRandom01(); }

void ChromOccupancyElement::Mutate(double relStepSize) {
  double absStepSize = relStepSize * m_spRefData->GetStepSize();
  double delta;
  if (absStepSize > 0) {
    delta = 2.0 * absStepSize * GetRand().GetRandom01() - absStepSize;
    m_value = StandardisedValue(m_value + delta);
  }
}

void ChromOccupancyElement::SyncFromModel() {
  m_value = m_spRefData->GetModelValue();
}

void ChromOccupancyElement::SyncToModel() {
  m_spRefData->SetModelValue(m_value);
}

ChromElement *ChromOccupancyElement::clone() const {
  return new ChromOccupancyElement(m_spRefData, m_value);
}

void ChromOccupancyElement::GetVector(std::vector<double> &v) const {
  v.push_back(m_value);
}

void ChromOccupancyElement::GetVector(XOverList &v) const {
  XOverElement occupancyElement;
  occupancyElement.push_back(m_value);
  v.push_back(occupancyElement);
}

void ChromOccupancyElement::SetVector(const std::vector<double> &v, int &i) {
  if (VectorOK(v, i)) {
    m_value = StandardisedValue(v[i++]);
  } else {
    throw BadArgument(_WHERE_,
                      "Index out of range or insufficient elements remaining");
  }
}

void ChromOccupancyElement::SetVector(const XOverList &v, int &i) {
  if (VectorOK(v, i)) {
    XOverElement occupancyElement(v[i++]);
    if (occupancyElement.size() == 1) {
      m_value = occupancyElement[0];
    } else {
      throw BadArgument(_WHERE_,
                        "occupancyElement vector is of incorrect length");
    }
  } else {
    throw BadArgument(_WHERE_,
                      "Index out of range or insufficient elements remaining");
  }
}

void ChromOccupancyElement::GetStepVector(std::vector<double> &v) const {
  v.push_back(m_spRefData->GetStepSize());
}

double ChromOccupancyElement::CompareVector(const std::vector<double> &v,
                                            int &i) const {
  double retVal(0.0);
  if (!VectorOK(v, i)) {
    retVal = -1.0;
  } else {
    double otherOccupancy = v[i++];
    double stepSize = m_spRefData->GetStepSize();
    if (stepSize > 0.0) {
      double absDiff = std::fabs(m_value - otherOccupancy);
      retVal = absDiff / stepSize;
    }
  }
  return retVal;
}

void ChromOccupancyElement::Print(std::ostream &s) const {
  s << "OCCUPANCY " << m_value << std::endl;
}

double ChromOccupancyElement::StandardisedValue(double occupancy) {
  if (occupancy < 0.0) {
    occupancy = 0.0;
  } else if (occupancy > 1.0) {
    occupancy = 1.0;
  }
  return occupancy;
}
