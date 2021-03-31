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

#include "rxdock/ChromOccupancyRefData.h"
#include "rxdock/Model.h"

using namespace rxdock;

const std::string ChromOccupancyRefData::_CT = "ChromOccupancyRefData";

ChromOccupancyRefData::ChromOccupancyRefData(Model *pModel, double stepSize,
                                             double threshold)
    : m_pModel(pModel), m_stepSize(stepSize), m_threshold(threshold) {
  m_initialValue = GetModelValue();
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

ChromOccupancyRefData::~ChromOccupancyRefData() {
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

double ChromOccupancyRefData::GetModelValue() const {
  return (m_pModel) ? m_pModel->GetOccupancy() : 0.0;
}

void ChromOccupancyRefData::SetModelValue(double occupancy) {
  if (m_pModel) {
    m_pModel->SetOccupancy(occupancy, m_threshold);
  }
}
