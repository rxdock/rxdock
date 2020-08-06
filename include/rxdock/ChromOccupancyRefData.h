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

// Manages the fixed reference data for a model occupancy chromosome element
// A single instance is designed to be shared between all clones of a given
// element
#ifndef RBTCHROMOCCUPANCYREFDATA_H_
#define RBTCHROMOCCUPANCYREFDATA_H_

#include "rxdock/Config.h"

namespace rxdock {

class Model;

class ChromOccupancyRefData {
public:
  // Class type string
  static std::string _CT;
  // Sole constructor
  ChromOccupancyRefData(
      Model *pModel,
      double stepSize,   // mutation step size
      double threshold); // threshold for enabling/disabling atoms
  virtual ~ChromOccupancyRefData();

  double GetStepSize() const { return m_stepSize; }
  double GetModelValue() const;
  void SetModelValue(double occupancy);
  double GetInitialValue() const { return m_initialValue; }

private:
  Model *m_pModel;
  double m_stepSize;
  double m_threshold;
  double m_initialValue;
};

typedef SmartPtr<ChromOccupancyRefData>
    ChromOccupancyRefDataPtr; // Smart pointer

} // namespace rxdock

#endif /*RBTCHROMOCCUPANCYREFDATA_H_*/
