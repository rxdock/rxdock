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

// Randomises all rotable bonds of the ligand, within specified +/- range
// No bump check on final conformation

#ifndef _RBTRANDLIGTRANSFORM_H_
#define _RBTRANDLIGTRANSFORM_H_

#include "rxdock/BaseUniMolTransform.h"
#include "rxdock/Rand.h"

#include <nlohmann/json.hpp>

using json = nlohmann::json;

namespace rxdock {

class RandLigTransform : public BaseUniMolTransform {
public:
  // Static data member for class type
  static const std::string _CT;
  // Parameter names
  static const std::string _TORS_STEP;

  ////////////////////////////////////////
  // Constructors/destructors
  RandLigTransform(const std::string &strName = "RANDLIG");
  virtual ~RandLigTransform();

  friend void to_json(json &j, const RandLigTransform &randLigTrans);
  friend void from_json(const json &j, RandLigTransform &randLigTrans);

  ////////////////////////////////////////
  // Public methods
  ////////////////

protected:
  ////////////////////////////////////////
  // Protected methods
  ///////////////////
  virtual void SetupTransform(); // Called by Update when model has changed
  virtual void Execute();

private:
  ////////////////////////////////////////
  // Private methods
  /////////////////
  RandLigTransform(
      const RandLigTransform &); // Copy constructor disabled by default
  RandLigTransform &
  operator=(const RandLigTransform &); // Copy assignment disabled by default

protected:
  ////////////////////////////////////////
  // Protected data
  ////////////////

private:
  ////////////////////////////////////////
  // Private data
  //////////////
  Rand &m_rand; // keep a reference to the singleton random number generator
  BondList m_rotableBonds;
};

void to_json(json &j, const RandLigTransform &randLigTrans);
void from_json(const json &j, RandLigTransform &randLigTrans);

// Useful typedefs
typedef SmartPtr<RandLigTransform> RandLigTransformPtr; // Smart pointer

} // namespace rxdock

#endif //_RBTRANDLIGTRANSFORM_H_
