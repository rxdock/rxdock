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

// Null transform which does nothing except fire off any requests to the
// scoring function

#ifndef _RBTNULLTRANSFORM_H_
#define _RBTNULLTRANSFORM_H_

#include "rxdock/BaseTransform.h"

namespace rxdock {

class NullTransform : public BaseTransform {
public:
  // Static data member for class type
  static const std::string _CT;

  ////////////////////////////////////////
  // Constructors/destructors
  NullTransform(const std::string &strName = "NULL");
  virtual ~NullTransform();

  ////////////////////////////////////////
  // Public methods
  ////////////////

  virtual void Update(Subject *theChangedSubject); // Does nothing

protected:
  ////////////////////////////////////////
  // Protected methods
  ///////////////////
  virtual void Execute(); // Does nothing

private:
  ////////////////////////////////////////
  // Private methods
  /////////////////
  NullTransform(const NullTransform &); // Copy constructor disabled by default
  NullTransform &
  operator=(const NullTransform &); // Copy assignment disabled by default

protected:
  ////////////////////////////////////////
  // Protected data
  ////////////////

private:
  ////////////////////////////////////////
  // Private data
  //////////////
};

// Useful typedefs
typedef SmartPtr<NullTransform> NullTransformPtr; // Smart pointer

} // namespace rxdock

#endif //_RBTNULLTRANSFORM_H_
