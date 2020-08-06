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

// Base class for unimolecular transforms. Overrides Update()

#ifndef _RBTBASEUNIMOLTRANSFORM_H_
#define _RBTBASEUNIMOLTRANSFORM_H_

#include "BaseTransform.h"
#include "Model.h"

namespace rxdock {

class BaseUniMolTransform : public BaseTransform {
public:
  // Class type string
  static std::string _CT;
  ////////////////////////////////////////
  // Constructors/destructors
  virtual ~BaseUniMolTransform();

  ////////////////////////////////////////
  // Public methods
  ////////////////

  ModelPtr GetLigand() const;

  // Override Observer pure virtual
  // Notify observer that subject has changed
  virtual void Update(Subject *theChangedSubject);

protected:
  ////////////////////////////////////////
  // Protected methods
  ///////////////////
  BaseUniMolTransform(const std::string &strClass, const std::string &strName);

  // PURE VIRTUAL - Derived classes must override
  virtual void SetupTransform() = 0; // Called by Update when model has changed

private:
  ////////////////////////////////////////
  // Private methods
  /////////////////

protected:
  ////////////////////////////////////////
  // Protected data
  ////////////////

private:
  ////////////////////////////////////////
  // Private data
  //////////////
  ModelPtr m_spLigand;
};

} // namespace rxdock

#endif //_RBTBASEUNIMOLTRANSFORM_H_
