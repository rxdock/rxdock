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

// Base class for docking site mapper function classes
// Designed in a similar fashion to scoring functions, in that mapper objects
// are attached to a workspace and can respond to changes in the defined
// receptor and ligand models
// Sub classes should implement operator() which returns a list of cavity smart
// pointers

#ifndef _RBTSITEMAPPER_H_
#define _RBTSITEMAPPER_H_

#include "rxdock/BaseObject.h"
#include "rxdock/Cavity.h"
#include "rxdock/Config.h"
#include "rxdock/Model.h"

namespace rxdock {

class SiteMapper : public BaseObject {
public:
  // Class type string
  static std::string _CT;

  ////////////////////////////////////////
  // Constructors/destructors
  virtual ~SiteMapper();
  ////////////////////////////////////////
  // Public methods
  ////////////////
  ModelPtr GetReceptor() const { return m_spReceptor; }

  // PURE VIRTUAL - subclasses must override
  // NB - subclasses should also override Observer::Update pure virtual
  virtual CavityList operator()() = 0;

  // Override Observer pure virtual
  // Notify observer that subject has changed
  virtual void Update(Subject *theChangedSubject);

protected:
  ////////////////////////////////////////
  // Protected methods
  ///////////////////
  SiteMapper(const std::string &strClass, const std::string &strName);

private:
  ////////////////////////////////////////
  // Private methods
  /////////////////
  SiteMapper();                   // default constructor disabled
  SiteMapper(const SiteMapper &); // Copy constructor disabled by default
  SiteMapper &
  operator=(const SiteMapper &); // Copy assignment disabled by default

protected:
  ////////////////////////////////////////
  // Protected data
  ////////////////

private:
  ////////////////////////////////////////
  // Private data
  //////////////
  ModelPtr m_spReceptor;
};

// Useful typedefs
typedef SmartPtr<SiteMapper> SiteMapperPtr; // Smart pointer

} // namespace rxdock

#endif //_RBTSITEMAPPER_H_
