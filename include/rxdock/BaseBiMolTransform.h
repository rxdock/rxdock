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

// Base class for bimolecular transforms. Overrides Update()

#ifndef _RBTBASEBIMOLTRANSFORM_H_
#define _RBTBASEBIMOLTRANSFORM_H_

#include "rxdock/BaseTransform.h"
#include "rxdock/Model.h"

namespace rxdock {

class BaseBiMolTransform : public BaseTransform {
public:
  // Class type string
  static std::string _CT;
  ////////////////////////////////////////
  // Constructors/destructors
  virtual ~BaseBiMolTransform();

  ////////////////////////////////////////
  // Public methods
  ////////////////

  ModelPtr GetReceptor() const;
  ModelPtr GetLigand() const;
  ModelList GetSolvent() const;

  // Override Observer pure virtual
  // Notify observer that subject has changed
  virtual void Update(Subject *theChangedSubject);

protected:
  ////////////////////////////////////////
  // Protected methods
  ///////////////////
  BaseBiMolTransform(const std::string &strClass, const std::string &strName);

  // PURE VIRTUAL - Derived classes must override
  virtual void SetupReceptor() = 0; // Called by Update when receptor is changed
  virtual void SetupLigand() = 0;   // Called by Update when ligand is changed
  virtual void SetupSolvent() {}    // Called by Update when ligand is changed
  virtual void
  SetupTransform() = 0; // Called by Update when either model has changed

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
  ModelPtr m_spReceptor;
  ModelPtr m_spLigand;
  ModelList m_solventList;
};

} // namespace rxdock

#endif //_RBTBASEBIMOLTRANSFORM_H_
