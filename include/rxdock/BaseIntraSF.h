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

// Virtual base class for intramolecular scoring functions. Overrides Update()

#ifndef _RBTBASEINTRASF_H_
#define _RBTBASEINTRASF_H_

#include "rxdock/BaseSF.h"
#include "rxdock/Model.h"

namespace rxdock {

class BaseIntraSF : public virtual BaseSF {
public:
  // Class type string
  static const std::string _CT;

  ////////////////////////////////////////
  // Constructors/destructors
  virtual ~BaseIntraSF();

  ////////////////////////////////////////
  // Public methods
  ////////////////

  ModelPtr GetLigand() const;

  // Override Observer pure virtual
  // Notify observer that subject has changed
  virtual void Update(Subject *theChangedSubject);

  // Override BaseSF::ScoreMap to provide additional raw descriptors
  virtual void ScoreMap(StringVariantMap &scoreMap) const;

protected:
  ////////////////////////////////////////
  // Protected methods
  ///////////////////
  BaseIntraSF();

  // PURE VIRTUAL - Derived classes must override
  virtual void SetupScore() = 0; // Called by Update when model has changed

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
  // 26 Mar 2003 (DM) Remember the raw score for the initial ligand conformation
  // This becomes the zero point for all subsequent score reporting
  // i.e. all intramolecular scores are reported relative to the initial score
  double m_zero;
};

} // namespace rxdock

#endif //_RBTBASEINTRASF_H_
