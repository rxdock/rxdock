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

// Generates a randomised population for the GA
#ifndef _RBTRANDPOPTRANSFORM_H_
#define _RBTRANDPOPTRANSFORM_H_

#include "BaseBiMolTransform.h"
#include "ChromElement.h"

namespace rxdock {

class RandPopTransform : public BaseBiMolTransform {
public:
  static std::string _CT;
  static std::string _POP_SIZE;
  static std::string _SCALE_CHROM_LENGTH;

  ////////////////////////////////////////
  // Constructors/destructors
  RBTDLL_EXPORT RandPopTransform(const std::string &strName = "RANDPOP");
  virtual ~RandPopTransform();

  ////////////////////////////////////////
  // Public methods
  ////////////////

protected:
  ////////////////////////////////////////
  // Protected methods
  ///////////////////
  virtual void
  SetupTransform(); // Called by Update when either model has changed
  virtual void SetupReceptor(); // Called by Update when receptor is changed
  virtual void SetupLigand();   // Called by Update when ligand is changed
  virtual void SetupSolvent();  // Called by Update when solvent is changed
  virtual void Execute();

private:
  ////////////////////////////////////////
  // Private methods
  /////////////////
  RandPopTransform(
      const RandPopTransform &); // Copy constructor disabled by default
  RandPopTransform &
  operator=(const RandPopTransform &); // Copy assignment disabled by default

protected:
  ////////////////////////////////////////
  // Protected data
  ////////////////

private:
  ////////////////////////////////////////
  // Private data
  //////////////
  ChromElementPtr m_chrom;
};

// Useful typedefs
typedef SmartPtr<RandPopTransform> RandPopTransformPtr; // Smart pointer

} // namespace rxdock

#endif //_RBTRANDPOPTRANSFORM_H_
