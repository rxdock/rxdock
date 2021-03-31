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

// Simplex Minimiser
#ifndef _RBTSIMPLEXTRANSFORM_H_
#define _RBTSIMPLEXTRANSFORM_H_

#include "rxdock/BaseBiMolTransform.h"
#include "rxdock/ChromElement.h"

namespace rxdock {

class SimplexTransform : public BaseBiMolTransform {
public:
  // Static data member for class type
  static const std::string _CT;
  // Parameter names
  static const std::string _MAX_CALLS;
  static const std::string _NCYCLES;
  static const std::string _STOPPING_STEP_LENGTH;
  static const std::string _PARTITION_DIST;
  static const std::string _STEP_SIZE;
  // Stop once score improves by less than convergence value
  // between cycles
  static const std::string _CONVERGENCE;

  RBTDLL_EXPORT static const std::string &GetMaxCalls();
  RBTDLL_EXPORT static const std::string &GetNCycles();
  RBTDLL_EXPORT static const std::string &GetStepSize();

  ////////////////////////////////////////
  // Constructors/destructors
  RBTDLL_EXPORT SimplexTransform(const std::string &strName = "SIMPLEX");
  virtual ~SimplexTransform();

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
  SimplexTransform(
      const SimplexTransform &); // Copy constructor disabled by default
  SimplexTransform &
  operator=(const SimplexTransform &); // Copy assignment disabled by default
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
typedef SmartPtr<SimplexTransform> SimplexTransformPtr; // Smart pointer

} // namespace rxdock

#endif //_RBTSIMPLEXTRANSFORM_H_
