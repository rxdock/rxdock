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

// Transform aggregate class. Manages collection of child transforms

#ifndef _RBTTRANSFORMAGG_H_
#define _RBTTRANSFORMAGG_H_

#include "BaseTransform.h"

namespace rxdock {

// Only check transform aggregate assertions in debug build
#ifdef _NDEBUG
const Bool TRANSFORMAGG_CHECK = false;
#else
const bool TRANSFORMAGG_CHECK = true;
#endif //_NDEBUG

// Useful typedefs

class TransformAgg : public BaseTransform {
public:
  // Static data member for class type (i.e. "TransformAgg")
  static std::string _CT;

  ////////////////////////////////////////
  // Constructors/destructors
  RBTDLL_EXPORT TransformAgg(const std::string &strName = "DOCK");
  virtual ~TransformAgg();

  ////////////////////////////////////////
  // Public methods
  ////////////////

  // Aggregate handling methods
  virtual void Add(BaseTransform *);
  virtual void Remove(BaseTransform *);
  virtual bool isAgg() const;
  virtual unsigned int GetNumTransforms() const;
  virtual BaseTransform *GetTransform(unsigned int iTransform) const;

  // WorkSpace handling methods
  // Register scoring function with a workspace
  // Aggregate version registers all children, but NOT itself
  //(Aggregates are just containers, and have no need for model information
  virtual void Register(WorkSpace *);
  // Unregister with a workspace
  // Aggregate version unregisters all children, but NOT itself
  virtual void Unregister();

  // Override Observer pure virtual
  // Notify observer that subject has changed
  // Does nothing in TransformAgg as aggregates do not require updating
  virtual void Update(Subject *theChangedSubject);

  // Request Handling method
  virtual void HandleRequest(RequestPtr spRequest);

  // Virtual function for dumping transform details to an output stream
  // Called by operator <<
  virtual void Print(std::ostream &s) const;

protected:
  ////////////////////////////////////////
  // Protected methods
  ///////////////////
  // Actually apply the transform
  // Aggregate version loops over all child transforms
  virtual void Execute();

private:
  ////////////////////////////////////////
  // Private methods
  /////////////////
  TransformAgg(const TransformAgg &); // Copy constructor disabled by default
  TransformAgg &
  operator=(const TransformAgg &); // Copy assignment disabled by default

protected:
  ////////////////////////////////////////
  // Protected data
  ////////////////

private:
  ////////////////////////////////////////
  // Private data
  //////////////
  BaseTransformList m_transforms;
};

// Useful typedefs
typedef SmartPtr<TransformAgg> TransformAggPtr; // Smart pointer
typedef std::vector<TransformAggPtr>
    TransformAggList; // Vector of smart pointers
typedef TransformAggList::iterator TransformAggListIter;
typedef TransformAggList::const_iterator TransformAggListConstIter;

} // namespace rxdock

#endif //_RBTTRANSFORMAGG_H_
