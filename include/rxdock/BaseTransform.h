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

// Abstract base class for all classes which manipulate molecular coords

#ifndef _RBTBASETRANSFORM_H_
#define _RBTBASETRANSFORM_H_

#include "rxdock/BaseObject.h"
#include "rxdock/Config.h"

namespace rxdock {

class TransformAgg; // forward declaration

class BaseTransform : public BaseObject {
public:
  // Class type string
  static const std::string _CT;
  // Parameter names

  ////////////////////////////////////////
  // Constructors/destructors
  virtual ~BaseTransform();

  // Give aggregates access to BaseSF private methods/data
  friend class TransformAgg;

  ////////////////////////////////////////
  // Public methods
  ////////////////
  // Fully qualified name, prefixed by all ancestors
  std::string GetFullName() const;

  // Main public method - actually apply the transform
  // Not virtual. Base class method checks if transform is enabled,
  // sends SFRequests to reconfigure the scoring function
  // then calls private virtual method Execute() to apply the transform
  void Go();

  // Aggregate handling methods
  virtual void Add(BaseTransform *);
  virtual void Remove(BaseTransform *);
  virtual bool isAgg() const;
  virtual unsigned int GetNumTransforms() const;
  virtual BaseTransform *GetTransform(unsigned int) const;
  void Orphan(); // Force removal from the parent aggregate
  BaseTransform *GetParentTransform() const;

  // Scoring function request handling
  // Transforms can store up a list of requests to send to the workspace
  // scoring function each time Go() is executed
  void AddSFRequest(RequestPtr);
  void ClearSFRequests();
  void SendSFRequests();

protected:
  ////////////////////////////////////////
  // Protected methods
  ///////////////////
  BaseTransform(const std::string &strClass, const std::string &strName);

private:
  ////////////////////////////////////////
  // Private methods
  /////////////////
  BaseTransform();
  BaseTransform(const BaseTransform &); // Copy constructor disabled by default
  BaseTransform &
  operator=(const BaseTransform &); // Copy assignment disabled by default

  // PURE VIRTUAL - subclasses should override. Applies the transform
  virtual void Execute() = 0;

protected:
  ////////////////////////////////////////
  // Protected data
  ////////////////

private:
  ////////////////////////////////////////
  // Private data
  //////////////
  BaseTransform *m_parent;
  RequestList m_SFRequests;
};

// Useful typedefs
typedef std::vector<BaseTransform *>
    BaseTransformList; // Vector of smart pointers
typedef BaseTransformList::iterator BaseTransformListIter;
typedef BaseTransformList::const_iterator BaseTransformListConstIter;

} // namespace rxdock

#endif //_RBTBASETRANSFORM_H_
