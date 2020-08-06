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

// Base class for all workspace observers which have named parameters and
// handle arbitrary requests (i.e. scoring functions and transforms)

#ifndef _RBTBASEOBJECT_H_
#define _RBTBASEOBJECT_H_

#include "rxdock/Config.h"
#include "rxdock/Observer.h"
#include "rxdock/ParamHandler.h"
#include "rxdock/RequestHandler.h"

namespace rxdock {

class WorkSpace; // forward definition

class BaseObject : public ParamHandler, public Observer, public RequestHandler {
public:
  // Class type string
  static std::string _CT;
  // Parameter names
  static std::string _CLASS;
  static std::string _NAME;
  static std::string _ENABLED;

  ////////////////////////////////////////
  // Constructors/destructors
  virtual ~BaseObject();

  ////////////////////////////////////////
  // Public methods
  ////////////////
  // Class name (e.g. ConstSF)
  std::string GetClass() const;
  // Short name (e.g. HBOND)
  std::string GetName() const;
  void SetName(const std::string &);
  // Fully qualified name (should be overridden by subclasses which can be
  // aggregated to prefix the name with the parent's name)
  virtual std::string GetFullName() const;
  void Enable();
  void Disable();
  bool isEnabled() const;

  // WorkSpace handling methods
  // Register scoring function with a workspace
  // Base class version just registers itself
  virtual void Register(WorkSpace *);
  // Unregister with a workspace
  // Base class version just unregisters itself
  virtual void Unregister();
  // Get workspace pointer
  WorkSpace *GetWorkSpace() const;

  // Override Observer method
  // Notify observer that subject is about to be deleted
  virtual void Deleted(Subject *theDeletedSubject);

  // Request Handling method
  virtual void HandleRequest(RequestPtr spRequest);

  // Virtual function for dumping parameters to an output stream
  // Called by operator <<
  virtual void Print(std::ostream &s) const;

protected:
  ////////////////////////////////////////
  // Protected methods
  ///////////////////
  BaseObject(const std::string &strClass, const std::string &strName);
  BaseObject();

  // DM 25 Oct 2000 - track changes to parameter values in local data members
  // ParameterUpdated is invoked by ParamHandler::SetParameter
  void ParameterUpdated(const std::string &strName);

private:
  ////////////////////////////////////////
  // Private methods
  /////////////////
  BaseObject(const BaseObject &); // Copy constructor disabled by default
  BaseObject &
  operator=(const BaseObject &); // Copy assignment disabled by default

protected:
  ////////////////////////////////////////
  // Protected data
  ////////////////

private:
  ////////////////////////////////////////
  // Private data
  //////////////
  WorkSpace *m_workspace;
  bool m_enabled;
};

// Useful typedefs
typedef std::vector<BaseObject *> BaseObjectList; // Vector of smart pointers
typedef BaseObjectList::iterator BaseObjectListIter;
typedef BaseObjectList::const_iterator BaseObjectListConstIter;

} // namespace rxdock

#endif //_RBTBASEOBJECT_H_
