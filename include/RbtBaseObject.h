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

#include "RbtConfig.h"
#include "RbtObserver.h"
#include "RbtParamHandler.h"
#include "RbtRequestHandler.h"

class RbtWorkSpace; // forward definition

class RbtBaseObject : public RbtParamHandler,
                      public RbtObserver,
                      public RbtRequestHandler {
public:
  // Class type string
  static std::string _CT;
  // Parameter names
  static std::string _CLASS;
  static std::string _NAME;
  static std::string _ENABLED;
  static std::string _TRACE; // DM 1 Mar 2002 - move from RbtBaseTransform

  ////////////////////////////////////////
  // Constructors/destructors
  virtual ~RbtBaseObject();

  ////////////////////////////////////////
  // Public methods
  ////////////////
  // Class name (e.g. RbtConstSF)
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

  // DM 1 Mar 2002 - move from RbtBaseTransform
  // get/set the trace level for debugging
  int GetTrace() const;
  void SetTrace(int);

  // WorkSpace handling methods
  // Register scoring function with a workspace
  // Base class version just registers itself
  virtual void Register(RbtWorkSpace *);
  // Unregister with a workspace
  // Base class version just unregisters itself
  virtual void Unregister();
  // Get workspace pointer
  RbtWorkSpace *GetWorkSpace() const;

  // Override Observer method
  // Notify observer that subject is about to be deleted
  virtual void Deleted(RbtSubject *theDeletedSubject);

  // Request Handling method
  virtual void HandleRequest(RbtRequestPtr spRequest);

  // Virtual function for dumping parameters to an output stream
  // Called by operator <<
  virtual void Print(std::ostream &s) const;

protected:
  ////////////////////////////////////////
  // Protected methods
  ///////////////////
  RbtBaseObject(const std::string &strClass, const std::string &strName);
  RbtBaseObject();

  // DM 25 Oct 2000 - track changes to parameter values in local data members
  // ParameterUpdated is invoked by RbtParamHandler::SetParameter
  void ParameterUpdated(const std::string &strName);

private:
  ////////////////////////////////////////
  // Private methods
  /////////////////
  RbtBaseObject(const RbtBaseObject &); // Copy constructor disabled by default
  RbtBaseObject &
  operator=(const RbtBaseObject &); // Copy assignment disabled by default

protected:
  ////////////////////////////////////////
  // Protected data
  ////////////////

private:
  ////////////////////////////////////////
  // Private data
  //////////////
  RbtWorkSpace *m_workspace;
  bool m_enabled;
  int m_trace;
};

// Useful typedefs
typedef vector<RbtBaseObject *> RbtBaseObjectList; // Vector of smart pointers
typedef RbtBaseObjectList::iterator RbtBaseObjectListIter;
typedef RbtBaseObjectList::const_iterator RbtBaseObjectListConstIter;

#endif //_RBTBASEOBJECT_H_
