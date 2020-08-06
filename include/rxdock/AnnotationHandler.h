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

// Base implementation class for managing an annotation list. Scoring function
// classes wishing to store annotations for later retrieval should derive from
// this class.

#ifndef _RBTANNOTATIONHANDLER_H_
#define _RBTANNOTATIONHANDLER_H_

#include "rxdock/Annotation.h"

namespace rxdock {

class AnnotationHandler {
public:
  static std::string _ANNOTATION_FIELD;

  ////////////////////////////////////////
  // Constructors/destructors
  virtual ~AnnotationHandler(); // Default destructor

  ////////////////////////////////////////
  // Public methods
  ////////////////
  bool isAnnotationEnabled() const;
  // Get a const ref to the annotation list (read-only)
  const AnnotationList &GetAnnotationList() const;
  int GetNumAnnotations() const;
  void RenderAnnotationList(const std::string &strName,
                            std::vector<std::string> &retVal) const;

protected:
  ////////////////////////////////////////
  // Protected methods
  ///////////////////
  RBTDLL_EXPORT AnnotationHandler();
  void AddAnnotation(AnnotationPtr spAnnotation) const;
  void ClearAnnotationList() const;
  void EnableAnnotations(bool bEnabled) const;

private:
  ////////////////////////////////////////
  // Private methods
  /////////////////
  AnnotationHandler(
      const AnnotationHandler &); // Copy constructor disabled by default
  AnnotationHandler &
  operator=(const AnnotationHandler &); // Copy assignment disabled by default

protected:
  ////////////////////////////////////////
  // Protected data
  ////////////////

private:
  ////////////////////////////////////////
  // Private data
  //////////////
  mutable bool m_bEnabled;
  mutable AnnotationList m_annotationList;
};

} // namespace rxdock

#endif //_RBTANNOTATIONHANDLER_H_
