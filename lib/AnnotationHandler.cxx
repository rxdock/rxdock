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

#include "rxdock/AnnotationHandler.h"

using namespace rxdock;

// Static data members
const std::string AnnotationHandler::_ANNOTATION_FIELD = "annotations";

////////////////////////////////////////
// Constructors/destructors
AnnotationHandler::AnnotationHandler() : m_bEnabled(false) {
  _RBTOBJECTCOUNTER_CONSTR_("AnnotationHandler");
}

AnnotationHandler::~AnnotationHandler() {
  _RBTOBJECTCOUNTER_DESTR_("AnnotationHandler");
}

////////////////////////////////////////
// Public methods
////////////////
bool AnnotationHandler::isAnnotationEnabled() const { return m_bEnabled; }

const AnnotationList &AnnotationHandler::GetAnnotationList() const {
  return m_annotationList;
}

int AnnotationHandler::GetNumAnnotations() const {
  return m_annotationList.size();
}

void AnnotationHandler::RenderAnnotationList(
    const std::string &strName, std::vector<std::string> &retVal) const {
  // It is callers responsibility to clear retVal before calling
  // RenderAnnotationList Allows for appending to existing string list
  for (AnnotationListConstIter iter = m_annotationList.begin();
       iter != m_annotationList.end(); iter++) {
    retVal.push_back(strName + "," + (*iter)->Render());
  }
}

////////////////////////////////////////
// Protected methods
///////////////////
void AnnotationHandler::AddAnnotation(AnnotationPtr spAnnotation) const {
  if (m_bEnabled) {
    m_annotationList.push_back(spAnnotation);
  }
}

void AnnotationHandler::ClearAnnotationList() const {
  m_annotationList.clear();
}

void AnnotationHandler::EnableAnnotations(bool bEnabled) const {
  m_bEnabled = bEnabled;
}

void rxdock::to_json(json &j, const AnnotationHandler &annotationHandler) {
  json annotationList;
  for (const auto &aIter : annotationHandler.m_annotationList) {
    json annotation = *aIter;
    annotationList.push_back(annotation);
  }
  j = json{{"enabled", annotationHandler.m_bEnabled},
           {"annotations", annotationList}};
}

void rxdock::from_json(const json &j, AnnotationHandler &annotationHandler) {
  j.at("enabled").get_to(annotationHandler.m_bEnabled);
  for (auto &annotation : j.at("annotations")) {
    AnnotationPtr spAnnotation = AnnotationPtr(new Annotation(annotation));
    annotationHandler.m_annotationList.push_back(spAnnotation);
  }
}