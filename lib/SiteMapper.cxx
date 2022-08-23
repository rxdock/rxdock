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

#include "rxdock/SiteMapper.h"
#include "rxdock/WorkSpace.h"

#include <loguru.hpp>

using namespace rxdock;

// Static data members
const std::string SiteMapper::_CT = "SiteMapper";

////////////////////////////////////////
// Constructors/destructors
SiteMapper::SiteMapper(const std::string &strClass, const std::string &strName)
    : BaseObject(strClass, strName) {
  LOG_F(2, "SiteMapper parameterised constructor");
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

SiteMapper::~SiteMapper() {
  LOG_F(2, "SiteMapper destructor");
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

// Override Observer pure virtual
// Notify observer that subject has changed
void SiteMapper::Update(Subject *theChangedSubject) {
  WorkSpace *pWorkSpace = GetWorkSpace();
  if (theChangedSubject == pWorkSpace) {
    // Check if receptor has been updated (model #0)
    if (pWorkSpace->GetNumModels() >= 1) {
      ModelPtr spReceptor = GetWorkSpace()->GetModel(0);
      if (spReceptor != m_spReceptor) {
        m_spReceptor = spReceptor;
        LOG_F(2, "SiteMapper::Update(): Receptor has been updated");
      }
    }
  }
}

void rxdock::to_json(json &j, const SiteMapper &siteMapper) {
  j = json{{"sp-receptor", *siteMapper.m_spReceptor}};
}

void rxdock::from_json(const json &j, SiteMapper &siteMapper) {
  j.at("sp-receptor").get_to(*siteMapper.m_spReceptor);
}