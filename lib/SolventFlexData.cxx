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

#include "rxdock/SolventFlexData.h"

using namespace rxdock;

const std::string &SolventFlexData::_OCCUPANCY = "occupancy";
const std::string &SolventFlexData::_OCCUPANCY_STEP = "occupancy-step";

SolventFlexData::SolventFlexData(DockingSite *pDockSite)
    : LigandFlexData(pDockSite) {
  AddParameter(_OCCUPANCY, 1.0);
  AddParameter(_OCCUPANCY_STEP, 1.0);
}
