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

#include "LigandFlexData.h"
#include "ChromElement.h"

using namespace rxdock;

const std::string &LigandFlexData::_TRANS_STEP = "TRANS_STEP";
const std::string &LigandFlexData::_ROT_STEP = "ROT_STEP";
const std::string &LigandFlexData::_DIHEDRAL_STEP = "DIHEDRAL_STEP";
const std::string &LigandFlexData::_TRANS_MODE = "TRANS_MODE";
const std::string &LigandFlexData::_ROT_MODE = "ROT_MODE";
const std::string &LigandFlexData::_DIHEDRAL_MODE = "DIHEDRAL_MODE";
const std::string &LigandFlexData::_MAX_TRANS = "MAX_TRANS";
const std::string &LigandFlexData::_MAX_ROT = "MAX_ROT";
const std::string &LigandFlexData::_MAX_DIHEDRAL = "MAX_DIHEDRAL";

const std::string &LigandFlexData::GetTransMode() { return _TRANS_MODE; }

const std::string &LigandFlexData::GetRotMode() { return _ROT_MODE; }

LigandFlexData::LigandFlexData(DockingSite *pDockSite) : FlexData(pDockSite) {
  AddParameter(_TRANS_STEP, 2.0);
  AddParameter(_ROT_STEP, 30.0);
  AddParameter(_DIHEDRAL_STEP, 30.0);
  AddParameter(_TRANS_MODE, ChromElement::ModeToStr(ChromElement::FREE));
  AddParameter(_ROT_MODE, ChromElement::ModeToStr(ChromElement::FREE));
  AddParameter(_DIHEDRAL_MODE, ChromElement::ModeToStr(ChromElement::FREE));
  AddParameter(_MAX_TRANS, 1.0);
  AddParameter(_MAX_ROT, 30.0);
  AddParameter(_MAX_DIHEDRAL, 30.0);
}
