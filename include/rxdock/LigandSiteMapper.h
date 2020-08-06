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

// Reference ligand site mapper function class

#ifndef _RBTLIGANDSITEMAPPER_H_
#define _RBTLIGANDSITEMAPPER_H_

#include "SiteMapper.h"

namespace rxdock {

class LigandSiteMapper : public SiteMapper {
public:
  // Static data member for class type
  static std::string _CT;
  // Parameter names
  static std::string _REF_MOL;
  static std::string _VOL_INCR;
  static std::string _SMALL_SPHERE;
  static std::string _GRIDSTEP;
  static std::string _RADIUS;
  static std::string _MIN_VOLUME;
  static std::string _MAX_CAVITIES;

  LigandSiteMapper(const std::string &strName = "LIGAND_MAPPER");
  virtual ~LigandSiteMapper();

  // Override SiteMapper pure virtual
  // This is the function which actually does the mapping
  virtual CavityList operator()();
};

// Useful typedefs
typedef SmartPtr<LigandSiteMapper> LigandSiteMapperPtr; // Smart pointer

} // namespace rxdock

#endif //_RBTLIGANDSITEMAPPER_H_
