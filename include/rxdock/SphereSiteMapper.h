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

// 2-sphere site mapper function class

#ifndef _RBTSPHERESITEMAPPER_H_
#define _RBTSPHERESITEMAPPER_H_

#include "rxdock/SiteMapper.h"

namespace rxdock {

class SphereSiteMapper : public SiteMapper {
public:
  // Static data member for class type
  static std::string _CT;
  // Parameter names
  static std::string _VOL_INCR;
  static std::string _SMALL_SPHERE;
  static std::string _LARGE_SPHERE;
  static std::string _GRIDSTEP;
  static std::string _CENTER;
  static std::string _RADIUS;
  static std::string _MIN_VOLUME;
  static std::string _MAX_CAVITIES;

  SphereSiteMapper(const std::string &strName = "SPHERE_MAPPER");
  virtual ~SphereSiteMapper();

  // Override SiteMapper pure virtual
  // This is the function which actually does the mapping
  virtual CavityList operator()();
};

// Useful typedefs
typedef SmartPtr<SphereSiteMapper> SphereSiteMapperPtr; // Smart pointer

} // namespace rxdock

#endif //_RBTSPHERESITEMAPPER_H_
