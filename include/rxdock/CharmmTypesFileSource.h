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

// File reader for Charmm MASSES.RTF atom type file
// Provides the atom type list as a vector of structs

#ifndef _RBTCHARMMTYPESFILESOURCE_H_
#define _RBTCHARMMTYPESFILESOURCE_H_

#include "BaseFileSource.h"

namespace rxdock {

// Simple struct for holding the type info in masses.rtf
class CharmmType {
public:
  CharmmType()
      : nAtomType(0), strAtomType(""), mass(0.0), element(""), comment("") {}
  int nAtomType;
  std::string strAtomType;
  double mass;
  std::string element;
  std::string comment;
};

typedef std::vector<CharmmType> CharmmTypeList;
typedef CharmmTypeList::iterator CharmmTypeListIter;

class CharmmTypesFileSource : public BaseFileSource {
public:
  // Constructors
  CharmmTypesFileSource(const char *fileName);
  CharmmTypesFileSource(const std::string fileName);

  // Destructor
  virtual ~CharmmTypesFileSource();

  ////////////////////////////////////////
  // Public methods
  int GetNumTypes();
  CharmmTypeList GetTypeList();

protected:
  // Protected methods

private:
  // Private methods
  CharmmTypesFileSource(); // Disable default constructor
  CharmmTypesFileSource(
      const CharmmTypesFileSource &); // Copy constructor disabled by default
  CharmmTypesFileSource &operator=(
      const CharmmTypesFileSource &); // Copy assignment disabled by default

  // Pure virtual in BaseFileSource - needs to be defined here
  virtual void Parse();
  void ClearTypesCache();

protected:
  // Protected data

private:
  // Private data
  CharmmTypeList m_typesList;
};

// useful typedefs
typedef SmartPtr<CharmmTypesFileSource>
    CharmmTypesFileSourcePtr; // Smart pointer

} // namespace rxdock

#endif //_RBTCHARMMTYPESFILESOURCE_H_
