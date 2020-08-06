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

// Handles retrieval of atomic (coordinate) info from Charmm CRD files

#ifndef _RBTCRDFILESOURCE_H_
#define _RBTCRDFILESOURCE_H_

#include "BaseMolecularFileSource.h"

namespace rxdock {

class CrdFileSource : public BaseMolecularFileSource {
public:
  // Constructors
  CrdFileSource(const char *fileName);
  CrdFileSource(const std::string fileName);

  // Default destructor
  virtual ~CrdFileSource();

  ////////////////////////////////////////
  // Override public methods from BaseMolecularDataSource
  virtual bool isTitleListSupported() { return true; }
  virtual bool isAtomListSupported() { return true; }
  virtual bool isCoordinatesSupported() { return true; }
  virtual bool isBondListSupported() { return false; }

protected:
  // Pure virtual in BaseFileSource - needs to be defined here
  virtual void Parse();

private:
  // Private methods
  CrdFileSource();                      // Disable default constructor
  CrdFileSource(const CrdFileSource &); // Copy constructor disabled by default
  CrdFileSource &
  operator=(const CrdFileSource &); // Copy assignment disabled by default
};

// useful typedefs
typedef SmartPtr<CrdFileSource> CrdFileSourcePtr; // Smart pointer

} // namespace rxdock

#endif //_RBTCRDFILESOURCE_H_
