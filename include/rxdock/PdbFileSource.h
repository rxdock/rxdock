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

// Handles retrieval of atomic (coordinate) info
// from XPLOR PDB files (with segment name in last column).

#ifndef _RBTPDBFILESOURCE_H_
#define _RBTPDBFILESOURCE_H_

#include "rxdock/BaseMolecularFileSource.h"

namespace rxdock {

class PdbFileSource : public BaseMolecularFileSource {
public:
  // Constructors
  PdbFileSource(const char *fileName);
  PdbFileSource(const std::string fileName);

  // Default destructor
  virtual ~PdbFileSource();

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
  PdbFileSource();                      // Disable default constructor
  PdbFileSource(const PdbFileSource &); // Copy constructor disabled by default
  PdbFileSource &
  operator=(const PdbFileSource &); // Copy assignment disabled by default
};

// useful typedefs
typedef SmartPtr<PdbFileSource> PdbFileSourcePtr; // Smart pointer

} // namespace rxdock

#endif //_RBTPDBFILESOURCE_H_
