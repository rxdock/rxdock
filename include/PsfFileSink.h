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

// Class for writing Models to PSF files

#ifndef _RBTPSFFILESINK_H_
#define _RBTPSFFILESINK_H_

#include "BaseMolecularFileSink.h"

namespace rxdock {

class PsfFileSink : public BaseMolecularFileSink {
public:
  ////////////////////////////////////////
  // Constructors/destructors
  // PsfFileSink(const char* fileName, ModelPtr spModel);
  RBTDLL_EXPORT PsfFileSink(const std::string &fileName, ModelPtr spModel);

  virtual ~PsfFileSink(); // Default destructor

  ////////////////////////////////////////
  // Public methods
  ////////////////
  // Override public methods from BaseFileSink
  virtual void Render();

protected:
  ////////////////////////////////////////
  // Protected methods
  ///////////////////

private:
  ////////////////////////////////////////
  // Private methods
  /////////////////

  PsfFileSink();                    // Disable default constructor
  PsfFileSink(const PsfFileSink &); // Copy constructor disabled by
                                    // default
  PsfFileSink &
  operator=(const PsfFileSink &); // Copy assignment disabled by default

protected:
  ////////////////////////////////////////
  // Protected data
  ////////////////

private:
  ////////////////////////////////////////
  // Private data
  //////////////
};

// Useful typedefs
typedef SmartPtr<PsfFileSink> PsfFileSinkPtr; // Smart pointer

} // namespace rxdock

#endif //_RBTPSFFILESINK_H_
