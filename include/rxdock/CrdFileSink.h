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

// Class for writing Models to CRD files

#ifndef _RBTCRDFILESINK_H_
#define _RBTCRDFILESINK_H_

#include "rxdock/BaseMolecularFileSink.h"

#include <nlohmann/json.hpp>

using json = nlohmann::json;

namespace rxdock {

class CrdFileSink : public BaseMolecularFileSink {
public:
  ////////////////////////////////////////
  // Constructors/destructors
  // CrdFileSink(const char* fileName, ModelPtr spModel);
  RBTDLL_EXPORT CrdFileSink(const std::string &fileName, ModelPtr spModel);

  virtual ~CrdFileSink(); // Default destructor

  friend void to_json(json &j, const CrdFileSink &crdFileSink);
  friend void from_json(const json &j, CrdFileSink &crdFileSink);

  ////////////////////////////////////////
  // Public methods
  ////////////////
  //
  // Derived classes should override
  // Controls whether file type can support the writing of multiple
  // conformations/models to a single file
  virtual bool isMultiConfSupported() { return true; }
  //
  ////////////////////////////////////////
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

  CrdFileSink();                    // Disable default constructor
  CrdFileSink(const CrdFileSink &); // Copy constructor disabled by
                                    // default
  CrdFileSink &
  operator=(const CrdFileSink &); // Copy assignment disabled by default

protected:
  ////////////////////////////////////////
  // Protected data
  ////////////////

private:
  ////////////////////////////////////////
  // Private data
  //////////////
  unsigned int m_numAtomsLineRec;
};

void to_json(json &j, const CrdFileSink &crdFileSink);
void from_json(const json &j, CrdFileSink &crdFileSink);

// Useful typedefs
typedef SmartPtr<CrdFileSink> CrdFileSinkPtr; // Smart pointer

} // namespace rxdock

#endif //_RBTCRDFILESINK_H_
