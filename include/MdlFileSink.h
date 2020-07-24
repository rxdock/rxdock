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

// Class for writing Model's to MDL SD and MOL files

#ifndef _RBTMDLFILESINK_H_
#define _RBTMDLFILESINK_H_

#include "BaseMolecularFileSink.h"
#include "ElementFileSource.h"

namespace rxdock {

// DM 19 June 2006
// Map to keep track of logical atom IDs as they are rendered to file
// This avoids the need for the atoms in the Model to have consecutive
// atom IDs. The rendering of enabled solvent models is now supported
// correctly.
// Key = Atom* pointer
// Value = logical atom ID as rendered to file
// Logical atom IDs are consecutive for all atoms rendered
// The actual atom ID (Atom::GetAtomId()) is unchanged
typedef std::map<Atom *, unsigned int, AtomPtrCmp_Ptr> AtomIdMap;

class MdlFileSink : public BaseMolecularFileSink {
public:
  ////////////////////////////////////////
  // Constructors/destructors
  RBTDLL_EXPORT MdlFileSink(const std::string &fileName, ModelPtr spModel);

  virtual ~MdlFileSink(); // Default destructor

  ////////////////////////////////////////
  // Public methods
  ////////////////
  //
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
  void RenderAtomList(const AtomList &atomList);
  void RenderBondList(const BondList &bondList);
  void RenderData(const StringVariantMap &dataMap);

  MdlFileSink();                    // Disable default constructor
  MdlFileSink(const MdlFileSink &); // Copy constructor disabled by
                                    // default
  MdlFileSink &
  operator=(const MdlFileSink &); // Copy assignment disabled by default

protected:
  ////////////////////////////////////////
  // Protected data
  ////////////////

private:
  ////////////////////////////////////////
  // Private data
  //////////////
  ElementFileSourcePtr m_spElementData;
  // DM 27 Apr 1999 - default behaviour is for the first Render to overwrite any
  // existing file then subsequent Renders to append.
  bool m_bFirstRender;
  AtomIdMap m_atomIdMap; // Keep track of logical atom IDs as rendered to
                         // file
};

// Useful typedefs
typedef SmartPtr<MdlFileSink> MdlFileSinkPtr; // Smart pointer

void RenumberScaffold(AtomList &atomList);

} // namespace rxdock

#endif //_RBTMDLFILESINK_H_
