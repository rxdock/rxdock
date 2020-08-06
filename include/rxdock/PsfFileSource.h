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

// Handles retrieval of molecular info from Charmm PSF files.

#ifndef _RBTPSFFILESOURCE_H_
#define _RBTPSFFILESOURCE_H_

#include "rxdock/BaseMolecularFileSource.h"
#include "rxdock/CharmmDataSource.h"
#include "rxdock/ElementFileSource.h"
#include "rxdock/ParameterFileSource.h"

namespace rxdock {

class PsfFileSource : public BaseMolecularFileSource {
public:
  // Constructors
  PsfFileSource(
      const std::string &fileName,
      const std::string &strMassesFile = GetDataFileName("data", "masses.rtf"),
      bool bImplHydrogens = true);

  // Default destructor
  virtual ~PsfFileSource();

  ////////////////////////////////////////
  // Override public methods from BaseMolecularDataSource
  virtual bool isTitleListSupported() { return true; }
  virtual bool isAtomListSupported() { return true; }
  virtual bool isCoordinatesSupported() { return false; }
  virtual bool isBondListSupported() { return true; }

protected:
  // Pure virtual in BaseFileSource - needs to be defined here
  virtual void Parse();

private:
  // Private methods
  PsfFileSource();                      // Disable default constructor
  PsfFileSource(const PsfFileSource &); // Copy constructor disabled by default
  PsfFileSource &
  operator=(const PsfFileSource &); // Copy assignment disabled by default

  // Sets up all the atomic attributes that are not explicitly stored in the PSF
  // file
  void SetupAtomParams();
  void SetupVdWRadii();
  void SetupPartialIonicGroups();
  void RemoveNonPolarHydrogens();
  // Is atom in same substructure
  class isSS_eq : public std::function<bool(Atom *)> {
    const Atom *a1;

  public:
    explicit isSS_eq(const Atom *aa) : a1(aa) {}
    bool operator()(const Atom *a2) const {
      return (a1->GetSubunitName() == a2->GetSubunitName()) &&
             (a1->GetSubunitId() == a2->GetSubunitId()) &&
             (a1->GetSegmentName() == a2->GetSegmentName());
    }
  };

  // Private data
  ElementFileSourcePtr
      m_spElementData;                // Elemental data source for vdW radii etc
  CharmmDataSourcePtr m_spCharmmData; // Charmm data source for fftype lookup
  ParameterFileSourcePtr
      m_spParamSource;   // Parameter file source for partial charges etc
  bool m_bImplHydrogens; // If true, non-hydrogen bonding hydrogens are removed
};

// useful typedefs
typedef SmartPtr<PsfFileSource> PsfFileSourcePtr; // Smart pointer

} // namespace rxdock

#endif //_RBTPSFFILESOURCE_H_
