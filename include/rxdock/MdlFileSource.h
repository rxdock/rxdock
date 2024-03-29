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

// Handles retrieval of molecular info from MDL Mol and SD files.

#ifndef _RBTMDLFILESOURCE_H_
#define _RBTMDLFILESOURCE_H_

#include "rxdock/BaseMolecularFileSource.h"
#include "rxdock/ElementFileSource.h"

#include <nlohmann/json.hpp>

using json = nlohmann::json;

namespace rxdock {

const std::string IDS_MDL_RECDELIM = "$$$$";

class MdlFileSource : public BaseMolecularFileSource {
public:
  // Constructors
  // MdlFileSource(const char* fileName);
  RBTDLL_EXPORT MdlFileSource(const std::string &fileName,
                              bool bPosIonisable = true,
                              bool bNegIonisable = true,
                              bool bImplHydrogens = true);

  // Default destructor
  virtual ~MdlFileSource();

  friend void to_json(json &j, const MdlFileSource &mdlFileSrc);
  friend void from_json(const json &j, MdlFileSource &mdlFileSrc);

  ////////////////////////////////////////
  // Override public methods from BaseMolecularDataSource
  virtual bool isTitleListSupported() { return true; }
  virtual bool isAtomListSupported() { return true; }
  virtual bool isCoordinatesSupported() { return true; }
  virtual bool isBondListSupported() { return true; }
  // DM 12 May 1999 - support for data records
  virtual bool isDataSupported() { return true; }

protected:
  // Pure virtual in BaseFileSource - needs to be defined here
  virtual void Parse();

private:
  // Private methods
  MdlFileSource();                      // Disable default constructor
  MdlFileSource(const MdlFileSource &); // Copy constructor disabled by default
  MdlFileSource &
  operator=(const MdlFileSource &); // Copy assignment disabled by default

  // Sets up all the atomic attributes that are not explicitly stored in the MDL
  // file
  void SetupAtomParams();

  // Helper functions for SetupAtomParams()

  // Sets hybridisation state and checks for Nsp3 next to sp2 special case
  void SetupHybridState();

  // Checks for valency within range
  // Adds implicit hydrogens to C,N,S
  // Defines vdW radius, correcting for extended atoms and H-bond donor
  // hydrogens Defines formal "force-field" type string
  void SetupTheRest();

  void SetupPosIonisableGroups();
  void SetupNegIonisableGroups();

  void AddHydrogen(AtomPtr spAtom);

  // SetupIonicGroups sets the "group charge" attribute of the atoms. The idea
  // is that the "group charge" is file-format independent, and eliminates the
  // need for the scoring function setup to have any knowledge of the input file
  // format. For example, PSF files represent COO- as OC-C-OC (both oxygens
  // charged) whereas SD files represent it as the formal resonance structure
  // O=C-O-. For rDock we actually want the negative charge on the central
  // carbon!
  //
  // It is envisaged each molecular file source will have a version of
  // SetupIonicGroups to convert from the native representation to the rDock
  // representation. The FormalCharge and PartialCharge attributes should be
  // left untouched to allow the model to be rewritten back in the same format.
  void SetupIonicGroups();
  void SetupNSP3Plus(); // Helper function to set up protonated amines
  void
  SetupNSP2Plus(); // Helper function to set up guanidines, imidazoles, amidines
  void SetupOTRIMinus(); // Helper function to set up deprotonated acids

  // DM 4 June 1999 - Search for the separate fragments in the record
  // and set distinct segment names for each one
  void SetupSegmentNames();

  // DM 2 Aug 1999 - remove all non-polar hydrogens
  void RemoveNonPolarHydrogens();

  // DM 8 Feb 2000 - setup atom and bond cyclic flags (previously in Model)
  // Set the atom and bond cyclic flags for all atoms and bonds in the source
  void SetCyclicFlags();

  // Private data
  ElementFileSourcePtr
      m_spElementData;   // Elemental data source for vdW radii etc
  bool m_bPosIonisable;  // If true, protonate Nsp3, and Nsp2 in imidazole and
                         // guanidinium
  bool m_bNegIonisable;  // If true, deprotonate Otri-Hsp3
  bool m_bImplHydrogens; // If true, non-hydrogen bonding hydrogens are removed
};

void to_json(json &j, const MdlFileSource &mdlFileSrc);
void from_json(const json &j, MdlFileSource &mdlFileSrc);

// useful typedefs
typedef SmartPtr<MdlFileSource> MdlFileSourcePtr; // Smart pointer

} // namespace rxdock

#endif //_RBTMDLFILESOURCE_H_
