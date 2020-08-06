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

// File reader for external NMR restraints. Supports NOE distance restraints
// and STD ligand-receptor contact surface restraints.
// NOE Syntax:
//<ATOM_LIST_1> <ATOM_LIST_2> <MAX DISTANCE>
//
// STD Syntax:
// STD <ATOM_LIST_1> <MAX DISTANCE> (this is the literal string "STD")
//
// Comment lines are supported (# first character)
// Atom names in list should be comma-separated, with NO SPACES
//
// For both types of restraint, brackets can be used around each atom list to
// denote various levels of ambiguity:
//
//(ATOM_LIST) = apply restraint to mean coords of atoms in list (MEAN)
//[ATOM_LIST] = apply restraint to all atoms in list simultaneously (AND)
// ATOM_LIST = restraint can be satisfied by any atom in list  (OR)
//
// e.g.
// MOL_1:C3 (_22:C5,_23:N1) 5.0 #MOL_1:C3 within 5A of the MEAN of _22:C5,_23:N1
// MOL_1:C3 [_22:C5,_23:N1] 5.0 #MOL_1:C3 within 5A of _22:C5 AND 5A of _23:N1
// MOL_1:C3 _22:C5,_23:N1   5.0 #MOL_1:C3 within 5A of _22:C5 OR 5A of _23:N1
// STD MOL_1:C3 5.0             #MOL_1:C3 within 5A of any part of the receptor

#ifndef _RBTNMRRESTRAINTFILESOURCE_H_
#define _RBTNMRRESTRAINTFILESOURCE_H_

#include "rxdock/BaseFileSource.h"
#include "rxdock/NoeRestraint.h"

namespace rxdock {

class NmrRestraintFileSource : public BaseFileSource {
public:
  ////////////////////////////////////////
  // Constructors/destructors
  NmrRestraintFileSource(const std::string &fileName);

  // Destructor
  virtual ~NmrRestraintFileSource();

  ////////////////////////////////////////
  // Public methods
  ////////////////
  unsigned int GetNumNoeRestraints(); // Number of NOE restraints in file
  NoeRestraintNamesList GetNoeRestraintList(); // List of NOE restraints
  unsigned int GetNumStdRestraints(); // Number of STD restraints in file
  StdRestraintNamesList GetStdRestraintList(); // List of STD restraints

protected:
  ////////////////////////////////////////
  // Protected methods
  ///////////////////

private:
  ////////////////////////////////////////
  // Private methods
  /////////////////

  NmrRestraintFileSource(); // Disable default constructor
  NmrRestraintFileSource(
      const NmrRestraintFileSource &); // Copy constructor disabled by default
  NmrRestraintFileSource &operator=(
      const NmrRestraintFileSource &); // Copy assignment disabled by default

  // Pure virtual in BaseFileSource - needs to be defined here
  virtual void Parse();
  void ClearRestraintCache();

  // Returns NOE restraint type and modifies the atom name string accordingly
  // Returns UNDEFINED if the atom name string has bad syntax
  eNoeType NoeRestraintType(std::string &strAtomNames);

protected:
  ////////////////////////////////////////
  // Protected data
  ////////////////

private:
  ////////////////////////////////////////
  // Private data
  //////////////
  NoeRestraintNamesList m_noeRestraintList;
  StdRestraintNamesList m_stdRestraintList;
};

// Useful typedefs
typedef SmartPtr<NmrRestraintFileSource>
    NmrRestraintFileSourcePtr; // Smart pointer
typedef std::vector<NmrRestraintFileSourcePtr>
    NmrRestraintFileSourceList; // Vector of smart pointers
typedef NmrRestraintFileSourceList::iterator NmrRestraintFileSourceListIter;
typedef NmrRestraintFileSourceList::const_iterator
    NmrRestraintFileSourceListConstIter;

} // namespace rxdock

#endif //_RBTNMRRESTRAINTFILESOURCE_H_
