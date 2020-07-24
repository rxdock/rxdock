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

// Legacy class, largely replaced by ChromElement subclasses
// The only remaining purpose of ModelMutator is to generate lists
// of flexible intramolecular interactions.
// i.e. intramolecular atom pairs across rotatable bonds
#ifndef _RBTMODELMUTATOR_H_
#define _RBTMODELMUTATOR_H_

#include "Atom.h"
#include "Bond.h"
#include "PrincipalAxes.h"

namespace rxdock {

class Model; // forward declaration

// Only check request assertions in debug build
#ifdef _NDEBUG
const Bool MUT_CHECK = true;
#else
const bool MUT_CHECK = true;
#endif //_NDEBUG

class ModelMutator {
public:
  ////////////////////////////////////////
  // Constructors/destructors
  ModelMutator(Model *pModel, const BondList &rotBonds,
               const AtomList &tetheredAtoms);
  virtual ~ModelMutator(); // Default destructor

  ////////////////////////////////////////
  // Public methods
  ////////////////
  const AtomRListList &GetFlexIntns() const;
  const AtomRListList &GetFlexAtoms() const;
  BondList GetFlexBonds() const;

protected:
  ////////////////////////////////////////
  // Protected methods
  ///////////////////

private:
  ////////////////////////////////////////
  // Private methods
  /////////////////

  ModelMutator(const ModelMutator &); // Copy constructor disabled by default
  ModelMutator &
  operator=(const ModelMutator &); // Copy assignment disabled by default

  void Setup();

protected:
  ////////////////////////////////////////
  // Protected data
  ////////////////

private:
  ////////////////////////////////////////
  // Private data
  //////////////
  Model *m_pModel;       // Parent model
  BondList m_rotBonds;   // A record of the defined rotatable bonds passed to
                         // the constructor
  AtomRList m_dih1Atoms; // Atom 1 in each rotable bond dihedral definition
  AtomRList m_dih2Atoms; // Atom 2 in each rotable bond dihedral definition
  AtomRList m_dih3Atoms; // Atom 3 in each rotable bond dihedral definition
  AtomRList m_dih4Atoms; // Atom 4 in each rotable bond dihedral definition
  AtomRListList
      m_rotAtoms; // List of atom lists on smallest side of each rotable bond
  // DM 25 Apr 2002
  // m_flexIntns is a vector of size N(atoms)
  // Each element is an AtomList (i.e. std::vector<AtomPtr>)
  //
  //*** The vector is indexed by AtomId
  //*** We assume all atom IDs in the molecule are assigned sequentially from 1
  //***
  //*** For molecules read from MdlFileSource and PsfFileSource this is
  // guaranteed
  //
  // m_flexIntns[i-1] is a vector of all the atoms which can interact with atom
  // i across at least one rotable bond: (these interaction distances are
  // conformation dependent) This is a symmetric matrix - For a flexible
  // interaction i,j, atom j will appear in the AtomList for atom i, and vice
  // versa. Atoms which are NOT in this list are therefore in the same rigid
  // fragment as atom i: (these interaction distances are fixed and can
  // therefore
  // be ignored by any scoring function)
  AtomRListList m_flexIntns;

  // DM 2 Jul 2002 - if model has tethered atoms we need to modify the behaviour
  // of ModelMutator Could be subclassed if necessary
  AtomList m_tetheredAtoms;
};

// Useful typedefs
typedef SmartPtr<ModelMutator> ModelMutatorPtr; // Smart pointer

} // namespace rxdock

#endif //_RBTMODELMUTATOR_H_
