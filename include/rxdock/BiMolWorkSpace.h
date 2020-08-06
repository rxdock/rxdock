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

// Subclass of WorkSpace for simulations involving
// Receptor = Model #0
// Ligand   = Model #1
// Solvent  = Model #2 upwards
#ifndef _RBTBIMOLWORKSPACE_H_
#define _RBTBIMOLWORKSPACE_H_

#include "Population.h"
#include "WorkSpace.h"

namespace rxdock {

class BiMolWorkSpace : public WorkSpace {
  enum eModelID {
    RECEPTOR = 0,
    LIGAND = 1,
    SOLVENT = 2 // Flexible waters are represented as separate models, from 2
  };

public:
  ////////////////////////////////////////
  // Constructors/destructors

  RBTDLL_EXPORT BiMolWorkSpace();
  virtual ~BiMolWorkSpace();

  ////////////////////////////////////////
  // Public methods
  ////////////////

  // Model handling
  RBTDLL_EXPORT ModelPtr GetReceptor() const;
  RBTDLL_EXPORT ModelPtr GetLigand() const;
  RBTDLL_EXPORT ModelList GetSolvent() const;
  RBTDLL_EXPORT bool hasSolvent() const;
  RBTDLL_EXPORT void SetReceptor(ModelPtr spReceptor);
  RBTDLL_EXPORT void SetLigand(ModelPtr spLigand);
  RBTDLL_EXPORT void SetSolvent(ModelList solventList);
  RBTDLL_EXPORT void RemoveSolvent();
  RBTDLL_EXPORT void
  UpdateModelCoordsFromChromRecords(BaseMolecularFileSource *pSource);

  // Model I/O
  // Saves ligand to file sink
  virtual void Save(bool bSaveScores = true);
  virtual void SaveHistory(bool bSaveScores = true);

protected:
  ////////////////////////////////////////
  // Protected methods
  ///////////////////

private:
  ////////////////////////////////////////
  // Private methods
  /////////////////

  BiMolWorkSpace(
      const BiMolWorkSpace &); // Copy constructor disabled by default
  BiMolWorkSpace &
  operator=(const BiMolWorkSpace &); // Copy assignment disabled by default
  void SaveLigand(MolecularFileSinkPtr spSink, bool bSaveScores = true);

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
typedef SmartPtr<BiMolWorkSpace> BiMolWorkSpacePtr; // Smart pointer

} // namespace rxdock

#endif //_RBTBIMOLWORKSPACE_H_
