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

// Factory class for generating workspace objects from <receptor>.prm file
// Responsible for generation of:
// Receptor, Ligand, Solvent,
#ifndef RBTPRMFACTORY_H_
#define RBTPRMFACTORY_H_

#include "rxdock/BaseMolecularFileSource.h"
#include "rxdock/Error.h"
#include "rxdock/Model.h"

namespace rxdock {

class ParameterFileSource;
class DockingSite;

class PRMFactory {
public:
  static const std::string &_CT;
  // Receptor parameters
  static const std::string &_REC_SECTION;
  static const std::string &_REC_FILE;
  static const std::string &_REC_TOPOL_FILE;
  static const std::string &_REC_COORD_FILE;
  static const std::string &_REC_NUM_COORD_FILES;
  static const std::string &_REC_FLEX_DISTANCE;
  static const std::string &_REC_DIHEDRAL_STEP;
  static const std::string &_REC_ALL_H;
  static const std::string &_REC_MASSES_FILE;
  static const std::string &_REC_SEGMENT_NAME;

  // Ligand parameters
  static const std::string &_LIG_SECTION;
  // Ligand file is not read here, so no need for _LIG_FILE
  // No need to specify ligand flexibility parameter names explicitly
  // as we assume that all parameters in SECTION LIGAND
  // are flexibility params.

  // Solvent parameters
  static const std::string &_SOLV_SECTION;
  static const std::string &_SOLV_FILE;
  // No need to specify solvent flexibility parameter names explicitly
  // as we assume that all parameters in SECTION SOLVENT
  // are flexibility params (except for _SOLV_FILE)

  RBTDLL_EXPORT PRMFactory(ParameterFileSource *pParamSource);
  RBTDLL_EXPORT PRMFactory(ParameterFileSource *pParamSource, DockingSite *pDS);

  DockingSite *GetDockingSite() const { return m_pDS; }
  void SetDockingSite(DockingSite *pDS) { m_pDS = pDS; }

  RBTDLL_EXPORT ModelPtr CreateReceptor();
  RBTDLL_EXPORT ModelPtr CreateLigand(BaseMolecularFileSource *pSource);
  RBTDLL_EXPORT ModelList CreateSolvent();

private:
  // Creates the appropriate source according to the file extension
  MolecularFileSourcePtr CreateMolFileSource(const std::string &fileName);
  void AttachReceptorFlexData(Model *pReceptor);
  void AttachLigandFlexData(Model *pLigand);
  void AttachSolventFlexData(Model *pSolvent);
  ParameterFileSource *m_pParamSource;
  DockingSite *m_pDS;
};

} // namespace rxdock

#endif /*RBTPRMFACTORY_H_*/
