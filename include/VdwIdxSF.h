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

// Indexed-grid-based intermolecular scoring function for vdw interactions
#ifndef _RBTVDWIDXSF_H_
#define _RBTVDWIDXSF_H_

#include "BaseIdxSF.h"
#include "BaseInterSF.h"
#include "VdwSF.h"

namespace rxdock {

class VdwIdxSF : public BaseInterSF, public BaseIdxSF, public VdwSF {
public:
  // Class type string
  static std::string _CT;
  // Parameter names
  // DM 12 Jun 2002 - score thresholds used for counting attractive and
  // repulsive interactions (thresholds represent partial vdW scores summed per
  // ligand atom, used for outputting ligand atom counts)
  static std::string _THRESHOLD_ATTR;
  static std::string _THRESHOLD_REP;
  // DM 03 Feb 2003 - score threshold for outputting individual (atom-atom)
  // lipophilic vdW interactions i.e between two non-polar carbons/hydrogens
  static std::string _ANNOTATION_LIPO;
  // DM 10 Apr 2003 - vdw annotation can be quite slow so provide boolean option
  // to disable
  static std::string _ANNOTATE;
  // DM 14 Jun 2006 - option to disable solvent performance enhancements, mainly
  // for testing
  static std::string _FAST_SOLVENT;

  RBTDLL_EXPORT VdwIdxSF(const std::string &strName = "VDW");
  virtual ~VdwIdxSF();

  // Override BaseSF::ScoreMap to provide additional raw descriptors
  virtual void ScoreMap(StringVariantMap &scoreMap) const;

protected:
  virtual void SetupReceptor();
  virtual void SetupLigand();
  virtual void SetupSolvent();
  virtual void SetupScore();
  virtual double RawScore() const;
  double InterScore() const;
  double ReceptorScore() const;
  double SolventScore() const;
  double ReceptorSolventScore() const;
  double LigandSolventScore() const;

  // DM 25 Oct 2000 - track changes to parameter values in local data members
  // ParameterUpdated is invoked by ParamHandler::SetParameter
  void ParameterUpdated(const std::string &strName);

private:
  void RenderAnnotationsByResidue(std::vector<std::string> &retVal) const;

  NonBondedGridPtr m_spGrid;        // Indexing grid for receptor
  NonBondedGridPtr m_spSolventGrid; // Indexing grid for fixed/tethered solvent
  AtomList m_recAtomList;
  AtomRList m_recRigidAtomList;
  AtomRList m_recFlexAtomList;
  AtomRList m_ligAtomList;
  AtomRList m_solventAtomList;        // All solvent atoms
  AtomRList m_solventFixTethAtomList; // Fixed or tethered solvent atoms
  AtomRList m_solventFreeAtomList;    // Free solvent atoms
  AtomRListList m_recFlexIntns;       // Intra-target flexible interactions
  AtomRListList
      m_recFlexPrtIntns; // Intra-target flexible partitioned interactions
  AtomRListList
      m_solventFixTethIntns; // Intra-solvent intns between fixed/tethered atoms
  AtomRListList m_solventFixTethPrtIntns; // Partitioned intns between
                                          // fixed/tethered solvent
  AtomRListList
      m_solventFreeIntns; // Intra-solvent intns between free solvent atoms
  // DM 12 Jun 2002 - keep track of number of ligand atoms involved in non-zero
  // vdW interactions
  mutable int m_nAttr; //#atoms with net attractive (-ve) vdw scores
  mutable int m_nRep;  //#atoms with net repulsive (+ve) vdw scores
  double m_attrThreshold;
  double m_repThreshold;
  double m_lipoAnnot;
  bool m_bAnnotate;
  bool m_bFlexRec;
  bool m_bFastSolvent;
};

} // namespace rxdock

#endif //_RBTVDWIDXSF_H_
