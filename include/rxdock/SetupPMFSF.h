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

#ifndef _RBTSETUPPMFSF_H_
#define _RBTSETUPPMFSF_H_

#include "BaseInterSF.h"

namespace rxdock {

class SetupPMFSF : public BaseInterSF {
  AtomList theLigandList;   // ligand typing
  AtomList theReceptorList; // receptor typing

public:
  static std::string _CT;

  SetupPMFSF(const std::string &strName = "SETUP_PMF");
  ~SetupPMFSF();

protected:
  virtual void SetupReceptor();
  virtual void SetupLigand();
  virtual void SetupScore();
  virtual double RawScore() const;

  void SetupReceptorPMFTypes(void);
  void SetupLigandPMFTypes(void);
  /**
   * Routines to figure out complex PMF types
   */
  PMFType GetPMFfor_lC(AtomPtr); // for ligands
  PMFType GetPMFfor_lN(AtomPtr);
  PMFType GetPMFfor_lO(AtomPtr);
  PMFType GetPMFfor_lS(AtomPtr);

  PMFType GetPMFfor_rC(AtomPtr); // for receptor
  PMFType GetPMFfor_rN(AtomPtr);
  PMFType GetPMFfor_rO(AtomPtr);
  PMFType GetPMFfor_rS(AtomPtr);

  bool IsChargedNitrogen(AtomPtr); // true for guanidino or other charged
};

} // namespace rxdock

#endif //_RBTSETUPPMFSF_H_
