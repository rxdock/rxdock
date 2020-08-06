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

// Chromosome factory class
// Creates an Chrom aggregate by visiting a number of
// FlexData subclasses.
// It is the users responsibility to delete the
// chromosome returned by GetChrom() when no longer required.
// No smart pointers are used within ChromFactory.
#ifndef _RBTCHROMFACTORY_H_
#define _RBTCHROMFACTORY_H_

#include "rxdock/FlexDataVisitor.h"
#include "rxdock/ModelMutator.h"

namespace rxdock {

class ChromElement;

class RBTDLL_EXPORT ChromFactory : public FlexDataVisitor {
public:
  ChromFactory();
  virtual void VisitReceptorFlexData(ReceptorFlexData *);
  virtual void VisitLigandFlexData(LigandFlexData *);
  virtual void VisitSolventFlexData(SolventFlexData *);

  ChromElement *GetChrom() const;
  // Temporary solution whilst we replace the ModelMutator class
  ModelMutatorPtr GetModelMutator() const { return m_spMutator; }

private:
  ChromElement *m_pChrom;
  ModelMutatorPtr m_spMutator;
};

} // namespace rxdock

#endif //_RBTCHROMFACTORY_H_
