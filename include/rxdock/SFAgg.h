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

// Scoring Function aggregate class. Manages collection of child scoring
// functions.

#ifndef _RBTSFAGG_H_
#define _RBTSFAGG_H_

#include "rxdock/BaseSF.h"

namespace rxdock {

// Only check SF aggregate assertions in debug build
#ifdef _NDEBUG
const Bool SFAGG_CHECK = false;
#else
const bool SFAGG_CHECK = true;
#endif //_NDEBUG

// Useful typedefs

class SFAgg : public BaseSF {
public:
  // Static data member for class type (i.e. "SFAgg")
  static const std::string _CT;

  ////////////////////////////////////////
  // Constructors/destructors
  RBTDLL_EXPORT SFAgg(const std::string &strName = GetMetaDataPrefix() +
                                                   "score");
  virtual ~SFAgg();

  ////////////////////////////////////////
  // Public methods
  ////////////////

  // Returns all child component scores as a string-variant map
  // Key = fully qualified component name, value = weighted score
  //(for saving in a Model's data fields)
  virtual void ScoreMap(StringVariantMap &scoreMap) const;

  // Aggregate handling methods
  virtual void Add(BaseSF *);
  virtual void Remove(BaseSF *);
  virtual bool isAgg() const;
  virtual unsigned int GetNumSF() const;
  virtual BaseSF *GetSF(unsigned int iSF) const;

  // WorkSpace handling methods
  // Register scoring function with a workspace
  // Aggregate version registers all children, but NOT itself
  //(Aggregates are just containers, and have no need for model information
  virtual void Register(WorkSpace *);
  // Unregister with a workspace
  // Aggregate version unregisters all children, but NOT itself
  virtual void Unregister();

  // Override Observer pure virtual
  // Notify observer that subject has changed
  // Does nothing in SFAgg as aggregates do not require updating
  virtual void Update(Subject *theChangedSubject);

  // Request Handling method
  virtual void HandleRequest(RequestPtr spRequest);

  // Virtual function for dumping scoring function details to an output stream
  // Called by operator <<
  virtual void Print(std::ostream &s) const;

protected:
  ////////////////////////////////////////
  // Protected methods
  ///////////////////
  virtual double RawScore() const;

private:
  ////////////////////////////////////////
  // Private methods
  /////////////////
  SFAgg(const SFAgg &);            // Copy constructor disabled by default
  SFAgg &operator=(const SFAgg &); // Copy assignment disabled by default

protected:
  ////////////////////////////////////////
  // Protected data
  ////////////////

private:
  ////////////////////////////////////////
  // Private data
  //////////////
  BaseSFList m_sf;
  int m_nNonHLigandAtoms; // for normalised scores (score / non-H ligand atoms)
};

// Useful typedefs
typedef SmartPtr<SFAgg> SFAggPtr;        // Smart pointer
typedef std::vector<SFAggPtr> SFAggList; // Vector of smart pointers
typedef SFAggList::iterator SFAggListIter;
typedef SFAggList::const_iterator SFAggListConstIter;

} // namespace rxdock

#endif //_RBTSFAGG_H_
