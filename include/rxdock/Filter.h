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

// Reads a file that contains a list of filters that determine
// when to terminate with a given ligand, and another filter to determine
// when to write the current conformation to the output.

#ifndef _RBT_FILTER_H_
#define _RBT_FILTER_H_

#include "rxdock/BaseObject.h"
#include "rxdock/Context.h"
#include "rxdock/FilterExpression.h"

namespace rxdock {

class FilterExpressionVisitor;

class Filter : public BaseObject {
public:
  static std::string _CT;
  RBTDLL_EXPORT Filter(std::string strfilter, bool filter = false);
  ///////////////////
  // Destructor
  //////////////////
  virtual ~Filter();
  // Override Observer pure virtual
  // Notify observer that subject has changed
  virtual void Update(Subject *theChangedSubject);

  void SetupReceptor();       // Called by Update when receptor is changed
  void SetupLigand();         // Called by Update when ligand is changed
  void SetupScore();          // Called by Update when either model has changed
  RBTDLL_EXPORT bool Write(); // Output conformation?
  RBTDLL_EXPORT bool Terminate(); // Finished with ligand?
  ModelPtr GetReceptor() const;
  ModelPtr GetLigand() const;
  void SetMaxNRuns(int n) { maxnruns = n; }

  ////////////////////
  // Private methods
  ////////////////////

protected:
  Filter(); // Default constructor disabled

private:
  void ReadFilters();

  int filteridx, nTermFilters, nWriteFilters, nruns, maxnruns;
  FilterExpressionList terminationFilters;
  FilterExpressionList writtingFilter;
  ModelPtr m_spReceptor;
  ModelPtr m_spLigand;
  ContextPtr contextp;
};

// Useful typedefs
typedef SmartPtr<Filter> FilterPtr; // Smart pointer
typedef std::vector<FilterPtr> FilterList;
// Vector of smart pointers
typedef FilterList::iterator FilterListIter;
typedef FilterList::const_iterator FilterListConstIter;

} // namespace rxdock

#endif //_Filter_H_
