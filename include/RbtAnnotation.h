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

// Simple struct for holding annotation information associated with
// intermolecular scores

#ifndef _RBTANNOTATION_H_
#define _RBTANNOTATION_H_

#include "RbtAtom.h"
#include "RbtConfig.h"

namespace rxdock {

class RbtAnnotation {
public:
  ////////////////////////////////////////
  // Constructors/destructors
  RbtAnnotation(const RbtAtom *pAtom1, const RbtAtom *pAtom2, double dist,
                double score);
  virtual ~RbtAnnotation();

  ////////////////////////////////////////
  // Public methods
  ////////////////
  // Get
  const RbtAtom *GetAtom1Ptr() const;
  const RbtAtom *GetAtom2Ptr() const;
  double GetDistance() const;
  double GetScore() const;
  // Get the fully qualified (FQ) residue name for atom 2 (target atom)
  std::string GetFQResName() const;

  // Set
  void SetAtom1Ptr(const RbtAtom *pAt1);
  void SetAtom2Ptr(const RbtAtom *pAt2);
  void SetDistance(double d);
  void SetScore(double s);

  // Render annotation into string in rDock Viewer format
  std::string Render() const;

  // Operators
  // Special meaning of operator+ for accumulating annotations by residue
  void operator+=(const RbtAnnotation &ann);

protected:
  ////////////////////////////////////////
  // Protected methods
  ///////////////////

private:
  ////////////////////////////////////////
  // Private methods
  /////////////////
  RbtAnnotation(); // Disable default constructor
  // RbtAnnotation(const RbtAnnotation&);//Copy constructor disabled by default
  // RbtAnnotation& operator=(const RbtAnnotation&);//Copy assignment disabled
  // by default

protected:
  ////////////////////////////////////////
  // Protected data
  ////////////////

private:
  ////////////////////////////////////////
  // Private data
  //////////////
  const RbtAtom *m_pAtom1;
  const RbtAtom *m_pAtom2;
  double m_dist;
  double m_score;
};

// Useful typedefs
typedef SmartPtr<RbtAnnotation> RbtAnnotationPtr; // Smart pointer
typedef std::vector<RbtAnnotationPtr>
    RbtAnnotationList; // Vector of smart pointers
typedef RbtAnnotationList::iterator RbtAnnotationListIter;
typedef RbtAnnotationList::const_iterator RbtAnnotationListConstIter;

///////////////////////////////////////////////
// Non-member functions (in rxdock namespace)
//////////////////////////////////////////

////////////////////////////////////////////////////////
// Comparison functions for sorting RbtAnnotation* containers
// For use by STL sort algorithms
////////////////////////////////////////////////////////

// Less than operator for sorting RbtAnnotation*'s by the atom ID of atom 2
// (target)
class RbtAnn_Cmp_AtomId2 {
public:
  bool operator()(RbtAnnotation *pAnn1, RbtAnnotation *pAnn2) const {
    return pAnn1->GetAtom2Ptr()->GetAtomId() <
           pAnn2->GetAtom2Ptr()->GetAtomId();
  }
};

} // namespace rxdock

#endif //_RBTANNOTATION_H_
