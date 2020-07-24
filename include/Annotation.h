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

#include "Atom.h"
#include "Config.h"

namespace rxdock {

class Annotation {
public:
  ////////////////////////////////////////
  // Constructors/destructors
  Annotation(const Atom *pAtom1, const Atom *pAtom2, double dist, double score);
  virtual ~Annotation();

  ////////////////////////////////////////
  // Public methods
  ////////////////
  // Get
  const Atom *GetAtom1Ptr() const;
  const Atom *GetAtom2Ptr() const;
  double GetDistance() const;
  double GetScore() const;
  // Get the fully qualified (FQ) residue name for atom 2 (target atom)
  std::string GetFQResName() const;

  // Set
  void SetAtom1Ptr(const Atom *pAt1);
  void SetAtom2Ptr(const Atom *pAt2);
  void SetDistance(double d);
  void SetScore(double s);

  // Render annotation into string in rDock Viewer format
  std::string Render() const;

  // Operators
  // Special meaning of operator+ for accumulating annotations by residue
  void operator+=(const Annotation &ann);

protected:
  ////////////////////////////////////////
  // Protected methods
  ///////////////////

private:
  ////////////////////////////////////////
  // Private methods
  /////////////////
  Annotation(); // Disable default constructor
  // Annotation(const Annotation&);//Copy constructor disabled by default
  // Annotation& operator=(const Annotation&);//Copy assignment disabled
  // by default

protected:
  ////////////////////////////////////////
  // Protected data
  ////////////////

private:
  ////////////////////////////////////////
  // Private data
  //////////////
  const Atom *m_pAtom1;
  const Atom *m_pAtom2;
  double m_dist;
  double m_score;
};

// Useful typedefs
typedef SmartPtr<Annotation> AnnotationPtr;        // Smart pointer
typedef std::vector<AnnotationPtr> AnnotationList; // Vector of smart pointers
typedef AnnotationList::iterator AnnotationListIter;
typedef AnnotationList::const_iterator AnnotationListConstIter;

///////////////////////////////////////////////
// Non-member functions (in rxdock namespace)
//////////////////////////////////////////

////////////////////////////////////////////////////////
// Comparison functions for sorting Annotation* containers
// For use by STL sort algorithms
////////////////////////////////////////////////////////

// Less than operator for sorting Annotation*'s by the atom ID of atom 2
// (target)
class Ann_Cmp_AtomId2 {
public:
  bool operator()(Annotation *pAnn1, Annotation *pAnn2) const {
    return pAnn1->GetAtom2Ptr()->GetAtomId() <
           pAnn2->GetAtom2Ptr()->GetAtomId();
  }
};

} // namespace rxdock

#endif //_RBTANNOTATION_H_
