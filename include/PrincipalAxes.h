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

// Principal axes calculation routines (in rxdock namespace)

#ifndef _RBTPRINCIPALAXES_H_
#define _RBTPRINCIPALAXES_H_

#include "Atom.h"
#include "Config.h"

namespace rxdock {

// Struct for holding principal axes info
// Center of mass, three principal axis vectors, three principal moments
class PrincipalAxes {
public:
  // Default is X,Y,Z cartesian axes centered at origin
  PrincipalAxes()
      : com(0.0, 0.0, 0.0), axis1(1.0, 0.0, 0.0), axis2(0.0, 1.0, 0.0),
        axis3(0.0, 0.0, 1.0), moment1(1.0), moment2(1.0), moment3(1.0) {}
  Coord com;
  Vector axis1;
  Vector axis2;
  Vector axis3;
  double moment1;
  double moment2;
  double moment3;
};

typedef std::vector<PrincipalAxes> PrincipalAxesList;
typedef PrincipalAxesList::iterator PrincipalAxesListIter;
typedef PrincipalAxesList::const_iterator PrincipalAxesListConstIter;

// Calculates principal axes and center of mass for the atoms in the atom list
RBTDLL_EXPORT PrincipalAxes GetPrincipalAxesOfAtoms(const AtomList &atomList);
// Calculates principal axes and center of mass for the coords in the coord list
// (assumes all masses=1)
RBTDLL_EXPORT PrincipalAxes GetPrincipalAxesOfAtoms(const CoordList &coordList);
// Special case for water
PrincipalAxes GetSolventPrincipalAxes(const AtomPtr &oAtom,
                                      const AtomPtr &h1Atom,
                                      const AtomPtr &h2Atom);
// DM 17 Jul 2001 - returns the quaternion used to effect the transformation
Quat AlignPrincipalAxesOfAtoms(AtomList &atomList,
                               const PrincipalAxes &alignAxes = PrincipalAxes(),
                               bool bAlignCOM = true);
// Returns the quaternion required to align principal axes with reference axes
RBTDLL_EXPORT Quat GetQuatFromAlignAxes(const PrincipalAxes &prAxes,
                                        const PrincipalAxes &refAxes);
// Returns the quaternion required to effect an alignment of v onto ref vector.
// Vectors do not have to be unit length.
// BadArgument exception is thrown if either v or ref is zero length.
Quat GetQuatFromAlignVectors(const Vector &v, const Vector &ref);

} // namespace rxdock

#endif //_RBTPRINCIPALAXES_H_
