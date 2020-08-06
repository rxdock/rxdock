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

// Manages the fixed reference data for a single dihedral angle chromosome
// element Also provides methods to map the genotype (dihedral angle value) onto
// the phenotype (model coords) A single instance is designed to be shared
// between all clones of a given element
#ifndef RBTCHROMDIHEDRALREFDATA_H_
#define RBTCHROMDIHEDRALREFDATA_H_

#include "Atom.h"
#include "Bond.h"
#include "ChromElement.h"

namespace rxdock {

class ChromDihedralRefData {
public:
  // Class type string
  static std::string _CT;
  // Sole constructor
  // If the tetheredAtoms list is empty, then
  //  the end of the bond with the fewest pendant atoms is rotated (other half
  //  remains fixed)
  // else if the tetheredAtoms list is not empty, then
  //  the end of the bond with the fewest tethered atoms is rotated (other half
  //  remains fixed)
  ChromDihedralRefData(
      BondPtr spBond,         // Rotatable bond
      AtomList tetheredAtoms, // Tethered atom list
      double stepSize,        // maximum mutation step size (degrees)
      ChromElement::eMode mode = ChromElement::FREE, // sampling mode
      double maxDihedral =
          0.0); // max deviation from reference (tethered mode only)
  virtual ~ChromDihedralRefData();

  // Gets the maximum step size for this bond
  double GetStepSize() const { return m_stepSize; }
  // Gets the sampling mode for this bond
  ChromElement::eMode GetMode() const { return m_mode; }
  // Gets the maximum deviation from the reference dihedral (tethered mode only)
  double GetMaxDihedral() const { return m_maxDihedral; }
  // Gets the current dihedral angle (phenotype) for this bond
  // from the model coords
  double GetModelValue() const;
  // Sets the phenotype (model coords) for this bond
  // to a given dihedral angle
  void SetModelValue(double dihedralAngle);
  // Gets the initial dihedral angle for this bond
  //(initialised from model coords in ChromDihedralRefData constructor)
  double GetInitialValue() const { return m_initialValue; }

private:
  // Sets up the movable atom list for this bond
  void Setup(BondPtr spBond, const AtomList &tetheredAtoms);
  Atom *m_atom1;
  Atom *m_atom2;
  Atom *m_atom3;
  Atom *m_atom4;
  AtomRList m_rotAtoms;
  double m_stepSize;
  double m_initialValue;
  ChromElement::eMode m_mode;
  double m_maxDihedral; // max deviation from reference (tethered mode only)
};

typedef SmartPtr<ChromDihedralRefData> ChromDihedralRefDataPtr; // Smart pointer

} // namespace rxdock

#endif /*RBTCHROMDIHEDRALREFDATA_H_*/
