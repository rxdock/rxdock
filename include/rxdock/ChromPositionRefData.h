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

// Manages the fixed reference data for a position chromosome element
// Also provides methods to map the genotype (COM and Euler angles) onto the
// phenotype (model coords)
// A single instance is designed to be shared between all clones of a given
// element
#ifndef RBTCHROMPOSITIONREFDATA_H_
#define RBTCHROMPOSITIONREFDATA_H_

#include "rxdock/ChromElement.h"
#include "rxdock/DockingSite.h"
#include "rxdock/Euler.h"
#include "rxdock/Model.h"

namespace rxdock {

class ChromPositionRefData {
public:
  // Class type string
  static std::string _CT;
  // Reference Cartesian axes
  static const PrincipalAxes CARTESIAN_AXES;
  ChromPositionRefData(const Model *pModel, const DockingSite *pDockSite,
                       double transStepSize, // Angstroms
                       double rotStepSize,   // radians
                       ChromElement::eMode transMode = ChromElement::FREE,
                       ChromElement::eMode rotMode = ChromElement::FREE,
                       double maxTrans = 0.0, // Angstroms
                       double maxRot = 0.0);  // radians
  virtual ~ChromPositionRefData();

  int GetNumStartCoords() const { return m_startCoords.size(); }
  const Coord &GetStartCoord(int iCoord) const { return m_startCoords[iCoord]; }
  double GetTransStepSize() const { return m_transStepSize; }
  double GetRotStepSize() const { return m_rotStepSize; }
  ChromElement::eMode GetTransMode() const { return m_transMode; }
  ChromElement::eMode GetRotMode() const { return m_rotMode; }
  // Chromosome length, excluding FIXED modes (0, 3 or 6)
  int GetLength() const { return m_length; }
  // Chromosome length for crossover, excluding FIXED modes (0, 1 or 2)
  int GetXOverLength() const { return m_xOverLength; }
  bool IsTransFixed() const { return m_transMode == ChromElement::FIXED; }
  bool IsRotFixed() const { return m_rotMode == ChromElement::FIXED; }
  double GetMaxTrans() const { return m_maxTrans; }
  double GetMaxRot() const { return m_maxRot; }
  const Coord &GetInitialCOM() const { return m_initialCom; }
  const Euler &GetInitialOrientation() const { return m_initialOrientation; }
  const Quat &GetInitialQuat() const { return m_initialQuat; }

  void GetModelValue(Coord &com, Euler &orientation) const;
  void SetModelValue(const Coord &com, const Euler &orientation);

private:
  AtomList m_refAtoms;
  AtomRList m_movableAtoms;
  CoordList m_startCoords;
  double m_transStepSize;
  double m_rotStepSize;
  Coord m_initialCom;
  Euler m_initialOrientation;
  Quat m_initialQuat;
  ChromElement::eMode m_transMode;
  ChromElement::eMode m_rotMode;
  int m_length;
  int m_xOverLength;
  // Max distance allowed from starting coord
  // Only used if m_transMode == TETHERED
  double m_maxTrans;
  // Max rot allowed from starting orientation
  // Only used if m_rotMode == TETHERED
  double m_maxRot;
};

typedef SmartPtr<ChromPositionRefData> ChromPositionRefDataPtr; // Smart pointer

} // namespace rxdock

#endif /*RBTCHROMPOSITIONREFDATA_H_*/
