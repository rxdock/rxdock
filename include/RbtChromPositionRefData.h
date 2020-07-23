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

#include "RbtChromElement.h"
#include "RbtDockingSite.h"
#include "RbtEuler.h"
#include "RbtModel.h"

namespace rxdock {

class RbtChromPositionRefData {
public:
  // Class type string
  static std::string _CT;
  // Reference Cartesian axes
  static const RbtPrincipalAxes CARTESIAN_AXES;
  RbtChromPositionRefData(
      const RbtModel *pModel, const RbtDockingSite *pDockSite,
      double transStepSize, // Angstroms
      double rotStepSize,   // radians
      RbtChromElement::eMode transMode = RbtChromElement::FREE,
      RbtChromElement::eMode rotMode = RbtChromElement::FREE,
      double maxTrans = 0.0, // Angstroms
      double maxRot = 0.0);  // radians
  virtual ~RbtChromPositionRefData();

  int GetNumStartCoords() const { return m_startCoords.size(); }
  const RbtCoord &GetStartCoord(int iCoord) const {
    return m_startCoords[iCoord];
  }
  double GetTransStepSize() const { return m_transStepSize; }
  double GetRotStepSize() const { return m_rotStepSize; }
  RbtChromElement::eMode GetTransMode() const { return m_transMode; }
  RbtChromElement::eMode GetRotMode() const { return m_rotMode; }
  // Chromosome length, excluding FIXED modes (0, 3 or 6)
  int GetLength() const { return m_length; }
  // Chromosome length for crossover, excluding FIXED modes (0, 1 or 2)
  int GetXOverLength() const { return m_xOverLength; }
  bool IsTransFixed() const { return m_transMode == RbtChromElement::FIXED; }
  bool IsRotFixed() const { return m_rotMode == RbtChromElement::FIXED; }
  double GetMaxTrans() const { return m_maxTrans; }
  double GetMaxRot() const { return m_maxRot; }
  const RbtCoord &GetInitialCOM() const { return m_initialCom; }
  const RbtEuler &GetInitialOrientation() const { return m_initialOrientation; }
  const RbtQuat &GetInitialQuat() const { return m_initialQuat; }

  void GetModelValue(RbtCoord &com, RbtEuler &orientation) const;
  void SetModelValue(const RbtCoord &com, const RbtEuler &orientation);

private:
  RbtAtomList m_refAtoms;
  RbtAtomRList m_movableAtoms;
  RbtCoordList m_startCoords;
  double m_transStepSize;
  double m_rotStepSize;
  RbtCoord m_initialCom;
  RbtEuler m_initialOrientation;
  RbtQuat m_initialQuat;
  RbtChromElement::eMode m_transMode;
  RbtChromElement::eMode m_rotMode;
  int m_length;
  int m_xOverLength;
  // Max distance allowed from starting coord
  // Only used if m_transMode == TETHERED
  double m_maxTrans;
  // Max rot allowed from starting orientation
  // Only used if m_rotMode == TETHERED
  double m_maxRot;
};

typedef SmartPtr<RbtChromPositionRefData>
    RbtChromPositionRefDataPtr; // Smart pointer

} // namespace rxdock

#endif /*RBTCHROMPOSITIONREFDATA_H_*/
