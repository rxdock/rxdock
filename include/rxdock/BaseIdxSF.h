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

// Base class for indexed-grid-based scoring functions
// Provides protected methods for subclasses to use to create indexing grids

#ifndef _RBTBASEIDXSF_H_
#define _RBTBASEIDXSF_H_

#include "rxdock/BaseSF.h"
#include "rxdock/InteractionGrid.h"
#include "rxdock/NonBondedGrid.h"
#include "rxdock/NonBondedHHSGrid.h"

namespace rxdock {

class BaseIdxSF : public virtual BaseSF {
public:
  // Class type string
  static std::string _CT;
  // Parameter names
  static std::string _GRIDSTEP;
  static std::string _BORDER;

  ////////////////////////////////////////
  // Constructors/destructors
  virtual ~BaseIdxSF();

  ////////////////////////////////////////
  // Public methods
  ////////////////
  double GetGridStep() const;
  void SetGridStep(double step);
  double GetBorder() const;
  void SetBorder(double border);

protected:
  ////////////////////////////////////////
  // Protected methods
  ///////////////////
  BaseIdxSF();

  InteractionGridPtr CreateInteractionGrid() const;
  NonBondedGridPtr CreateNonBondedGrid() const;
  NonBondedHHSGridPtr CreateNonBondedHHSGrid() const;
  double GetMaxError() const;
  // DM 12 Apr 2002
  // Returns the maximum range of the scoring function,
  // corrected for max grid error, and grid border around docking site
  // This should be used by subclasses for selecting the receptor atoms to index
  // GetCorrectedRange() = GetRange() + GetMaxError() + GetBorder()
  double GetCorrectedRange() const;

  // As this has a virtual base class we need a separate OwnParameterUpdated
  // which can be called by concrete subclass ParameterUpdated methods
  // See Stroustrup C++ 3rd edition, p395, on programming virtual base classes
  void OwnParameterUpdated(const std::string &strName);

private:
  ////////////////////////////////////////
  // Private methods
  /////////////////

protected:
  ////////////////////////////////////////
  // Protected data
  ////////////////

private:
  ////////////////////////////////////////
  // Private data
  //////////////
  double m_gridStep;
  double m_border;
};

} // namespace rxdock

#endif //_RBTBASEIDXSF_H_
