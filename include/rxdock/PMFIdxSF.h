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

// Indexed-grid-based PMF intermolecular scoring function

#ifndef _RBTPMFIDXSF_H_
#define _RBTPMFIDXSF_H_

#include "rxdock/AnnotationHandler.h"
#include "rxdock/BaseIdxSF.h"
#include "rxdock/BaseInterSF.h"
#include "rxdock/PMF.h"
// for PMF pseudogrids where
// X is distance, X and Y are receptor and ligand
// distances respectively
#include "rxdock/RealGrid.h"

namespace rxdock {

/** PMFIdxSF class for PMF scoring.
 */
class PMFIdxSF : public BaseInterSF,
                 public BaseIdxSF,
                 public virtual AnnotationHandler {
  // by default "private" members
  /**
   * vGrid: vector of grids. every PMF type needs a grid
   */
  AtomList theLigandList;   // ligand typing
  AtomList theReceptorList; // recepor typing
  // regular lists
  AtomRList theLigandRList;   // ligand
  AtomRList theReceptorRList; // recepor
  std::vector<RealGridPtr>
      theTypeGrid; // grids for PMF values for different atom types in receptor
  std::vector<PMFType> theLigandTypes; // type values in theTypeGrid
  NonBondedGridPtr theSurround;        // atoms arond a gridpoint
  RealGridPtr thePMFGrid;              // grid for X-distance Y
                                       // this is the representation of the PMFs
  RealGridPtr theSlopeGrid; // grid to store values where the plateaus starts

public:
  PMFIdxSF(const std::string &strName = "PMF"); /**< The only one constructor */
  virtual ~PMFIdxSF();                          /**< The virtual destructor */

  /**
   * This string serves as a name for the weight factor
   * of the PMF score. Should be read from the .prm
   * file.
   */
  static std::string _CT;
  static std::string _PMFDIR; /**< String pointing to the dir with PMF tables*/
  static std::string _CC_CUTOFF;
  static std::string _SLOPE;
  /**
   *  Observer functionalities pure virtual, so it must be
   *  overloaded.
   */
  virtual void Update(Subject *theChangedSubject);

protected:
  /**
   * Setting up receptor means atom typing for all the
   * receptor atoms and create NonBondedGrid for NN
   */
  virtual void SetupReceptor();
  /**
   * determine PMF types for all ligand atoms
   */
  virtual void SetupLigand();
  /**
   * Setting up scores
   */
  virtual void SetupScore();
  /**
   * RawScore(void) returning with the actual PMF score
   */
  virtual double RawScore(void) const;
  /**
   * Estimate value for short distances instead of using plateau in PMFs
   */
  double GetLinearCloseRangeValue(double aDist, PMFType aRecType,
                                  PMFType aLigType) const;
};

} // namespace rxdock

#endif // _RBTPMFIDXSF_H_
