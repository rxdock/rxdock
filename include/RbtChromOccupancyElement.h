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

// Chromosome element for solvent model occupancy
#ifndef RBTCHROMOCCUPANCYELEMENT_H_
#define RBTCHROMOCCUPANCYELEMENT_H_

#include "RbtChromElement.h"
#include "RbtChromOccupancyRefData.h"
#include "RbtRand.h"

namespace rxdock {

class RbtChromOccupancyElement : public RbtChromElement {
public:
  // Class type string
  static std::string _CT;
  // Sole constructor
  RBTDLL_EXPORT RbtChromOccupancyElement(RbtModel *pModel, double stepSize,
                                         double threshold);

  virtual ~RbtChromOccupancyElement();
  virtual void Reset();
  virtual void Randomise();
  virtual void Mutate(double relStepSize);
  virtual void SyncFromModel();
  virtual void SyncToModel();
  virtual RbtChromElement *clone() const;
  virtual int GetLength() const { return 1; }
  virtual int GetXOverLength() const { return 1; }
  virtual void GetVector(std::vector<double> &v) const;
  virtual void GetVector(RbtXOverList &v) const;
  virtual void SetVector(const std::vector<double> &v, int &i);
  virtual void SetVector(const RbtXOverList &v, int &i);
  virtual void GetStepVector(std::vector<double> &v) const;
  virtual double CompareVector(const std::vector<double> &v, int &i) const;
  virtual void Print(std::ostream &s) const;

  // Returns a standardised occupancy value in the range [0,1]
  static double StandardisedValue(double occupancy);

protected:
  // For use by clone()
  RbtChromOccupancyElement(RbtChromOccupancyRefDataPtr spRefData, double value);
  RbtChromOccupancyElement();

private:
  RbtChromOccupancyRefDataPtr m_spRefData; // Fixed reference data
  double m_value;                          // The genotype value
};

} // namespace rxdock

#endif /*RBTCHROMOCCUPANCYELEMENT_H_*/
