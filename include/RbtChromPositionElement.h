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

// Chromosome element for model centre of mass and orientation
#ifndef RBTCHROMPOSITIONELEMENT_H_
#define RBTCHROMPOSITIONELEMENT_H_

#include "RbtChromElement.h"
#include "RbtChromPositionRefData.h"
#include "RbtRand.h"

namespace rxdock {

class RbtChromPositionElement : public RbtChromElement {
public:
  // Class type string
  static std::string _CT;

  RBTDLL_EXPORT RbtChromPositionElement(
      const RbtModel *pModel, const RbtDockingSite *pDockSite,
      double transStepSize, // Angstroms
      double rotStepSize,   // Radians
      RbtChromElement::eMode transMode = RbtChromElement::FREE,
      RbtChromElement::eMode rotMode = RbtChromElement::FREE,
      double maxTrans = 0.0, // Angstroms
      double maxRot = 0.0);  // radians
  virtual ~RbtChromPositionElement();
  virtual void Reset();
  virtual void Randomise();
  virtual void Mutate(double relStepSize);
  virtual void SyncFromModel();
  virtual void SyncToModel();
  virtual RbtChromElement *clone() const;
  virtual int GetLength() const { return m_spRefData->GetLength(); }
  virtual int GetXOverLength() const { return m_spRefData->GetXOverLength(); }
  virtual void GetVector(std::vector<double> &v) const;
  virtual void GetVector(RbtXOverList &v) const;
  virtual void SetVector(const std::vector<double> &v, int &i);
  virtual void SetVector(const RbtXOverList &v, int &i);
  virtual void GetStepVector(std::vector<double> &v) const;
  virtual double CompareVector(const std::vector<double> &v, int &i) const;
  virtual void Print(std::ostream &s) const;

  // Returns a standardised rotation angle in the range [-M_PI, +M_PI}
  // This function operates in radians
  static double StandardisedValue(double rotationAngle);

protected:
  // For use by clone()
  RbtChromPositionElement(RbtChromPositionRefDataPtr spRefData,
                          const RbtCoord &com, const RbtEuler &orientation);
  RbtChromPositionElement();
  void RandomiseCOM();
  void RandomiseOrientation();
  void MutateCOM(double relStepSize);
  void MutateOrientation(double relStepSize);
  void CorrectTetheredCOM();
  void CorrectTetheredOrientation();

private:
  RbtChromPositionRefDataPtr m_spRefData; // Fixed reference data
  RbtCoord m_com;                         // Centre of mass genotype value
  RbtEuler m_orientation; // Euler angle orientation genotype value
};

} // namespace rxdock

#endif /*RBTCHROMPOSITIONELEMENT_H_*/
