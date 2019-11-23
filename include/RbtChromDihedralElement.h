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

// Chromosome element for a single dihedral angle
#ifndef RBTCHROMDIHEDRALELEMENT_H_
#define RBTCHROMDIHEDRALELEMENT_H_

#include "RbtChromDihedralRefData.h"
#include "RbtChromElement.h"
#include "RbtRand.h"

class RbtChromDihedralElement : public RbtChromElement {
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
  RBTDLL_EXPORT RbtChromDihedralElement(
      RbtBondPtr spBond,         // Rotatable bond
      RbtAtomList tetheredAtoms, // Tethered atom list
      double stepSize,           // maximum mutation step size (degrees)
      RbtChromElement::eMode mode = RbtChromElement::FREE, // sampling mode
      double maxDihedral =
          0.0); // max deviation from reference (tethered mode only)

  virtual ~RbtChromDihedralElement();
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

  // Returns a standardised dihedral angle in the range [-180, +180}
  // This function operates in degrees
  static double StandardisedValue(double dihedralAngle);

protected:
  // For use by clone()
  RbtChromDihedralElement(RbtChromDihedralRefDataPtr spRefData, double value);
  RbtChromDihedralElement();
  void CorrectTetheredDihedral();

private:
  RbtChromDihedralRefDataPtr m_spRefData; // Fixed reference data
  double m_value;                         // The genotype value
};

#endif /*RBTCHROMDIHEDRALELEMENT_H_*/
