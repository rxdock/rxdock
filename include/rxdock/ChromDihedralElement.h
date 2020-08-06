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

#include "rxdock/ChromDihedralRefData.h"
#include "rxdock/ChromElement.h"
#include "rxdock/Rand.h"

namespace rxdock {

class ChromDihedralElement : public ChromElement {
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
  RBTDLL_EXPORT ChromDihedralElement(
      BondPtr spBond,         // Rotatable bond
      AtomList tetheredAtoms, // Tethered atom list
      double stepSize,        // maximum mutation step size (degrees)
      ChromElement::eMode mode = ChromElement::FREE, // sampling mode
      double maxDihedral =
          0.0); // max deviation from reference (tethered mode only)

  virtual ~ChromDihedralElement();
  virtual void Reset();
  virtual void Randomise();
  virtual void Mutate(double relStepSize);
  virtual void SyncFromModel();
  virtual void SyncToModel();
  virtual ChromElement *clone() const;
  virtual int GetLength() const { return 1; }
  virtual int GetXOverLength() const { return 1; }
  virtual void GetVector(std::vector<double> &v) const;
  virtual void GetVector(XOverList &v) const;
  virtual void SetVector(const std::vector<double> &v, int &i);
  virtual void SetVector(const XOverList &v, int &i);
  virtual void GetStepVector(std::vector<double> &v) const;
  virtual double CompareVector(const std::vector<double> &v, int &i) const;
  virtual void Print(std::ostream &s) const;

  // Returns a standardised dihedral angle in the range [-180, +180}
  // This function operates in degrees
  static double StandardisedValue(double dihedralAngle);

protected:
  // For use by clone()
  ChromDihedralElement(ChromDihedralRefDataPtr spRefData, double value);
  ChromDihedralElement();
  void CorrectTetheredDihedral();

private:
  ChromDihedralRefDataPtr m_spRefData; // Fixed reference data
  double m_value;                      // The genotype value
};

} // namespace rxdock

#endif /*RBTCHROMDIHEDRALELEMENT_H_*/
