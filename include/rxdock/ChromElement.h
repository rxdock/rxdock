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

// Abstract interface for all chromosome elements (degrees of freedom)
#ifndef RBTCHROMELEMENT_H_
#define RBTCHROMELEMENT_H_

#include "Config.h"
#include "Rand.h"

namespace rxdock {

// Typedefs for crossover data.
// To prevent splitting (for example) an orientation or position vector in two
// during crossover, we convert the chromosome to a vector of vector of doubles
// An orientation vector (3 Euler angles) is a single XOverElement
// GetVector(std::vector<double>&) returns a flat vector of doubles (for use by
// Simplex) GetVector(XOverList&) returns a vector of XOverElements (for
// use by crossover)
typedef std::vector<double> XOverElement;
typedef XOverElement::iterator XOverElementIter;
typedef XOverElement::const_iterator XOverElementConstIter;
typedef std::vector<XOverElement> XOverList;
typedef XOverList::iterator XOverListIter;
typedef XOverList::const_iterator XOverListConstIter;

class ChromElement {
public:
  // Class type string
  static std::string _CT;
  // Threshold used to assess equality of two chromosome elements
  static double _THRESHOLD;

  RBTDLL_EXPORT static double &GetThreshold();

  // General enum for defining types of chromosome element flexibility.
  // For use by subclasses of ChromElement if desired.
  // There is no functionality in the base class that uses this enum.
  enum eMode { FIXED = 0, TETHERED = 1, FREE = 2 };
  // Static methods to convert from mode enum to string and vice versa
  static eMode StrToMode(const std::string &modeStr); // case insensitive
  static std::string RBTDLL_EXPORT
  ModeToStr(eMode mode); // returns "FIXED", "TETHERED" or "FREE"ss

  virtual ~ChromElement();
  //
  // INTERFACE METHODS - must be implemented in subclasses
  //
  // Resets the chromosome element to its initial value (extracted from
  // model coords at time of construction)
  virtual void Reset() = 0;
  // Randomises the chromosome element within the context of any restraints
  virtual void Randomise() = 0;
  // Mutates the chromosome element
  // relStepSize = fraction (0.0->1.0) of maximum defined step size
  virtual void Mutate(double relStepSize) = 0;
  // Updates the chromosome element to match the current model coords
  virtual void SyncFromModel() = 0;
  // Updates the model coords to match the chromosome element
  virtual void SyncToModel() = 0;
  // Creates an independent clone
  // The current value of the cloned element should match the current value
  // of this element.
  virtual ChromElement *clone() const = 0;
  // Gets the number of double values needed to represent this
  // chromosome element
  virtual int GetLength() const = 0;
  // Gets the number of vectors of double values needed to represent this
  // chromosome element for crossover purposes
  // e.g. an orientation vector consists of 3 doubles, but should be crossed
  // over as a single intact entity
  virtual int GetXOverLength() const = 0;
  // Converts chromosome element to a vector of double values,
  // Number of double values appended should match GetLength().
  // v = vector of doubles to append to.
  virtual void GetVector(std::vector<double> &v) const = 0;
  // Converts chromosome element to a vector of XOverElements
  // where each XOverElement is itself a vector of double values
  // v = vector of XOverElements to append to.
  virtual void GetVector(XOverList &v) const = 0;
  // Updates chromosome element from a vector of double values,
  // Number of double values read should match GetLength().
  // v = vector of doubles to extract from
  // i = index of next vector element to read (should be updated by method)
  virtual void SetVector(const std::vector<double> &v, int &i) = 0;
  // Updates chromosome element from a vector of XOverElements,
  // Number of XOverElements read should match GetXOverLength().
  // v = vector of XOverElements to extract from
  // i = index of next vector element to read (should be updated by method)
  virtual void SetVector(const XOverList &v, int &i) = 0;
  // Gets the vector of absolute step sizes that correspond to each double
  // value.
  virtual void GetStepVector(std::vector<double> &v) const = 0;
  // Gets the maximum relative difference between this element and another
  // element as represented by the vector. Differences should be normalised by
  // the step size. This method is provided to allow better comparisons between
  // values that are difficult to compare by simple numerical differences. e.g.
  // Dihedral angles are cyclical, therefore -180 and + 179 only differ by 1
  // deg.
  virtual double CompareVector(const std::vector<double> &v, int &i) const = 0;
  //
  // IMPLEMENTED VIRTUAL METHODS
  //
  // Invalid operation in base class
  virtual void Add(ChromElement *pChromElement);
  // Prints details of element to stream (null implementation in base class)
  virtual void Print(std::ostream &s) const {}
  //
  // NON-VIRTUAL METHODS
  //
  Rand &GetRand() const { return m_rand; }
  void CauchyMutate(double mean, double variance);
  // Compares two chromosome elements. Returns -1 if the comparison is invalid
  //(unequal lengths), else returns the maximum relative pair-wise difference
  // as returned by CompareVector.
  RBTDLL_EXPORT double Compare(const ChromElement &c) const;
  // Returns true if chromosome elements have near-equal values.
  // i.e. if Compare(c) < threshold
  // Returns false if comparison is invalid (unequal lengths)
  bool Equals(const ChromElement &c, double threshold) const;
  // Convenience method that calls SetVector(const std::vector<double>& v,
  // Int& i) with i initialised to zero
  void SetVector(const std::vector<double> &v);
  // Convenience method that calls SetVector(const XOverList& v, Int& i)
  // with i initialised to zero
  void SetVector(const XOverList &v);
  // operator== and operator!= are implemented by calling Equals with the
  // static _THRESHOLD value
  RBTDLL_EXPORT friend bool operator==(const ChromElement &c1,
                                       const ChromElement &c2);
  RBTDLL_EXPORT friend bool operator!=(const ChromElement &c1,
                                       const ChromElement &c2);
  // operator<< is implemented by calling Print()
  RBTDLL_EXPORT friend std::ostream &operator<<(std::ostream &s,
                                                const ChromElement &c);

protected:
  RBTDLL_EXPORT ChromElement();
  RBTDLL_EXPORT ChromElement(const ChromElement &c);
  RBTDLL_EXPORT ChromElement &operator=(const ChromElement &c);
  // Helper method to check that index i is in range [0,v.size()}
  // and that v has sufficient elements remaining to satisfy
  // GetLength()
  bool VectorOK(const std::vector<double> &v, unsigned int i) const;
  // Helper method to check that index i is in range [0,v.size()}
  // and that v has sufficient elements remaining to satisfy
  // GetXOverLength()
  bool VectorOK(const XOverList &v, unsigned int i) const;

private:
  Rand &m_rand; // Reference to singleton random number generator
};

typedef SmartPtr<ChromElement> ChromElementPtr;
typedef std::vector<ChromElement *>
    ChromElementList; // Vector of regular pointers
typedef ChromElementList::iterator ChromElementListIter;
typedef ChromElementList::const_iterator ChromElementListConstIter;

bool operator==(const ChromElement &c1, const ChromElement &c2);
bool operator!=(const ChromElement &c1, const ChromElement &c2);
std::ostream &operator<<(std::ostream &s, const ChromElement &c);

// 2-point crossover
RBTDLL_EXPORT void Crossover(ChromElement *pChr1, ChromElement *pChr2,
                             ChromElement *pChr3, ChromElement *pChr4);

} // namespace rxdock

#endif /*RBTCHROMELEMENT_H_*/
