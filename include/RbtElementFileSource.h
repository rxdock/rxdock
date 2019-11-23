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

// File reader for Rbt elemental data file (atomic no,element name,vdw radii)
// Provides the data as a vector of structs

#ifndef _RBTELEMENTFILESOURCE_H_
#define _RBTELEMENTFILESOURCE_H_

#include "RbtBaseFileSource.h"

// Simple struct for holding the element data
class RbtElementData {
public:
  RbtElementData()
      : atomicNo(0), element(""), minVal(0), maxVal(0), commonVal(0), mass(0.0),
        vdwRadius(0.0) {}
  int atomicNo;
  std::string element;
  int minVal;
  int maxVal;
  int commonVal;
  double mass;
  double vdwRadius; // Regular vdW radius
};

// Map with element data indexed by element name
typedef std::map<std::string, RbtElementData> RbtStringElementDataMap;
typedef RbtStringElementDataMap::iterator RbtStringElementDataMapIter;
typedef RbtStringElementDataMap::const_iterator
    RbtStringElementDataMapConstIter;
// Map with element data indexed by atomic number
typedef std::map<int, RbtElementData> RbtIntElementDataMap;
typedef RbtIntElementDataMap::iterator RbtIntElementDataMapIter;
typedef RbtIntElementDataMap::const_iterator RbtIntElementDataMapConstIter;

class RbtElementFileSource : public RbtBaseFileSource {
public:
  // Constructors
  // RbtElementFileSource(const char* fileName);
  RbtElementFileSource(const std::string &fileName);

  // Destructor
  virtual ~RbtElementFileSource();

  ////////////////////////////////////////
  // Public methods
  RBTDLL_EXPORT std::string GetTitle();
  std::string GetVersion();
  unsigned int GetNumElements();
  std::vector<std::string> GetElementNameList(); // List of element names
  std::vector<int> GetAtomicNumberList();        // List of atomic numbers
  // Get element data for a given element name, throws error if not found
  RbtElementData GetElementData(const std::string &strElementName);
  // Get element data for a given atomic number, throws error if not found
  RbtElementData GetElementData(int nAtomicNumber);
  // Check if given element name is present
  bool isElementNamePresent(const std::string &strElementName);
  // Check if given atomic number is present
  bool isAtomicNumberPresent(int nAtomicNumber);

  // Get vdw radius increment for hydrogen-bonding donors
  double GetHBondRadiusIncr();
  // Get vdw radius increment for atoms with implicit hydrogens
  double GetImplicitRadiusIncr();

protected:
  // Protected methods

private:
  // Private methods
  RbtElementFileSource(); // Disable default constructor
  RbtElementFileSource(
      const RbtElementFileSource &); // Copy constructor disabled by default
  RbtElementFileSource &operator=(
      const RbtElementFileSource &); // Copy assignment disabled by default

  // Pure virtual in RbtBaseFileSource - needs to be defined here
  virtual void Parse();
  void ClearElementDataCache();

protected:
  // Protected data

private:
  // Private data
  std::string m_strTitle;
  std::string m_strVersion;
  double m_dHBondRadiusIncr;    // Increment to add to vdW radius for H-bonding
                                // hydrogens
  double m_dImplicitRadiusIncr; // Increment to add to vdW radius for atoms
                                // with implicit hydrogens
  RbtStringElementDataMap m_elementNameMap;
  RbtIntElementDataMap m_atomicNumberMap;
};

// useful typedefs
typedef SmartPtr<RbtElementFileSource> RbtElementFileSourcePtr; // Smart pointer

#endif //_RBTELEMENTFILESOURCE_H_
