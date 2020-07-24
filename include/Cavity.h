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

// Struct for holding a coordinate list representing an active site cavity
// Properties are read-only once set in the constructor

#ifndef _RBTCAVITY_H_
#define _RBTCAVITY_H_

#include "Config.h"
#include "Coord.h"
#include "PrincipalAxes.h"
#include "RealGrid.h"

namespace rxdock {

class Cavity {
public:
  ////////////////////////////////////////
  // Constructors/destructors

  // Constructor
  Cavity(const CoordList &coordList, const Vector gridStep)
      : m_coordList(coordList), m_gridStep(gridStep),
        m_prAxes(GetPrincipalAxesOfAtoms(coordList)),
        m_minCoord(Min(coordList)), m_maxCoord(Max(coordList)) {}

  // Constructor reading from binary stream
  Cavity(std::istream &istr) { Read(istr); }

  // Destructor
  virtual ~Cavity() {}

  ///////////////////////////////////////////////
  // Stream functions
  ///////////////////////////////////////////////

  // Insertion operator
  friend std::ostream &operator<<(std::ostream &s, const Cavity &cavity) {
    cavity.Print(s);
    return s;
  }

  // Virtual function for dumping cavity details to an output stream
  // Derived classes can override if required
  virtual void Print(std::ostream &s) const {
    s << "Size=" << m_coordList.size() << " points; Vol=" << GetVolume()
      << " A^3; Min=" << m_minCoord << "; Max=" << m_maxCoord
      << "; Center=" << m_prAxes.com << "; Extent=" << m_maxCoord - m_minCoord;
  }

  // Binary output (serialisation)
  virtual void Write(std::ostream &ostr) const {
    // DM 4 Apr 2002 - Write grid step
    m_gridStep.Write(ostr);
    // Write number of coords
    int nCoords = m_coordList.size();
    WriteWithThrow(ostr, (const char *)&nCoords, sizeof(nCoords));
    for (CoordListConstIter cIter = m_coordList.begin();
         cIter != m_coordList.end(); cIter++) {
      (*cIter).Write(ostr); // Write each coord
    }
  }

  // Binary input, replaces existing cavity
  virtual void Read(std::istream &istr) {
    // Clear the existing cavity
    m_coordList.clear();
    // DM 4 Apr 2002 - Read grid step
    m_gridStep.Read(istr);
    // Read number of coords
    int nCoords;
    ReadWithThrow(istr, (char *)&nCoords, sizeof(nCoords));
    m_coordList.reserve(nCoords);
    Coord c;
    // Read each coord and add it to the cavity list
    for (int i = 0; i < nCoords; i++) {
      c.Read(istr);
      m_coordList.push_back(c);
    }
    // Recalculate the other properties
    m_prAxes = GetPrincipalAxesOfAtoms(m_coordList);
    m_minCoord = Min(m_coordList);
    m_maxCoord = Max(m_coordList);
  }

  // DM 4 Apr 2002 - return a grid with all cavity points set to 1.0
  RealGridPtr GetGrid() const {
    Vector extent = m_maxCoord - m_minCoord;
    // Add +3 to allow a small border
    Eigen::Vector3d nXYZ = extent.xyz.array() / m_gridStep.xyz.array();
    unsigned int nX = static_cast<unsigned int>(nXYZ(0)) + 3;
    unsigned int nY = static_cast<unsigned int>(nXYZ(1)) + 3;
    unsigned int nZ = static_cast<unsigned int>(nXYZ(2)) + 3;
    RealGridPtr spGrid = RealGridPtr(
        new RealGrid(m_minCoord - m_gridStep, m_gridStep, nX, nY, nZ));
    for (CoordListConstIter iter = m_coordList.begin();
         iter != m_coordList.end(); iter++) {
      spGrid->SetValue(*iter, 1.0);
    }
    return spGrid;
  }

  ////////////////////////////////////////
  // Public methods
  ////////////////
  const Coord &GetCenterOfMass() const { return m_prAxes.com; }
  const PrincipalAxes &GetPrincipalAxes() const { return m_prAxes; }
  int GetNumCoords() const { return m_coordList.size(); }
  const CoordList &GetCoordList() const { return m_coordList; }
  const Coord &GetMinCoord() const { return m_minCoord; }
  const Coord &GetMaxCoord() const { return m_maxCoord; }
  const Vector &GetGridStep() const { return m_gridStep; }
  double GetVolume() const {
    return m_coordList.size() * m_gridStep.xyz.prod();
  }

protected:
  ////////////////////////////////////////
  // Protected methods
  ///////////////////

private:
  ////////////////////////////////////////
  // Private methods
  /////////////////
  Cavity();                          // Disable default constructor
  Cavity(const Cavity &);            // Copy constructor disabled by default
  Cavity &operator=(const Cavity &); // Copy assignment disabled by
                                     // default

protected:
  ////////////////////////////////////////
  // Protected data
  ////////////////

private:
  ////////////////////////////////////////
  // Private data
  //////////////
  CoordList m_coordList;
  Vector m_gridStep;
  PrincipalAxes m_prAxes;
  Coord m_minCoord;
  Coord m_maxCoord;
};

// Useful typedefs
typedef SmartPtr<Cavity> CavityPtr;        // Smart pointer
typedef std::vector<CavityPtr> CavityList; // Vector of smart pointers
typedef CavityList::iterator CavityListIter;
typedef CavityList::const_iterator CavityListConstIter;

////////////////////////////////////////////////////////
// Comparison functions for sorting CavityPtr containers
// For use by STL sort algorithms
////////////////////////////////////////////////////////
// Less than operator for sorting CavityPtrs into ascending order by center
// of mass distance from a given coord
class CavityPtrCmp_Distance {
  const Coord &c;

public:
  CavityPtrCmp_Distance(const Coord &cc) : c(cc) {}
  bool operator()(CavityPtr spCav1, CavityPtr spCav2) const {
    return Length2(spCav1->GetCenterOfMass() - c) <
           Length2(spCav2->GetCenterOfMass() - c);
  }
};
// Less than operator for sorting CavityPtrs into descending order by volume
class CavityPtrCmp_Volume {
public:
  CavityPtrCmp_Volume() {}
  bool operator()(CavityPtr spCav1, CavityPtr spCav2) const {
    return spCav1->GetNumCoords() > spCav2->GetNumCoords();
  }
};
////////////////////////////////////////////////////////
// Predicate functions for CavityPtr
// For use by STL algorithms
////////////////////////////////////////////////////////
class isCavityNearCoord : public std::function<bool(CavityPtr)> {
  const Coord &c;
  double r2; // radius squared (to avoid taking square roots)
public:
  explicit isCavityNearCoord(const Coord &cc, double rr) : c(cc), r2(rr * rr) {}
  bool operator()(CavityPtr spCavity) const {
    return Length2(spCavity->GetCenterOfMass() - c) <= r2;
  }
};

////////////////////////////////////////////
// Functions objects for performing actions on cavities
// For use by STL algorithms
////////////////////////////////////////////
inline const Coord &ExtractCavityMinCoord(CavityPtr spCav) {
  return spCav->GetMinCoord();
}
inline const Coord &ExtractCavityMaxCoord(CavityPtr spCav) {
  return spCav->GetMaxCoord();
}

} // namespace rxdock

#endif //_RBTCAVITY_H_
