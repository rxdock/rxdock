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

#include "rxdock/NonBondedGrid.h"
#include "rxdock/FileError.h"

#include <loguru.hpp>

using namespace rxdock;

// Static data members
std::string NonBondedGrid::_CT("NonBondedGrid");

////////////////////////////////////////
// Constructors/destructors
// Construct a NXxNYxNZ grid running from gridMin at gridStep resolution
NonBondedGrid::NonBondedGrid(const Coord &gridMin, const Coord &gridStep,
                             unsigned int NX, unsigned int NY, unsigned int NZ,
                             unsigned int NPad)
    : BaseGrid(gridMin, gridStep, NX, NY, NZ, NPad) {
  CreateMap();
  _RBTOBJECTCOUNTER_CONSTR_("NonBondedGrid");
}

// Constructor reading params from binary stream
NonBondedGrid::NonBondedGrid(json j) : BaseGrid(j) {
  CreateMap();
  j.get_to(*this);
  _RBTOBJECTCOUNTER_CONSTR_("NonBondedGrid");
}

// Default destructor
NonBondedGrid::~NonBondedGrid() {
  ClearAtomLists();
  _RBTOBJECTCOUNTER_DESTR_("NonBondedGrid");
}

// Copy constructor
NonBondedGrid::NonBondedGrid(const NonBondedGrid &grid) : BaseGrid(grid) {
  CreateMap();
  CopyGrid(grid);
  _RBTOBJECTCOUNTER_COPYCONSTR_("NonBondedGrid");
}

// Copy constructor taking base class argument
NonBondedGrid::NonBondedGrid(const BaseGrid &grid) : BaseGrid(grid) {
  CreateMap();
  _RBTOBJECTCOUNTER_COPYCONSTR_("NonBondedGrid");
}

// Copy assignment
NonBondedGrid &NonBondedGrid::operator=(const NonBondedGrid &grid) {
  if (this != &grid) {
    // Clear the current atom lists
    ClearAtomLists();
    // In this case we need to explicitly call the base class operator=
    BaseGrid::operator=(grid);
    CopyGrid(grid);
  }
  return *this;
}
// Copy assignment taking base class argument
NonBondedGrid &NonBondedGrid::operator=(const BaseGrid &grid) {
  if (this != &grid) {
    // Clear the current atom lists
    ClearAtomLists();
    // In this case we need to explicitly call the base class operator=
    BaseGrid::operator=(grid);
  }
  return *this;
}

////////////////////////////////////////
// Virtual functions for reading/writing grid data to streams in
// text and binary format
// Subclasses should provide their own private OwnPrint
// method to handle subclass data members, and override the public
// Print method

// Text output
void NonBondedGrid::Print(std::ostream &ostr) const {
  BaseGrid::Print(ostr);
  OwnPrint(ostr);
}

////////////////////////////////////////
// Public methods
////////////////

/////////////////////////
// Get attribute functions
/////////////////////////

const AtomRList &NonBondedGrid::GetAtomList(unsigned int iXYZ) const {
  // AtomListMapConstIter iter = m_atomMap.find(iXYZ);
  // if (iter != m_atomMap.end())
  //  return iter->second;
  // else
  //  return m_emptyList;

  // DM 6 Nov 2000 - map replaced by vector
  if (isValid(iXYZ)) {
    return m_atomMap[iXYZ];
  } else {
    return m_emptyList;
  }
}

const AtomRList &NonBondedGrid::GetAtomList(const Coord &c) const {
  // if (isValid(c))
  //  return GetAtomList(GetIXYZ(c));
  // else
  //  return m_emptyList;

  // DM 6 Nov 2000 - map replaced by vector
  if (isValid(c)) {
    return m_atomMap[GetIXYZ(c)];
  } else {
    return m_emptyList;
  }
}

/////////////////////////
// Set attribute functions
/////////////////////////
void NonBondedGrid::SetAtomLists(Atom *pAtom, double radius) {
  const Coord &c = pAtom->GetCoords();
  std::vector<unsigned int> sphereIndices;
  GetSphereIndices(c, radius, sphereIndices);

  for (std::vector<unsigned int>::const_iterator sphereIter =
           sphereIndices.begin();
       sphereIter != sphereIndices.end(); sphereIter++) {
    // AtomListMapIter mapIter = m_atomMap.find(*sphereIter);
    // if (mapIter != m_atomMap.end()) {
    //  (mapIter->second).push_back(pAtom);
    //}
    // else {
    //  m_atomMap[*sphereIter] = AtomRList(1,pAtom);
    //}
    m_atomMap[*sphereIter].push_back(pAtom);
  }
}

void NonBondedGrid::ClearAtomLists() {
  // m_atomMap.clear();
  // Clear each atom list separately, without clearing the whole map (vector)
  for (AtomListMapIter iter = m_atomMap.begin(); iter != m_atomMap.end();
       iter++) {
    (*iter).clear();
  }
}

void NonBondedGrid::UniqueAtomLists() {
  for (AtomListMapIter iter = m_atomMap.begin(); iter != m_atomMap.end();
       iter++) {
    LOG_F(1, "NonBondedGrid::UniqueAtomLists: before = {}", (*iter).size());
    std::sort((*iter).begin(), (*iter).end(), AtomPtrCmp_Ptr());
    AtomRListIter uniqIter = std::unique((*iter).begin(), (*iter).end());
    (*iter).erase(uniqIter, (*iter).end());
    LOG_F(1, "NonBondedGrid::UniqueAtomLists: after = {}", (*iter).size());
  }
}

///////////////////////////////////////////////////////////////////////////
// Protected methods

// Protected method for writing data members for this class to text stream
void NonBondedGrid::OwnPrint(std::ostream &ostr) const {
  ostr << std::endl << "Class\t" << _CT << std::endl;
  ostr << "No. of entries in the map: " << m_atomMap.size() << std::endl;
  // TO BE COMPLETED - no real need for dumping the atom list info
}

///////////////////////////////////////////////////////////////////////////
// Private methods
// Helper function called by copy constructor and assignment operator
void NonBondedGrid::CopyGrid(const NonBondedGrid &grid) {
  // This copies the atom lists, but of course the atoms themselves are not
  // copied
  m_atomMap = grid.m_atomMap;
}

// DM 6 Nov 2000 - create AtomListMap of the appropriate size
void NonBondedGrid::CreateMap() { m_atomMap = AtomListMap(GetN()); }
