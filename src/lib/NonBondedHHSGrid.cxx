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

#include "NonBondedHHSGrid.h"
#include "Atom.h"
#include "FileError.h"

using namespace rxdock;

std::string NonBondedHHSGrid::_CT("NonBondedHHSGrid");

NonBondedHHSGrid::NonBondedHHSGrid(const Coord &gridMin, const Coord &gridStep,
                                   unsigned int NX, unsigned int NY,
                                   unsigned int NZ, unsigned int NPad /*=0*/)
    : BaseGrid(gridMin, gridStep, NX, NY, NZ, NPad) {
  CreateMap();
  _RBTOBJECTCOUNTER_CONSTR_("NonBondedHHSGrid");
}

NonBondedHHSGrid::NonBondedHHSGrid(std::istream &istr) : BaseGrid(istr) {
  CreateMap();
  OwnRead(istr);
  _RBTOBJECTCOUNTER_CONSTR_("NonBondedHHSGrid");
}

NonBondedHHSGrid::~NonBondedHHSGrid() {
  ClearHHSLists();
  _RBTOBJECTCOUNTER_DESTR_("NonBondedHHSGrid");
}

NonBondedHHSGrid::NonBondedHHSGrid(const NonBondedHHSGrid &aGrid)
    : BaseGrid(aGrid) {
  CreateMap();
  CopyGrid(aGrid);
  _RBTOBJECTCOUNTER_COPYCONSTR_("NonBondedHHSGrid");
}

NonBondedHHSGrid::NonBondedHHSGrid(const BaseGrid &grid) : BaseGrid(grid) {
  CreateMap();
  _RBTOBJECTCOUNTER_COPYCONSTR_("NonBondedHHSGrid");
}

NonBondedHHSGrid &NonBondedHHSGrid::operator=(const NonBondedHHSGrid &grid) {
  if (this != &grid) {
    // Clear the current atom lists
    ClearHHSLists();
    // In this case we need to explicitly call the base class operator=
    BaseGrid::operator=(grid);
    CopyGrid(grid);
  }
  return *this;
}
// Copy assignment taking base class argument
NonBondedHHSGrid &NonBondedHHSGrid::operator=(const BaseGrid &grid) {
  if (this != &grid) {
    // Clear the current atom lists
    ClearHHSLists();
    // In this case we need to explicitly call the base class operator=
    BaseGrid::operator=(grid);
  }
  return *this;
}

// Text output
void NonBondedHHSGrid::Print(std::ostream &ostr) const {
  BaseGrid::Print(ostr);
  OwnPrint(ostr);
}

// Binary output
void NonBondedHHSGrid::Write(std::ostream &ostr) const {
  BaseGrid::Write(ostr);
  OwnWrite(ostr);
}

// Binary input
void NonBondedHHSGrid::Read(std::istream &istr) {
  ClearHHSLists();
  BaseGrid::Read(istr);
  OwnRead(istr);
}

const HHS_SolvationRList &
NonBondedHHSGrid::GetHHSList(unsigned int iXYZ) const {
  if (isValid(iXYZ)) {
    return m_hhsMap[iXYZ];
  } else {
    return m_emptyList;
  }
}

const HHS_SolvationRList &NonBondedHHSGrid::GetHHSList(const Coord &c) const {
  if (isValid(c)) {
    return m_hhsMap[GetIXYZ(c)];
  } else {
    return m_emptyList;
  }
}

void NonBondedHHSGrid::SetHHSLists(HHS_Solvation *pHHS, double radius) {
  const Coord &c = (pHHS->GetAtom())->GetCoords();
  std::vector<unsigned int> sphereIndices;
  GetSphereIndices(c, radius, sphereIndices);

  for (std::vector<unsigned int>::const_iterator sphereIter =
           sphereIndices.begin();
       sphereIter != sphereIndices.end(); sphereIter++) {
    m_hhsMap[*sphereIter].push_back(pHHS);
  }
}

void NonBondedHHSGrid::ClearHHSLists() {
  // m_atomMap.clear();
  // Clear each atom list separately, without clearing the whole map (vector)
  for (HHS_SolvationListMapIter iter = m_hhsMap.begin(); iter != m_hhsMap.end();
       iter++) {
    (*iter).clear();
  }
}

void NonBondedHHSGrid::OwnPrint(std::ostream &ostr) const {
  ostr << std::endl << "Class\t" << _CT << std::endl;
  ostr << "No. of entries in the map: " << m_hhsMap.size() << std::endl;
}

void NonBondedHHSGrid::OwnWrite(std::ostream &ostr) const {}

void NonBondedHHSGrid::OwnRead(std::istream &istr) {}

void NonBondedHHSGrid::CopyGrid(const NonBondedHHSGrid &grid) {
  m_hhsMap = grid.m_hhsMap;
}

void NonBondedHHSGrid::CreateMap() { m_hhsMap = HHS_SolvationListMap(GetN()); }
