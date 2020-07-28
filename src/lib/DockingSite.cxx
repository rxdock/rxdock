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

#include "DockingSite.h"
#include "FileError.h"
#include <cstring>

using namespace rxdock;

namespace rxdock {

// Less than operator for sorting coords
class CoordCmp {
public:
  bool operator()(const Coord &c1, const Coord &c2) const {
    if (c1.xyz(0) < c2.xyz(0))
      return true;
    else if (c1.xyz(0) == c2.xyz(0)) {
      if (c1.xyz(1) < c2.xyz(1))
        return true;
      else if (c1.xyz(1) == c2.xyz(1)) {
        if (c1.xyz(2) < c2.xyz(2))
          return true;
      }
    }
    return false;
  }
};

} // namespace rxdock

// Static data members
std::string DockingSite::_CT("DockingSite");

// STL predicate for selecting atoms within a defined distance range from
// nearest cavity coords Uses precalculated distance grid
bool DockingSite::isAtomInRange::operator()(Atom *pAtom) const {
  const Coord &c = pAtom->GetCoords();
  if (!m_pGrid->isValid(c)) {
    return false;
  }
  double dist = m_pGrid->GetSmoothedValue(c);
  return (dist >= m_minDist && dist <= m_maxDist);
}

////////////////////////////////////////
// Constructors/destructors
DockingSite::DockingSite(const CavityList &cavList, double border)
    : m_cavityList(cavList), m_border(border) {
  if (!m_cavityList.empty()) {
    CoordList minCoords;
    CoordList maxCoords;
    std::transform(m_cavityList.begin(), m_cavityList.end(),
                   std::back_inserter(minCoords), ExtractCavityMinCoord);
    std::transform(m_cavityList.begin(), m_cavityList.end(),
                   std::back_inserter(maxCoords), ExtractCavityMaxCoord);
    m_minCoord = Min(minCoords);
    m_maxCoord = Max(maxCoords);
    _RBTOBJECTCOUNTER_CONSTR_(_CT);
  }

  // Make sure grid has been calculated
  if (m_spGrid.Null()) {
    CreateGrid();
  }
}

DockingSite::DockingSite(json j) {
  j.get_to(*this);
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

DockingSite::~DockingSite() { _RBTOBJECTCOUNTER_DESTR_(_CT); }

// Insertion operator
std::ostream &rxdock::operator<<(std::ostream &s, const DockingSite &site) {
  site.Print(s);
  return s;
}

// Virtual function for dumping docking site details to an output stream
// Derived classes can override if required
void DockingSite::Print(std::ostream &s) const {
  s << "Total volume " << GetVolume() << " A^3" << std::endl;
  for (unsigned int i = 0; i < m_cavityList.size(); i++) {
    s << "Cavity #" << i + 1 << "\t" << *(m_cavityList[i]) << std::endl;
  }
}

// Return distance grid, calculated on demand if necessary
RealGridPtr DockingSite::GetGrid() {
  if (m_spGrid.Null()) {
    CreateGrid();
  }
  return m_spGrid;
}

// returns total volume of all cavities in A^3
double DockingSite::GetVolume() const {
  double vol(0.0);
  for (auto &iter : m_cavityList) {
    vol += iter->GetVolume();
  }
  return vol;
}

// Returns the combined coord lists of all the cavities
void DockingSite::GetCoordList(CoordList &retVal) const {
  retVal.clear();
  for (auto &iter : m_cavityList) {
    const CoordList &cavCoords = iter->GetCoordList();
    retVal.reserve(retVal.size() + cavCoords.size());
    std::copy(cavCoords.begin(), cavCoords.end(), std::back_inserter(retVal));
    // Sort the coords so we can remove any dups
    std::sort(retVal.begin(), retVal.end(), CoordCmp());
    CoordListIter uniqIter = std::unique(retVal.begin(), retVal.end());
    retVal.erase(uniqIter, retVal.end());
    // std::cout << "Cav = " << cavCoords.size() << "; total = " <<
    // retVal.size() << std::endl;
  }
}

// Filters an atom list according to distance from the cavity coords
// Only returns atoms within minDist and maxDist from cavity
AtomList DockingSite::GetAtomList(const AtomList &atomList, double minDist,
                                  double maxDist) {
  if (maxDist > m_border) {
    throw BadArgument(_WHERE_,
                      "maxDist is greater than grid border; recalculate grid");
  }
  AtomList newAtomList;
  if (m_spGrid.Null()) {
    CreateGrid();
  }
  if (m_spGrid.Null())
    return newAtomList;
  isAtomInRange bInRange(m_spGrid, minDist, maxDist);
  std::copy_if(atomList.begin(), atomList.end(),
               std::back_inserter(newAtomList), bInRange);
  return newAtomList;
}

// Filters an atom list according to distance from the cavity coords
// Only returns atoms within maxDist from cavity
// This version does not require the cavity grid
AtomList DockingSite::GetAtomList(const AtomList &atomList, double maxDist) {
  AtomList newAtomList;
  CoordList allCoords;
  GetCoordList(allCoords);
  isAtomNearCoordList bInRange(allCoords, maxDist);
  std::copy_if(atomList.begin(), atomList.end(),
               std::back_inserter(newAtomList), bInRange);
  return newAtomList;
}

// Returns the count of atoms within minDist and maxDist from cavity
unsigned int DockingSite::GetNumAtoms(const AtomList &atomList, double minDist,
                                      double maxDist) {
  if (maxDist > m_border) {
    throw BadArgument(_WHERE_,
                      "maxDist is greater than grid border; recalculate grid");
  }
  if (m_spGrid.Null()) {
    CreateGrid();
  }
  if (m_spGrid.Null())
    return 0;
  isAtomInRange bInRange(m_spGrid, minDist, maxDist);
  return std::count_if(atomList.begin(), atomList.end(), bInRange);
}

// Create a grid whose values represent the distance to the nearest cavity coord
void DockingSite::CreateGrid() {
  m_spGrid = RealGridPtr();
  if (m_cavityList.empty())
    return;

  // Get the grid step to use from the first cavity
  Vector gridStep = m_cavityList.front()->GetGridStep();
  Coord minCoord = m_minCoord - m_border;
  Coord maxCoord = m_maxCoord + m_border;
  Vector extent = maxCoord - minCoord;
  Eigen::Vector3d nXYZ = extent.xyz.array() / gridStep.xyz.array();
  unsigned int nX = static_cast<unsigned int>(nXYZ(0)) + 1;
  unsigned int nY = static_cast<unsigned int>(nXYZ(1)) + 1;
  unsigned int nZ = static_cast<unsigned int>(nXYZ(2)) + 1;
  m_spGrid = RealGridPtr(new RealGrid(minCoord, gridStep, nX, nY, nZ));
  m_spGrid->SetAllValues(999999.9);

  // Get the total list of cavity coords
  CoordList allCoords;
  for (auto &iter : m_cavityList) {
    const CoordList cavCoords = iter->GetCoordList();
    // Reserve enough space for appending the next cavity coord list
    allCoords.reserve(allCoords.size() + cavCoords.size());
    // Perform the append
    std::copy(cavCoords.begin(), cavCoords.end(),
              std::back_inserter(allCoords));
    // We can save a bit of time by initialising the distance^2 of the grid
    // points nearest to each cavity coord Normally zero, but we allow the
    // possibility that the grid step may be different i.e. the cavity coords
    // may not lie directly on a grid point.
    for (auto &cIter : cavCoords) {
      unsigned int i =
          m_spGrid->GetIXYZ(cIter); // Grid index of nearest grid point
      double dist2 = Length2(cIter, m_spGrid->GetCoord(i));
      m_spGrid->SetValue(i, dist2);
    }
    // Sort the coords so we can remove any dups

    std::sort(allCoords.begin(), allCoords.end(), CoordCmp());
    CoordListIter uniqIter = std::unique(allCoords.begin(), allCoords.end());
    allCoords.erase(uniqIter, allCoords.end());
    std::cout << "Cav = " << cavCoords.size()
              << "; total = " << allCoords.size() << std::endl;
  }

  // Loop over all grid points in the distance grid
  // Can terminate when distance^2 is less than or equal to mindist^2 (shortest
  // length of grid interval)
  double mindist2 = std::pow(gridStep.xyz.minCoeff(), 2);

  // parallelize over grid cuboids for multicore handling of time critical part
#pragma omp parallel for schedule(dynamic, 1)
  for (int i = 0; i < m_spGrid->GetN(); i++) {
    const Coord &c = m_spGrid->GetCoord(i);
    // Initialise dist^2 from initial grid value (999999.9 or 0.0 for cavity
    // coords)
    double dist2 = m_spGrid->GetValue(i);
    // Determine min distance to any of the cavity coords
    for (CoordListConstIter iter = allCoords.begin();
         iter != allCoords.end() && dist2 > mindist2; iter++) {
      dist2 = std::min(dist2, Length2(c - (*iter)));
    }
    m_spGrid->SetValue(i, std::sqrt(dist2));
  }
}

void rxdock::to_json(json &j, const DockingSite &site) {
  json cavityList;
  for (const auto &cIter : site.m_cavityList) {
    json cavity = *cIter;
    cavityList.push_back(cavity);
  }

  j = json{{"min-coord", site.m_minCoord},
           {"max-coord", site.m_maxCoord},
           {"border", site.m_border},
           {"cavities", cavityList}};
  if (cavityList.size() > 0 && !site.m_spGrid.Null()) {
    j["real-grid"] = *site.m_spGrid;
  }
}

void rxdock::from_json(const json &j, DockingSite &site) {
  site.m_cavityList.clear();
  site.m_cavityList.reserve(j.at("cavities").size());
  for (auto &cavity : j.at("cavities")) {
    CavityPtr spCavity = CavityPtr(new Cavity(cavity));
    site.m_cavityList.push_back(spCavity);
  }

  j.at("min-coord").get_to(site.m_minCoord);
  j.at("max-coord").get_to(site.m_maxCoord);
  if (site.m_cavityList.size() > 0) {
    site.m_spGrid = RealGridPtr(new RealGrid(j.at("real-grid")));
  }
  j.at("border").get_to(site.m_border);
}
