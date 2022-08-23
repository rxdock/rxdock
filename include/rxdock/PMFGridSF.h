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

#ifndef _RBTPMFGRIDSF_H_
#define _RBTPMFGRIDSF_H_

#include "rxdock/BaseInterSF.h"
#include "rxdock/RealGrid.h"

#include <nlohmann/json.hpp>

using json = nlohmann::json;

namespace rxdock {

class PMFGridSF : public BaseInterSF {
  bool m_bSmoothed;
  AtomList theLigandList;            // vector to store the ligand
  std::vector<PMFType> theTypeList;  // store PMF used types here
  std::vector<RealGridPtr> theGrids; // grids with PMF data

  void ReadGrids(json pmfGrids);

public:
  static const std::string _CT;   // class name
  static const std::string _GRID; // filename extension (.grd)
  static const std::string
      _SMOOTHED; // controls wether to smooth the grid values

  PMFGridSF(const std::string &strName = "pmf-grid");
  virtual ~PMFGridSF();

protected:
  virtual void SetupReceptor();
  virtual void SetupLigand();
  virtual void SetupScore() {}
  virtual double RawScore() const;
  unsigned int GetCorrectedType(PMFType aType) const;

private:
  friend void to_json(json &j, const PMFGridSF &sf) {
    json atomList;
    for (const auto &cIter : sf.theLigandList) {
      json atom = *cIter;
      atomList.push_back(atom);
    }
    json grids;
    for (const auto &cIter : sf.theGrids) {
      json grid = *cIter;
      grids.push_back(grid);
    }

    j = json{{"smoothed", sf.m_bSmoothed},
             {"atoms", atomList},
             {"pmf-types", sf.theTypeList},
             {"pmf-data", grids}};
  };
  friend void from_json(const json &j, PMFGridSF &sf) {
    j.at("smoothed").get_to(sf.m_bSmoothed);

    sf.theLigandList.clear();
    sf.theLigandList.reserve(j.at("atoms").size());
    for (auto &atom : j.at("atoms")) {
      AtomPtr spAtom = AtomPtr(new Atom(atom));
      sf.theLigandList.push_back(spAtom);
    }
    j.at("pmf-types").get_to(sf.theTypeList);

    sf.theGrids.clear();
    sf.theGrids.reserve(j.at("grids").size());
    for (auto &grid : j.at("grids")) {
      RealGridPtr spGrid = RealGridPtr(new RealGrid(grid));
      sf.theGrids.push_back(spGrid);
    }
  };
};

} // namespace rxdock

#endif // _RBTPMFGRIDSF_H_
