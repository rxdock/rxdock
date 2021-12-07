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

#include "rxdock/VdwGridSF.h"
#include "rxdock/FileError.h"
#include "rxdock/WorkSpace.h"

#include <loguru.hpp>

using namespace rxdock;

// Static data members
const std::string VdwGridSF::_CT = "VdwGridSF";
const std::string VdwGridSF::_GRID = "grid";
const std::string VdwGridSF::_SMOOTHED = "smoothed";

const std::string &VdwGridSF::GetCt() { return _CT; }

// NB - Virtual base class constructor (BaseSF) gets called first,
// implicit constructor for BaseInterSF is called second
VdwGridSF::VdwGridSF(const std::string &strName)
    : BaseSF(_CT, strName), m_bSmoothed(true) {
  LOG_F(2, "VdwGridSF parameterised constructor");
  // Add parameters
  AddParameter(_GRID, ".grd");
  AddParameter(_SMOOTHED, m_bSmoothed);
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

VdwGridSF::~VdwGridSF() {
  LOG_F(2, "VdwGridSF parameterised constructor");
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

void VdwGridSF::SetupReceptor() {
  m_grids.clear();
  if (GetReceptor().Null())
    return;

  // Trap multiple receptor conformations and flexible OH/NH3 here: this SF does
  // not support them yet
  bool bEnsemble = (GetReceptor()->GetNumSavedCoords() > 1);
  bool bFlexRec = GetReceptor()->isFlexible();
  if (bEnsemble || bFlexRec) {
    std::string message("Vdw grid scoring function does not support multiple "
                        "receptor conformations\n");
    message += "or flexible OH/NH3 groups yet";
    throw InvalidRequest(_WHERE_, message);
  }

  // Read grids
  // File names are composed of workspace name + grid suffix
  // Reasonably safe to assume that GetWorkSpace is not nullptr,
  // as otherwise SetupReceptor would not have been called
  std::string strWSName = GetWorkSpace()->GetName();

  std::string strSuffix = GetParameter(_GRID);
  std::string strFile = GetDataFileName("data/grids", strWSName + strSuffix);
  std::ifstream file(strFile.c_str());
  json vdwGrids;
  file >> vdwGrids;
  file.close();
  ReadGrids(vdwGrids.at("vdw-grids"));
}

void VdwGridSF::SetupLigand() {
  m_ligAtomList.clear();
  m_ligAtomTypes.clear();
  if (GetLigand().Null())
    return;

  AtomList tmpList = GetLigand()->GetAtomList();
  // Strip off the smart pointers
  std::copy(tmpList.begin(), tmpList.end(), std::back_inserter(m_ligAtomList));
}

void VdwGridSF::SetupSolvent() {
  ModelList solvent = GetSolvent();
  if (!solvent.empty()) {
    std::string message(
        "Vdw grid scoring function does not support explicit solvent yet\n");
    throw InvalidRequest(_WHERE_, message);
  }
}

void VdwGridSF::SetupScore() {
  // Determine probe grid type for each atom, based on comparing Tripos atom
  // type with probe atom types This needs to be in SetupScore as it is
  // dependent on both the ligand and receptor grid data
  m_ligAtomTypes.clear();
  if (m_ligAtomList.empty())
    return;

  m_ligAtomTypes.reserve(m_ligAtomList.size());
  // Check if we have a grid for the UNDEFINED type:
  // If so, we can use it if a particular atom type grid is missing
  // If not, then we have to throw an error if a particular atom type grid is
  // missing
  bool bHasUndefined = !m_grids[TriposAtomType::UNDEFINED].Null();
  TriposAtomType triposType;

  if (bHasUndefined) {
    LOG_F(INFO, "The grid for the UNDEFINED atom type will be used if a grid "
                "is missing for any atom type");
  } else {
    LOG_F(INFO, "There is no grid for the UNDEFINED atom type. An error will "
                "be thrown if a grid is missing for any atom type");
  }

  for (AtomRListConstIter iter = m_ligAtomList.begin();
       iter != m_ligAtomList.end(); iter++) {
    TriposAtomType::eType aType = (*iter)->GetTriposType();

    // If there is no grid for this atom type, revert to using the UNDEFINED
    // grid if available else throw an error
    if (m_grids[aType].Null()) {
      std::string strError = "No vdw grid available for " +
                             (*iter)->GetFullAtomName() + " (type " +
                             triposType.Type2Str(aType) + ")";
      LOG_F(ERROR, "{}", strError);
      if (bHasUndefined) {
        aType = TriposAtomType::UNDEFINED;
      } else {
        throw FileError(_WHERE_, strError);
      }
    }

    m_ligAtomTypes.push_back(aType);
    LOG_F(1, "Using grid #{} for {}", aType, (*iter)->GetFullAtomName());
  }
}

double VdwGridSF::RawScore() const {
  double score = 0.0;

  // Check grids are defined
  if (m_grids.empty())
    return score;

  // Loop over all ligand atoms
  AtomRListConstIter aIter = m_ligAtomList.begin();
  TriposAtomTypeListConstIter tIter = m_ligAtomTypes.begin();
  if (m_bSmoothed) {
    for (; aIter != m_ligAtomList.end(); aIter++, tIter++) {
      score += m_grids[*tIter]->GetSmoothedValue((*aIter)->GetCoords());
    }
  } else {
    for (; aIter != m_ligAtomList.end(); aIter++, tIter++) {
      score += m_grids[*tIter]->GetValue((*aIter)->GetCoords());
    }
  }
  return score;
}

// Read grids from input stream, checking that header string matches
// VdwGridSF
void VdwGridSF::ReadGrids(json vdwGrids) {
  m_grids.clear();

  // Now read number of grids
  LOG_F(1, "VdwGridSF: reading {} grids...", vdwGrids.size());

  // We actually create a vector of size TriposAtomType::MAXTYPES
  // Each grid is in the file is prefixed by the atom type string (e.g. C.2)
  // This string is converted to the corresponding TriposAtomType::eType
  // and the grid stored at m_grids[eType]
  // This is slightly more transferable as it means grids can be read correctly
  // even if the atom type enums change (as long as the atom type string stay
  // the same) It also means we do not have to have a grid for each and every
  // atom type if we don't want to
  m_grids = RealGridList(TriposAtomType::MAXTYPES);
  for (int i = 0; i < vdwGrids.size(); i++) {
    // Read the atom type string
    std::string strType;
    vdwGrids.at(i).at("tripos-type").get_to(strType);
    TriposAtomType triposType;
    TriposAtomType::eType aType = triposType.Str2Type(strType);
    LOG_F(1, "Grid# {} atom type={} (type #{})", i, strType, aType);
    // Now we can read the grid
    RealGridPtr spGrid(new RealGrid(vdwGrids.at(i).at("real-grid")));
    m_grids[aType] = spGrid;
  }
}

// DM 25 Oct 2000 - track changes to parameter values in local data members
// ParameterUpdated is invoked by ParamHandler::SetParameter
void VdwGridSF::ParameterUpdated(const std::string &strName) {
  // DM 25 Oct 2000 - heavily used params
  if (strName == _SMOOTHED) {
    m_bSmoothed = GetParameter(_SMOOTHED);
  } else {
    BaseSF::ParameterUpdated(strName);
  }
}
