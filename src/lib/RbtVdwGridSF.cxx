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

#include "RbtVdwGridSF.h"
#include "RbtFileError.h"
#include "RbtWorkSpace.h"

// Static data members
std::string RbtVdwGridSF::_CT("RbtVdwGridSF");
std::string RbtVdwGridSF::_GRID("GRID");
std::string RbtVdwGridSF::_SMOOTHED("SMOOTHED");

// NB - Virtual base class constructor (RbtBaseSF) gets called first,
// implicit constructor for RbtBaseInterSF is called second
RbtVdwGridSF::RbtVdwGridSF(const std::string &strName)
    : RbtBaseSF(_CT, strName), m_bSmoothed(true) {
  // Add parameters
  AddParameter(_GRID, ".grd");
  AddParameter(_SMOOTHED, m_bSmoothed);
#ifdef _DEBUG
  std::cout << _CT << " parameterised constructor" << std::endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtVdwGridSF::~RbtVdwGridSF() {
#ifdef _DEBUG
  std::cout << _CT << " destructor" << std::endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

void RbtVdwGridSF::SetupReceptor() {
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
    throw RbtInvalidRequest(_WHERE_, message);
  }

  // Read grids
  // File names are composed of workspace name + grid suffix
  // Reasonably safe to assume that GetWorkSpace is not NULL,
  // as otherwise SetupReceptor would not have been called
  std::string strWSName = GetWorkSpace()->GetName();

  std::string strSuffix = GetParameter(_GRID);
  std::string strFile =
      Rbt::GetRbtFileName("data/grids", strWSName + strSuffix);
  // DM 26 Sep 2000 - std::ios_base::binary qualifier doesn't appear to be valid
  // with IRIX CC
#ifdef __sgi
  std::ifstream istr(strFile.c_str(), std::ios_base::in);
#else
  std::ifstream istr(strFile.c_str(),
                     std::ios_base::in | std::ios_base::binary);
#endif
  ReadGrids(istr);
  istr.close();
}

void RbtVdwGridSF::SetupLigand() {
  m_ligAtomList.clear();
  m_ligAtomTypes.clear();
  if (GetLigand().Null())
    return;

  RbtAtomList tmpList = GetLigand()->GetAtomList();
  // Strip off the smart pointers
  std::copy(tmpList.begin(), tmpList.end(), std::back_inserter(m_ligAtomList));
}

void RbtVdwGridSF::SetupSolvent() {
  RbtModelList solvent = GetSolvent();
  if (!solvent.empty()) {
    std::string message(
        "Vdw grid scoring function does not support explicit solvent yet\n");
    throw RbtInvalidRequest(_WHERE_, message);
  }
}

void RbtVdwGridSF::SetupScore() {
  // Determine probe grid type for each atom, based on comparing Tripos atom
  // type with probe atom types This needs to be in SetupScore as it is
  // dependent on both the ligand and receptor grid data
  m_ligAtomTypes.clear();
  if (m_ligAtomList.empty())
    return;

  int iTrace = GetTrace();
  m_ligAtomTypes.reserve(m_ligAtomList.size());
  // Check if we have a grid for the UNDEFINED type:
  // If so, we can use it if a particular atom type grid is missing
  // If not, then we have to throw an error if a particular atom type grid is
  // missing
  bool bHasUndefined = !m_grids[RbtTriposAtomType::UNDEFINED].Null();
  RbtTriposAtomType triposType;

  if (iTrace > 1) {
    if (bHasUndefined) {
      std::cout
          << "The grid for the UNDEFINED atom type will be used if a grid is "
             "missing for any atom type"
          << std::endl;
    } else {
      std::cout
          << "There is no grid for the UNDEFINED atom type. An error will be "
             "thrown if a grid is missing for any atom type"
          << std::endl;
    }
  }

  for (RbtAtomRListConstIter iter = m_ligAtomList.begin();
       iter != m_ligAtomList.end(); iter++) {
    RbtTriposAtomType::eType aType = (*iter)->GetTriposType();

    // If there is no grid for this atom type, revert to using the UNDEFINED
    // grid if available else throw an error
    if (m_grids[aType].Null()) {
      std::string strError = "No vdw grid available for " +
                             (*iter)->GetFullAtomName() + " (type " +
                             triposType.Type2Str(aType) + ")";
      if (iTrace > 1) {
        std::cout << strError << std::endl;
      }
      if (bHasUndefined) {
        aType = RbtTriposAtomType::UNDEFINED;
      } else {
        throw(RbtFileError(_WHERE_, strError));
      }
    }

    m_ligAtomTypes.push_back(aType);
    if (iTrace > 1) {
      std::cout << "Using grid #" << aType << " for "
                << (*iter)->GetFullAtomName() << std::endl;
    }
  }
}

double RbtVdwGridSF::RawScore() const {
  double score = 0.0;

  // Check grids are defined
  if (m_grids.empty())
    return score;

  // Loop over all ligand atoms
  RbtAtomRListConstIter aIter = m_ligAtomList.begin();
  RbtTriposAtomTypeListConstIter tIter = m_ligAtomTypes.begin();
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
// RbtVdwGridSF
void RbtVdwGridSF::ReadGrids(std::istream &istr) throw(RbtError) {
  m_grids.clear();
  int iTrace = GetTrace();

  // Read header string
  int length;
  Rbt::ReadWithThrow(istr, (char *)&length, sizeof(length));
  char *header = new char[length + 1];
  Rbt::ReadWithThrow(istr, header, length);
  // Add null character to end of string
  header[length] = '\0';
  // Compare title with
  bool match = (_CT == header);
  delete[] header;
  if (!match) {
    throw RbtFileParseError(_WHERE_,
                            "Invalid title string in " + _CT + "::ReadGrids()");
  }

  // Now read number of grids
  int nGrids;
  Rbt::ReadWithThrow(istr, (char *)&nGrids, sizeof(nGrids));
  if (iTrace > 0) {
    std::cout << _CT << ": reading " << nGrids << " grids..." << std::endl;
  }

  // We actually create a vector of size RbtTriposAtomType::MAXTYPES
  // Each grid is in the file is prefixed by the atom type string (e.g. C.2)
  // This string is converted to the corresponding RbtTriposAtomType::eType
  // and the grid stored at m_grids[eType]
  // This is slightly more transferable as it means grids can be read correctly
  // even if the atom type enums change (as long as the atom type string stay
  // the same) It also means we do not have to have a grid for each and every
  // atom type if we don't want to
  m_grids = RbtRealGridList(RbtTriposAtomType::MAXTYPES);
  for (int i = 0; i < nGrids; i++) {
    // Read the atom type string
    Rbt::ReadWithThrow(istr, (char *)&length, sizeof(length));
    char *szType = new char[length + 1];
    Rbt::ReadWithThrow(istr, szType, length);
    // Add null character to end of string
    szType[length] = '\0';
    std::string strType(szType);
    delete[] szType;
    RbtTriposAtomType triposType;
    RbtTriposAtomType::eType aType = triposType.Str2Type(strType);
    if (iTrace > 0) {
      std::cout << "Grid# " << i << "\t"
                << "atom type=" << strType << " (type #" << aType << ")"
                << std::endl;
    }
    // Now we can read the grid
    RbtRealGridPtr spGrid(new RbtRealGrid(istr));
    m_grids[aType] = spGrid;
  }
}

// DM 25 Oct 2000 - track changes to parameter values in local data members
// ParameterUpdated is invoked by RbtParamHandler::SetParameter
void RbtVdwGridSF::ParameterUpdated(const std::string &strName) {
  // DM 25 Oct 2000 - heavily used params
  if (strName == _SMOOTHED) {
    m_bSmoothed = GetParameter(_SMOOTHED);
  } else {
    RbtBaseSF::ParameterUpdated(strName);
  }
}
