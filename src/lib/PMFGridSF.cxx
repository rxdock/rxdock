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

#include "PMFGridSF.h"
#include "FileError.h"
#include "WorkSpace.h"

#include <functional>

using namespace rxdock;

std::string PMFGridSF::_CT("PMFGridSF");
std::string PMFGridSF::_GRID("GRID");
std::string PMFGridSF::_SMOOTHED("SMOOTHED");

PMFGridSF::PMFGridSF(const std::string &strName)
    : BaseSF(_CT, strName), m_bSmoothed(true) {
  AddParameter(_GRID, ".grd");
  AddParameter(_SMOOTHED, m_bSmoothed);

  std::cout << _CT << " parameterised constructor" << std::endl;

  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

PMFGridSF::~PMFGridSF() {
  std::cout << _CT << " destructor" << std::endl;
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

void PMFGridSF::SetupReceptor() {
  std::cout << _CT << "::SetupReceptor()" << std::endl;
  theGrids.clear();

  if (GetReceptor().Null())
    return;

  std::string strWSName = GetWorkSpace()->GetName();

  std::string strSuffix = GetParameter(_GRID);
  std::string strFile = GetDataFileName("data/grids", strWSName + strSuffix);
  std::cout << _CT << " Reading PMF grid from " << strFile << std::endl;
  std::ifstream istr(strFile.c_str(),
                     std::ios_base::in | std::ios_base::binary);
  ReadGrids(istr);
  istr.close();
}
// Determine PMF grid type for each atom
void PMFGridSF::SetupLigand() {
  theLigandList.clear();
  if (GetLigand().Null())
    return;

  // get the  non-H atoms only
  theLigandList = GetAtomListWithPredicate(GetLigand()->GetAtomList(),
                                           std::not1(isAtomicNo_eq(1)));
#ifdef _DEBUG
  std::cout << _CT << "::SetupLigand(): #ATOMS = " << theLigandList.size()
            << std::endl;
#endif //_DEBUG
}

double PMFGridSF::RawScore() const {
  double theScore = 0.0;
  if (theGrids.empty()) // if no grids defined
    return theScore;

  // Loop over all ligand atoms
  AtomListConstIter iter = theLigandList.begin();
  if (m_bSmoothed) {
    for (unsigned int i = 0; iter != theLigandList.end(); iter++, i++) {
      unsigned int theType = GetCorrectedType((*iter)->GetPMFType());
      double score =
          theGrids[theType - 1]->GetSmoothedValue((*iter)->GetCoords());
      theScore += score;
    }
  } else {
    for (unsigned int i = 0; iter != theLigandList.end(); iter++, i++) {
      unsigned int theType = GetCorrectedType((*iter)->GetPMFType());
      double score = theGrids[theType - 1]->GetValue((*iter)->GetCoords());
      theScore += score;
    }
  }
  return theScore;
}

void PMFGridSF::ReadGrids(std::istream &istr) {
  std::cout << "**************************************************************"
            << std::endl;
  std::cout << "                    Reading grid..." << std::endl;
  std::cout << "**************************************************************"
            << std::endl;
  std::cout << _CT << "::ReadGrids(std::istream&)" << std::endl;
  theGrids.clear();

  // Read header string
  int length;
  ReadWithThrow(istr, (char *)&length, sizeof(length));
  char *header = new char[length + 1];
  ReadWithThrow(istr, header, length);
  // Add null character to end of string
  header[length] = '\0';
  // Compare title with
  bool match = (_CT == header);
  delete[] header;
  if (!match) {
    throw FileParseError(_WHERE_,
                         "Invalid title string in " + _CT + "::ReadGrids()");
  }

  // Now read number of grids
  int nGrids;
  ReadWithThrow(istr, (char *)&nGrids, sizeof(nGrids));
  std::cout << "Reading " << nGrids << " grids..." << std::endl;
  theGrids.reserve(nGrids);
  for (int i = CF; i <= nGrids; i++) {
    std::cout << "Grid# " << i << " ";
    // Read type
    PMFType theType;
    ReadWithThrow(istr, (char *)&theType, sizeof(theType));
    std::cout << "type " << PMFType2Str(theType);
    // Now we can read the grid
    RealGridPtr spGrid(new RealGrid(istr));
    std::cout << " done" << std::endl;
    theGrids.push_back(spGrid);
  }
}

// since there is no  HH,HL,Fe,V,Mn grid, we have to correct
unsigned int PMFGridSF::GetCorrectedType(PMFType aType) const {
  if (aType < HL)
    return (unsigned int)aType;
  else if (aType < Mn)
    return (unsigned int)aType -
           1; // HL and HH but list starts at 0 still PMFType starts at 1
  else if (aType < Fe)
    return (unsigned int)aType - 2; // HL,HH,Mn (3-1)
  else
    return (unsigned int)aType -
           3; // HL,HH,Mn,Fe - V is the very last in the list
}
