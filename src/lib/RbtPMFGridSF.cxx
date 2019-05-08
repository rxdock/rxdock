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

#include "RbtPMFGridSF.h"
#include "RbtFileError.h"
#include "RbtWorkSpace.h"

std::string RbtPMFGridSF::_CT("RbtPMFGridSF");
std::string RbtPMFGridSF::_GRID("GRID");
std::string RbtPMFGridSF::_SMOOTHED("SMOOTHED");

RbtPMFGridSF::RbtPMFGridSF(const std::string &strName)
    : RbtBaseSF(_CT, strName), m_bSmoothed(true) {
  AddParameter(_GRID, ".grd");
  AddParameter(_SMOOTHED, m_bSmoothed);

  cout << _CT << " parameterised constructor" << endl;

  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtPMFGridSF::~RbtPMFGridSF() {
  cout << _CT << " destructor" << endl;
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

void RbtPMFGridSF::SetupReceptor() {
  cout << _CT << "::SetupReceptor()" << endl;
  theGrids.clear();

  if (GetReceptor().Null())
    return;

  std::string strWSName = GetWorkSpace()->GetName();

  std::string strSuffix = GetParameter(_GRID);
  std::string strFile =
      Rbt::GetRbtFileName("data/grids", strWSName + strSuffix);
  cout << _CT << " Reading PMF grid from " << strFile << endl;
  ifstream istr(strFile.c_str(), ios_base::in | ios_base::binary);
  ReadGrids(istr);
  istr.close();
}
// Determine PMF grid type for each atom
void RbtPMFGridSF::SetupLigand() {
  theLigandList.clear();
  if (GetLigand().Null())
    return;

  // get the  non-H atoms only
  theLigandList = Rbt::GetAtomList(GetLigand()->GetAtomList(),
                                   std::not1(Rbt::isAtomicNo_eq(1)));
#ifdef _DEBUG
  cout << _CT << "::SetupLigand(): #ATOMS = " << theLigandList.size() << endl;
#endif //_DEBUG
}

double RbtPMFGridSF::RawScore() const {
  double theScore = 0.0;
  if (theGrids.empty()) // if no grids defined
    return theScore;

  // Loop over all ligand atoms
  RbtAtomListConstIter iter = theLigandList.begin();
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

void RbtPMFGridSF::ReadGrids(istream &istr) throw(RbtError) {
  cout << "**************************************************************"
       << endl;
  cout << "                    Reading grid..." << endl;
  cout << "**************************************************************"
       << endl;
  cout << _CT << "::ReadGrids(istream&)" << endl;
  theGrids.clear();

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
  cout << "Reading " << nGrids << " grids..." << endl;
  theGrids.reserve(nGrids);
  for (int i = CF; i <= nGrids; i++) {
    cout << "Grid# " << i << " ";
    // Read type
    RbtPMFType theType;
    Rbt::ReadWithThrow(istr, (char *)&theType, sizeof(theType));
    cout << "type " << Rbt::PMFType2Str(theType);
    // Now we can read the grid
    RbtRealGridPtr spGrid(new RbtRealGrid(istr));
    cout << " done" << endl;
    theGrids.push_back(spGrid);
  }
}

// since there is no  HH,HL,Fe,V,Mn grid, we have to correct
unsigned int RbtPMFGridSF::GetCorrectedType(RbtPMFType aType) const {
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
