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

#include "RbtDihedralSF.h"

#include <functional>

RbtDihedral::RbtDihedral(RbtAtom *pAtom1, RbtAtom *pAtom2, RbtAtom *pAtom3,
                         RbtAtom *pAtom4, const prms &dihprms)
    : m_pAtom1(pAtom1), m_pAtom2(pAtom2), m_pAtom3(pAtom3), m_pAtom4(pAtom4) {
  m_prms.push_back(dihprms);
}

void RbtDihedral::AddTerm(const prms &dihprms) { m_prms.push_back(dihprms); }

double RbtDihedral::operator()() const {
  // std::cout.precision(3);
  // std::cout.setf(std::ios_base::fixed,ios_base::floatfield);
  // std::cout.setf(std::ios_base::right,ios_base::adjustfield);
  double dih = Rbt::BondDihedral(m_pAtom1, m_pAtom2, m_pAtom3, m_pAtom4);
  double score(0.0);
  for (unsigned int i = 0; i != m_prms.size(); ++i) {
    // Subtract the implicit hydrogen offset from the actual dihedral angle
    double dih1 = dih - m_prms[i].offset;
    score +=
        m_prms[i].k *
        (1.0 + m_prms[i].sign * std::cos(m_prms[i].s * dih1 * M_PI / 180.0));
    // std::cout << m_pAtom1->GetName() << "," << m_pAtom2->GetName() <<
    // ","
    // << m_pAtom3->GetName() << ","
    //	 << m_pAtom4->GetName() << "\t" << m_prms[i].sign << "\t" <<
    // m_prms[i].s << "\t" << m_prms[i].k << "\t"
    //	 << dih << "\t" << dih1 << "\t" << score << std::endl;
  }
  return score;
}

// Static data members
std::string RbtDihedralSF::_CT("RbtDihedralSF");
std::string RbtDihedralSF::_IMPL_H_CORR("IMPL_H_CORR");

RbtDihedralSF::RbtDihedralSF() {
#ifdef _DEBUG
  std::cout << _CT << " default constructor" << std::endl;
#endif //_DEBUG
  // Add parameters
  AddParameter(_IMPL_H_CORR, false);
  m_spDihedralSource = RbtParameterFileSourcePtr(new RbtParameterFileSource(
      Rbt::GetRbtFileName("data/sf", "Tripos52_dihedrals.prm")));
  m_centralPairs = m_spDihedralSource->GetSectionList();
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtDihedralSF::~RbtDihedralSF() {
#ifdef _DEBUG
  std::cout << _CT << " destructor" << std::endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

RbtDihedralList RbtDihedralSF::CreateDihedralList(const RbtBondList &bondList) {
  RbtDihedralList dihList;
  int iTrace = GetTrace();
  bool bImplHCorr = GetParameter(_IMPL_H_CORR);

  // Loop over each rotable bond
  for (RbtBondListConstIter iter = bondList.begin(); iter != bondList.end();
       iter++) {
    RbtAtom *pAtom2 = (*iter)->GetAtom1Ptr();
    RbtAtom *pAtom3 = (*iter)->GetAtom2Ptr();
    RbtTriposAtomType::eType t2 = m_triposType(pAtom2);
    RbtTriposAtomType::eType t3 = m_triposType(pAtom3);
    RbtAtomList bondedAtoms2, bondedAtoms3;
    RbtDoubleList offset2, offset3;
    // Get the bonded atoms and implicit hydrogen offsets required to define all
    // the dihedrals across this bond
    CalcBondedAtoms(pAtom2, pAtom3, bondedAtoms2, offset2);
    CalcBondedAtoms(pAtom3, pAtom2, bondedAtoms3, offset3);
    for (RbtAtomListConstIter iter1 = bondedAtoms2.begin();
         iter1 != bondedAtoms2.end(); iter1++) {
      RbtAtom *pAtom1 = *iter1;
      RbtTriposAtomType::eType t1 = m_triposType(pAtom1);
      for (RbtAtomListConstIter iter4 = bondedAtoms3.begin();
           iter4 != bondedAtoms3.end(); iter4++) {
        RbtAtom *pAtom4 = *iter4;
        RbtTriposAtomType::eType t4 = m_triposType(pAtom4);
        RbtDihedral::prms dihprms = FindDihedralParams(t1, t2, t3, t4);
        // Create the base dihedral for this specifier
        RbtDihedral *pDih =
            new RbtDihedral(pAtom1, pAtom2, pAtom3, pAtom4, dihprms);
        if (bImplHCorr) {
          // Add all the ghost dihedral combinations for implicit Hs on pAtom2
          // but only for the first time through the outer loop
          if (iter1 == bondedAtoms2.begin()) {
            for (unsigned int i1 = 0; i1 != offset2.size(); i1++) {
              dihprms = FindDihedralParams(RbtTriposAtomType::H, t2, t3, t4);
              dihprms.offset = offset2[i1];
              if (iTrace > 1) {
                std::cout << _CT << ": <H" << i1 + 1 << ">-"
                          << pAtom2->GetName() << "-" << pAtom3->GetName()
                          << "-" << pAtom4->GetName()
                          << "\toffset = " << dihprms.offset << std::endl;
              }
              pDih->AddTerm(dihprms);
              // Add all the ghost dihedral combinations for implicit Hs on
              // pAtom2 AND on pAtom3 but only for the first time through the
              // outer loop AND inner loop
              if (iter4 == bondedAtoms3.begin()) {
                for (unsigned int i4 = 0; i4 != offset3.size(); i4++) {
                  dihprms = FindDihedralParams(RbtTriposAtomType::H, t2, t3,
                                               RbtTriposAtomType::H);
                  // Combined offset should be offset2 - offset3
                  // but offset3 was determined looking the other way along the
                  // bond so sign is inverted
                  dihprms.offset = offset2[i1] + offset3[i4];
                  if (dihprms.offset >= 180.0)
                    dihprms.offset -= 360.0;
                  if (dihprms.offset < -180.0)
                    dihprms.offset += 360.0;
                  if (iTrace > 1) {
                    std::cout << _CT << ": <H" << i1 + 1 << ">-"
                              << pAtom2->GetName() << "-" << pAtom3->GetName()
                              << "-<H" << i4 + 1
                              << ">\toffset = " << dihprms.offset << std::endl;
                  }
                  pDih->AddTerm(dihprms);
                }
              }
            }
          }
          // Add all the ghost dihedral combinations for implicit Hs on pAtom3
          // but only for the first time through the inner loop
          if (iter4 == bondedAtoms3.begin()) {
            for (unsigned int i4 = 0; i4 != offset3.size(); i4++) {
              dihprms = FindDihedralParams(t1, t2, t3, RbtTriposAtomType::H);
              // Offset should be -offset3
              // Remember offset3 was determined looking the other way along the
              // bond so sign is inverted
              dihprms.offset = offset3[i4];
              if (iTrace > 1) {
                std::cout << _CT << ": " << pAtom1->GetName() << "-"
                          << pAtom2->GetName() << "-" << pAtom3->GetName()
                          << "-" << pAtom4->GetName()
                          << "\toffset = " << dihprms.offset << std::endl;
              }
              pDih->AddTerm(dihprms);
            }
          }
        }
        dihList.push_back(pDih);
      }
    }
  }
  return dihList;
}

RbtDihedral::prms RbtDihedralSF::FindDihedralParams(
    RbtTriposAtomType::eType t1, RbtTriposAtomType::eType t2,
    RbtTriposAtomType::eType t3, RbtTriposAtomType::eType t4) {

  std::string str1 = m_triposType.Type2Str(t1);
  std::string str2 = m_triposType.Type2Str(t2);
  std::string str3 = m_triposType.Type2Str(t3);
  std::string str4 = m_triposType.Type2Str(t4);
  std::string wild = std::string("*");
  int iTrace = GetTrace();
  if (iTrace > 2) {
    std::cout << std::endl
              << _CT << ": Searching for " << str1 << "," << str2 << "," << str3
              << "," << str4 << std::endl;
  }

  std::string strParams;

  std::string str23 = str2 + "_" + str3;
  std::string str32 = str3 + "_" + str2;
  if (std::find(m_centralPairs.begin(), m_centralPairs.end(), str23) !=
      m_centralPairs.end()) {
    m_spDihedralSource->SetSection(str23);
    std::string str14 = str1 + "_" + str4;
    std::string str1w = str1 + "_" + wild;
    std::string strw4 = wild + "_" + str4;
    std::string strww = wild + "_" + wild;
    if (iTrace > 2) {
      std::cout << "Matched " << str23 << " central pair" << std::endl;
    }
    if (m_spDihedralSource->isParameterPresent(str14)) {
      strParams = m_spDihedralSource->GetParameterValueAsString(str14);
      if (iTrace > 2) {
        std::cout << "Matched " << str14 << " outer pair" << std::endl;
      }
    } else if (m_spDihedralSource->isParameterPresent(str1w)) {
      strParams = m_spDihedralSource->GetParameterValueAsString(str1w);
      if (iTrace > 2) {
        std::cout << "Matched " << str1w << " outer pair" << std::endl;
      }
    } else if (m_spDihedralSource->isParameterPresent(strw4)) {
      strParams = m_spDihedralSource->GetParameterValueAsString(strw4);
      if (iTrace > 2) {
        std::cout << "Matched " << strw4 << " outer pair" << std::endl;
      }
    } else if (m_spDihedralSource->isParameterPresent(strww)) {
      strParams = m_spDihedralSource->GetParameterValueAsString(strww);
      if (iTrace > 2) {
        std::cout << "Matched " << strww << " outer pair" << std::endl;
      }
    }
  } else if (std::find(m_centralPairs.begin(), m_centralPairs.end(), str32) !=
             m_centralPairs.end()) {
    m_spDihedralSource->SetSection(str32);
    std::string str41 = str4 + "_" + str1;
    std::string str4w = str4 + "_" + wild;
    std::string strw1 = wild + "_" + str1;
    std::string strww = wild + "_" + wild;
    if (iTrace > 2) {
      std::cout << "Matched " << str32 << " central pair" << std::endl;
    }
    if (m_spDihedralSource->isParameterPresent(str41)) {
      strParams = m_spDihedralSource->GetParameterValueAsString(str41);
      if (iTrace > 2) {
        std::cout << "Matched " << str41 << " outer pair" << std::endl;
      }
    } else if (m_spDihedralSource->isParameterPresent(str4w)) {
      strParams = m_spDihedralSource->GetParameterValueAsString(str4w);
      if (iTrace > 2) {
        std::cout << "Matched " << str4w << " outer pair" << std::endl;
      }
    } else if (m_spDihedralSource->isParameterPresent(strw1)) {
      strParams = m_spDihedralSource->GetParameterValueAsString(strw1);
      if (iTrace > 2) {
        std::cout << "Matched " << strw1 << " outer pair" << std::endl;
      }
    } else if (m_spDihedralSource->isParameterPresent(strww)) {
      strParams = m_spDihedralSource->GetParameterValueAsString(strww);
      if (iTrace > 2) {
        std::cout << "Matched " << strww << " outer pair" << std::endl;
      }
    }
  } else {
    m_spDihedralSource->SetSection();
    strParams = m_spDihedralSource->GetParameterValueAsString("DEFAULT");
    if (iTrace > 0) {
      std::cout << _CT << ": No match for " << str23 << ": using DEFAULT params"
                << std::endl;
    }
  }

  RbtStringList paramList = Rbt::ConvertDelimitedStringToList(strParams);
  // Add checks on #params
  double k = std::atof(paramList[0].c_str());
  double s = std::atof(paramList[1].c_str());
  if (iTrace > 1) {
    std::cout << _CT << ": Assigned " << str1 << "," << str2 << "," << str3
              << "," << str4 << "\tk=" << k << ", s=" << s << std::endl;
  }
  return RbtDihedral::prms(s, k);
}

// Determines the list of bonded atoms for pAtom1, and the list of dihedral
// offsets for any implicit hydrogens pAtom2 is the other central atom in the
// bond, is excluded from the returned list of bonded atoms
void RbtDihedralSF::CalcBondedAtoms(RbtAtom *pAtom1, RbtAtom *pAtom2,
                                    RbtAtomList &bondedAtoms,
                                    RbtDoubleList &offsets) {
  offsets.clear();
  bondedAtoms =
      Rbt::GetAtomList(Rbt::GetBondedAtomList(pAtom1),
                       std::not1(std::bind2nd(Rbt::isAtomPtr_eq(), pAtom2)));
  if (bondedAtoms.empty())
    return;
  int nH = pAtom1->GetNumImplicitHydrogens();
  if (nH == 0)
    return;
  if (GetTrace() > 2) {
    std::cout << _CT << ": determining impl H offsets for " << pAtom1->GetName()
              << " - " << pAtom2->GetName() << std::endl;
  }
  RbtAtom::eHybridState hybrid = pAtom1->GetHybridState();
  if (hybrid == RbtAtom::SP3) {
    // SP3: hardest case first: 2 heavy atoms plus 1 implicit H
    // offset will be approx +/- 120 deg, determine the sign from the improper
    // dihedral already present
    if ((nH == 1) && (bondedAtoms.size() == 2)) {
      double improper =
          Rbt::BondDihedral(bondedAtoms[0], pAtom1, pAtom2, bondedAtoms[1]);
      double offset = -improper;
      if (GetTrace() > 2) {
        std::cout << _CT << ": offset for SP3 atom " << pAtom1->GetName()
                  << " = " << offset << std::endl;
      }
      offsets.push_back(offset);
    }
    // SP3: 1 heavy atom plus 2 implicit H, return both offsets of +/- 120 deg
    else if ((nH == 2) && (bondedAtoms.size() == 1)) {
      if (GetTrace() > 2) {
        std::cout << _CT << ": offsets for SP3 atom " << pAtom1->GetName()
                  << " = -120,+120" << std::endl;
      }
      offsets.push_back(-120.0);
      offsets.push_back(+120.0);
    } else {
      std::cout << _CT << ": WARNING - inconsistent bonding around SP3 atom "
                << pAtom1->GetFullAtomName() << std::endl;
    }
  } else if (hybrid == RbtAtom::SP2) {
    // SP2: The only combination that makes sense is 1 heavy atom and 1 implicit
    // hydrogen Offset is always 180 deg.
    if ((nH == 1) && (bondedAtoms.size() == 1)) {
      if (GetTrace() > 2) {
        std::cout << _CT << ": offset for SP2 atom " << pAtom1->GetName()
                  << " = 180" << std::endl;
      }
      offsets.push_back(-180.0);
    } else {
      std::cout << _CT << ": WARNING - inconsistent bonding around SP2 atom "
                << pAtom1->GetFullAtomName() << std::endl;
    }
  } else {
    std::cout << _CT << ": WARNING - unexpected hybridisation state for "
              << pAtom1->GetFullAtomName() << std::endl;
  }
}
