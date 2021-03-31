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

#include "rxdock/DihedralSF.h"

#include <loguru.hpp>

#include <functional>

using namespace rxdock;

DihedralElement::DihedralElement(Atom *pAtom1, Atom *pAtom2, Atom *pAtom3,
                                 Atom *pAtom4, const prms &dihprms)
    : m_pAtom1(pAtom1), m_pAtom2(pAtom2), m_pAtom3(pAtom3), m_pAtom4(pAtom4) {
  m_prms.push_back(dihprms);
}

void DihedralElement::AddTerm(const prms &dihprms) {
  m_prms.push_back(dihprms);
}

double DihedralElement::operator()() const {
  double dih = BondDihedral(m_pAtom1, m_pAtom2, m_pAtom3, m_pAtom4);
  double score(0.0);
  for (unsigned int i = 0; i != m_prms.size(); ++i) {
    // Subtract the implicit hydrogen offset from the actual dihedral angle
    double dih1 = dih - m_prms[i].offset;
    score +=
        m_prms[i].k *
        (1.0 + m_prms[i].sign * std::cos(m_prms[i].s * dih1 * M_PI / 180.0));
    LOG_F(1, "{},{},{},{}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}",
          m_pAtom1->GetName(), m_pAtom2->GetName(), m_pAtom3->GetName(),
          m_pAtom4->GetName(), m_prms[i].sign, m_prms[i].s, m_prms[i].k, dih,
          dih1, score);
  }
  return score;
}

// Static data members
const std::string DihedralSF::_CT = "DihedralSF";
const std::string DihedralSF::_IMPL_H_CORR = "IMPL_H_CORR";

DihedralSF::DihedralSF() {
  LOG_F(2, "DihedralSF default constructor");
  // Add parameters
  AddParameter(_IMPL_H_CORR, false);
  m_spDihedralSource = ParameterFileSourcePtr(new ParameterFileSource(
      GetDataFileName("data/sf", "Tripos52_dihedrals.prm")));
  m_centralPairs = m_spDihedralSource->GetSectionList();
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

DihedralSF::~DihedralSF() {
  LOG_F(2, "DihedralSF destructor");
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

DihedralList DihedralSF::CreateDihedralList(const BondList &bondList) {
  DihedralList dihList;
  bool bImplHCorr = GetParameter(_IMPL_H_CORR);

  // Loop over each rotable bond
  for (BondListConstIter iter = bondList.begin(); iter != bondList.end();
       iter++) {
    Atom *pAtom2 = (*iter)->GetAtom1Ptr();
    Atom *pAtom3 = (*iter)->GetAtom2Ptr();
    TriposAtomType::eType t2 = m_triposType(pAtom2);
    TriposAtomType::eType t3 = m_triposType(pAtom3);
    AtomList bondedAtoms2, bondedAtoms3;
    std::vector<double> offset2, offset3;
    // Get the bonded atoms and implicit hydrogen offsets required to define all
    // the dihedrals across this bond
    CalcBondedAtoms(pAtom2, pAtom3, bondedAtoms2, offset2);
    CalcBondedAtoms(pAtom3, pAtom2, bondedAtoms3, offset3);
    for (AtomListConstIter iter1 = bondedAtoms2.begin();
         iter1 != bondedAtoms2.end(); iter1++) {
      Atom *pAtom1 = *iter1;
      TriposAtomType::eType t1 = m_triposType(pAtom1);
      for (AtomListConstIter iter4 = bondedAtoms3.begin();
           iter4 != bondedAtoms3.end(); iter4++) {
        Atom *pAtom4 = *iter4;
        TriposAtomType::eType t4 = m_triposType(pAtom4);
        DihedralElement::prms dihprms = FindDihedralParams(t1, t2, t3, t4);
        // Create the base dihedral for this specifier
        DihedralElement *pDih =
            new DihedralElement(pAtom1, pAtom2, pAtom3, pAtom4, dihprms);
        if (bImplHCorr) {
          // Add all the ghost dihedral combinations for implicit Hs on pAtom2
          // but only for the first time through the outer loop
          if (iter1 == bondedAtoms2.begin()) {
            for (unsigned int i1 = 0; i1 != offset2.size(); i1++) {
              dihprms = FindDihedralParams(TriposAtomType::H, t2, t3, t4);
              dihprms.offset = offset2[i1];
              LOG_F(
                  1,
                  "DihedralSF::CreateDihedralList: <H{}>-{}-{}-{}\toffset = {}",
                  i1 + 1, pAtom2->GetName(), pAtom3->GetName(),
                  pAtom4->GetName(), dihprms.offset);
              pDih->AddTerm(dihprms);
              // Add all the ghost dihedral combinations for implicit Hs on
              // pAtom2 AND on pAtom3 but only for the first time through the
              // outer loop AND inner loop
              if (iter4 == bondedAtoms3.begin()) {
                for (unsigned int i4 = 0; i4 != offset3.size(); i4++) {
                  dihprms = FindDihedralParams(TriposAtomType::H, t2, t3,
                                               TriposAtomType::H);
                  // Combined offset should be offset2 - offset3
                  // but offset3 was determined looking the other way along the
                  // bond so sign is inverted
                  dihprms.offset = offset2[i1] + offset3[i4];
                  if (dihprms.offset >= 180.0)
                    dihprms.offset -= 360.0;
                  if (dihprms.offset < -180.0)
                    dihprms.offset += 360.0;
                  LOG_F(1,
                        "DihedralSF::CreateDihedralList: "
                        "<H{}>-{}-{}-<H{}>\toffset = {}",
                        i1 + 1, pAtom2->GetName(), pAtom3->GetName(), i4 + 1,
                        dihprms.offset);
                  pDih->AddTerm(dihprms);
                }
              }
            }
          }
          // Add all the ghost dihedral combinations for implicit Hs on pAtom3
          // but only for the first time through the inner loop
          if (iter4 == bondedAtoms3.begin()) {
            for (unsigned int i4 = 0; i4 != offset3.size(); i4++) {
              dihprms = FindDihedralParams(t1, t2, t3, TriposAtomType::H);
              // Offset should be -offset3
              // Remember offset3 was determined looking the other way along the
              // bond so sign is inverted
              dihprms.offset = offset3[i4];
              LOG_F(1,
                    "DihedralSF::CreateDihedralList: {}-{}-{}-{}\toffset = {}",
                    pAtom1->GetName(), pAtom2->GetName(), pAtom3->GetName(),
                    pAtom4->GetName(), dihprms.offset);
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

DihedralElement::prms DihedralSF::FindDihedralParams(TriposAtomType::eType t1,
                                                     TriposAtomType::eType t2,
                                                     TriposAtomType::eType t3,
                                                     TriposAtomType::eType t4) {

  std::string str1 = m_triposType.Type2Str(t1);
  std::string str2 = m_triposType.Type2Str(t2);
  std::string str3 = m_triposType.Type2Str(t3);
  std::string str4 = m_triposType.Type2Str(t4);
  std::string wild = std::string("*");
  LOG_F(1, "DihedralSF::CreateDihedralList: Searching for {},{},{},{}", str1,
        str2, str3, str4);

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
    LOG_F(1, "Matched {} central pair", str23);
    if (m_spDihedralSource->isParameterPresent(str14)) {
      strParams = m_spDihedralSource->GetParameterValueAsString(str14);
      LOG_F(1, "Matched {} outer pair", str14);
    } else if (m_spDihedralSource->isParameterPresent(str1w)) {
      strParams = m_spDihedralSource->GetParameterValueAsString(str1w);
      LOG_F(1, "Matched {} outer pair", str1w);
    } else if (m_spDihedralSource->isParameterPresent(strw4)) {
      strParams = m_spDihedralSource->GetParameterValueAsString(strw4);
      LOG_F(1, "Matched {} outer pair", strw4);
    } else if (m_spDihedralSource->isParameterPresent(strww)) {
      strParams = m_spDihedralSource->GetParameterValueAsString(strww);
      LOG_F(1, "Matched {} outer pair", strww);
    }
  } else if (std::find(m_centralPairs.begin(), m_centralPairs.end(), str32) !=
             m_centralPairs.end()) {
    m_spDihedralSource->SetSection(str32);
    std::string str41 = str4 + "_" + str1;
    std::string str4w = str4 + "_" + wild;
    std::string strw1 = wild + "_" + str1;
    std::string strww = wild + "_" + wild;
    LOG_F(1, "Matched {} central pair", str32);
    if (m_spDihedralSource->isParameterPresent(str41)) {
      strParams = m_spDihedralSource->GetParameterValueAsString(str41);
      LOG_F(1, "Matched {} outer pair", str41);
    } else if (m_spDihedralSource->isParameterPresent(str4w)) {
      strParams = m_spDihedralSource->GetParameterValueAsString(str4w);
      LOG_F(1, "Matched {} outer pair", str4w);
    } else if (m_spDihedralSource->isParameterPresent(strw1)) {
      strParams = m_spDihedralSource->GetParameterValueAsString(strw1);
      LOG_F(1, "Matched {} outer pair", strw1);
    } else if (m_spDihedralSource->isParameterPresent(strww)) {
      strParams = m_spDihedralSource->GetParameterValueAsString(strww);
      LOG_F(1, "Matched {} outer pair", strww);
    }
  } else {
    m_spDihedralSource->SetSection();
    strParams = m_spDihedralSource->GetParameterValueAsString("DEFAULT");
    LOG_F(
        WARNING,
        "DihedralSF::FindDihedralParams: No match for {}: using DEFAULT params",
        str23);
  }

  std::vector<std::string> paramList = ConvertDelimitedStringToList(strParams);
  // Add checks on #params
  double k = std::atof(paramList[0].c_str());
  double s = std::atof(paramList[1].c_str());
  LOG_F(1, "DihedralSF::FindDihedralParams: Assigned {},{},{},{}\tk={}, s={}",
        str1, str2, str3, str4, k, s);
  return DihedralElement::prms(s, k);
}

// Determines the list of bonded atoms for pAtom1, and the list of dihedral
// offsets for any implicit hydrogens pAtom2 is the other central atom in the
// bond, is excluded from the returned list of bonded atoms
void DihedralSF::CalcBondedAtoms(Atom *pAtom1, Atom *pAtom2,
                                 AtomList &bondedAtoms,
                                 std::vector<double> &offsets) {
  offsets.clear();
  bondedAtoms = GetAtomListWithPredicate(
      GetBondedAtomList(pAtom1),
      std::bind(std::not2(isAtomPtr_eq()), std::placeholders::_1, pAtom2));
  if (bondedAtoms.empty())
    return;
  int nH = pAtom1->GetNumImplicitHydrogens();
  if (nH == 0)
    return;
  LOG_F(1,
        "DihedralSF::CalcBondedAtoms: determining impl H offsets for {} - {}",
        pAtom1->GetName(), pAtom2->GetName());
  Atom::eHybridState hybrid = pAtom1->GetHybridState();
  if (hybrid == Atom::SP3) {
    // SP3: hardest case first: 2 heavy atoms plus 1 implicit H
    // offset will be approx +/- 120 deg, determine the sign from the improper
    // dihedral already present
    if ((nH == 1) && (bondedAtoms.size() == 2)) {
      double improper =
          BondDihedral(bondedAtoms[0], pAtom1, pAtom2, bondedAtoms[1]);
      double offset = -improper;
      LOG_F(1, "DihedralSF::CalcBondedAtoms: offset for SP3 atom {} = {}",
            pAtom1->GetName(), offset);
      offsets.push_back(offset);
    }
    // SP3: 1 heavy atom plus 2 implicit H, return both offsets of +/- 120 deg
    else if ((nH == 2) && (bondedAtoms.size() == 1)) {
      LOG_F(1,
            "DihedralSF::CalcBondedAtoms: offsets for SP3 atom {} = -120,+120",
            pAtom1->GetName());
      offsets.push_back(-120.0);
      offsets.push_back(+120.0);
    } else {
      LOG_F(WARNING,
            "DihedralSF::CalcBondedAtoms: inconsistent bonding around SP3 atom "
            "{}",
            pAtom1->GetFullAtomName());
    }
  } else if (hybrid == Atom::SP2) {
    // SP2: The only combination that makes sense is 1 heavy atom and 1 implicit
    // hydrogen Offset is always 180 deg.
    if ((nH == 1) && (bondedAtoms.size() == 1)) {
      LOG_F(1, "DihedralSF::CalcBondedAtoms: offset for SP2 atom {} = 180",
            pAtom1->GetName());
      offsets.push_back(-180.0);
    } else {
      LOG_F(WARNING,
            "DihedralSF::CalcBondedAtoms: inconsistent bonding around SP2 atom "
            "{}",
            pAtom1->GetFullAtomName());
    }
  } else {
    LOG_F(WARNING,
          "DihedralSF::CalcBondedAtoms: unexpected hybridisation state for {}",
          pAtom1->GetFullAtomName());
  }
}
