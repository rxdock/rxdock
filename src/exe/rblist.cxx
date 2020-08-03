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

// SD file dump, listing number of atoms, bonds, interaction centers etc

#include "MdlFileSink.h"
#include "MdlFileSource.h"
#include <algorithm> //for min,max
#include <functional>
#include <iomanip>
#include <iostream>

using namespace rxdock;

namespace rxdock {

class isAtomNitroN : public std::function<bool(AtomPtr)> {
  // Useful predicates
  isFFType_eq bIsN_SP2P;
  isFFType_eq bIsO_SP2;
  isFFType_eq bIsO_SP3M;

  // DM 11 Jul 2001 - checks for nitro nitrogen (NSP2+ bonded to OSP3- and OSP2)
public:
  explicit isAtomNitroN()
      : bIsN_SP2P("N_SP2+"), bIsO_SP2("O_SP2"), bIsO_SP3M("O_SP3-") {}
  bool operator()(AtomPtr spAtom) const {
    if (!bIsN_SP2P(spAtom))
      return false;
    AtomList bondedAtomList = GetBondedAtomList(spAtom);
    return (GetNumAtomsWithPredicate(bondedAtomList, bIsO_SP2) == 1) &&
           (GetNumAtomsWithPredicate(bondedAtomList, bIsO_SP3M) == 1);
  }
};

void ListAtoms(const AtomList &atomList, std::string strTitle,
               int nMaxToPrint = 20) {
  int nSize = atomList.size();
  int nToPrint = std::min(nMaxToPrint, nSize);
  std::cout << std::endl
            << "1st " << nToPrint << " " << strTitle << " (out of " << nSize
            << ")" << std::endl;
  for (AtomListConstIter iter = atomList.begin();
       (iter != (atomList.begin() + nMaxToPrint)) && (iter != atomList.end());
       iter++)
    std::cout << (**iter) << std::endl;
}

void ListBonds(const BondList &bondList, std::string strTitle,
               int nMaxToPrint = 20) {
  int nSize = bondList.size();
  int nToPrint = std::min(nMaxToPrint, nSize);
  std::cout << std::endl
            << "1st " << nToPrint << " " << strTitle << " (out of " << nSize
            << ")" << std::endl;
  for (BondListConstIter iter = bondList.begin();
       (iter != (bondList.begin() + nMaxToPrint)) && (iter != bondList.end());
       iter++)
    std::cout << (**iter) << std::endl;
}

// DM 7 Jun 1999 - set primary and secondary amide bonds to 180.0 deg
void CheckAmideBonds(ModelPtr spModel) {
  // Useful predicates
  isFFType_eq bIsN_TRI("N_TRI");
  isFFType_eq bIsC_SP2("C_SP2");
  isFFType_eq bIsO_SP2("O_SP2");
  isAtomicNo_eq bIsH(1);

  BondList amideBondList =
      GetBondListWithPredicate(spModel->GetBondList(), isBondAmide());
  // DM 9 Oct 2000 - remove cyclic bonds from list
  amideBondList =
      GetBondListWithPredicate(amideBondList, std::not1(isBondCyclic()));
  int nAMIDE = amideBondList.size();
  std::cout << nAMIDE << " amide bonds:" << std::endl;
  for (BondListConstIter bIter = amideBondList.begin();
       bIter != amideBondList.end(); bIter++) {
    AtomPtr spAtom2 = (*bIter)->GetAtom1Ptr();
    AtomPtr spAtom3 = (*bIter)->GetAtom2Ptr();

    // Find N and C
    AtomPtr spN = (bIsN_TRI(spAtom2)) ? spAtom2 : spAtom3;
    AtomPtr spC = (bIsC_SP2(spAtom2)) ? spAtom2 : spAtom3;

    // Find bonded O_SP2 and H (if any)
    AtomList bondedAtomsN = GetBondedAtomList(spN);
    AtomList bondedAtomsC = GetBondedAtomList(spC);
    AtomListIter oIter = FindAtomInList(bondedAtomsC, bIsO_SP2);
    AtomListIter hIter = FindAtomInList(bondedAtomsN, bIsH);

    if ((oIter != bondedAtomsC.end()) && (hIter != bondedAtomsN.end())) {
      AtomPtr spH = (*hIter);
      AtomPtr spO = (*oIter);
      double phi = BondDihedral(spH, spN, spC, spO);
      std::cout << spH->GetName() << "(" << spH->GetFFType() << ") - "
                << spN->GetName() << "(" << spN->GetFFType() << ") - "
                << spC->GetName() << "(" << spC->GetFFType() << ") - "
                << spO->GetName() << "(" << spO->GetFFType() << "): " << phi
                << " deg" << std::endl;
      double deltaPhi = 180.0 - phi;
      std::cout << "Rotating bond by " << deltaPhi << " deg" << std::endl;
      spModel->RotateBond(*bIter, deltaPhi);
      phi = BondDihedral(spH, spN, spC, spO);
      std::cout << "Dihedral is now " << phi << " deg" << std::endl;
    } else {
      std::cout << spN->GetName() << "(" << spN->GetFFType() << ") - "
                << spC->GetName() << "(" << spC->GetFFType() << "): TERTIARY"
                << std::endl;
    }
  }
}

} // namespace rxdock

int main(int argc, char *argv[]) {
  std::cout.setf(std::ios_base::left, std::ios_base::adjustfield);

  // Strip off the path to the executable, leaving just the file name
  std::string strExeName(argv[0]);
  std::string::size_type i = strExeName.rfind("/");
  if (i != std::string::npos)
    strExeName.erase(0, i + 1);

  // Print a standard header
  PrintStdHeader(std::cout, strExeName);

  // Default values for optional arguments
  std::string strInputSDFile;
  std::string strOutputSDFile;
  bool bPosIonise(false);
  bool bNegIonise(false);
  bool bImplH(true); // if true, read only polar hydrogens from SD file, else
                     // read all H's present
  bool bCheckAmides(false); // if true, sets all 2ndry amide bonds to trans
  bool bList(false);

  // Display brief help message if no args
  if (argc == 1) {
    std::cout
        << std::endl
        << strExeName
        << " - output interaction center info for ligands in SD file (with "
           "optional autoionisation)"
        << std::endl;
    std::cout << std::endl
              << "Usage:\t" << strExeName
              << " -i<InputSDFile> [-o<OutputSDFile>] [-ap] [-an] [-allH]"
              << std::endl;
    std::cout << std::endl
              << "Options:\t-i<InputSDFile> - input ligand SD file"
              << std::endl;
    std::cout << "\t\t-o<OutputSDFile> - output SD file with descriptors "
                 "(default=no output)"
              << std::endl;
    std::cout
        << "\t\t-ap - protonate all neutral amines, guanidines, imidazoles "
           "(default=disabled)"
        << std::endl;
    std::cout
        << "\t\t-an - deprotonate all carboxylic, sulphur and phosphorous "
           "acid groups (default=disabled)"
        << std::endl;
    std::cout
        << "\t\t-allH - read all hydrogens present (default=polar hydrogens "
           "only)"
        << std::endl;
    std::cout
        << "\t\t-tr - rotate all 2ndry amides to trans (default=leave alone)"
        << std::endl;
    std::cout << "\t\t-l - verbose listing of ligand atoms and rotable bonds "
                 "(default = compact table format)"
              << std::endl;
    return 1;
  }

  // Check command line arguments
  std::cout << std::endl << "Command line args:" << std::endl;
  for (int iarg = 1; iarg < argc; iarg++) {
    std::cout << argv[iarg];
    std::string strArg(argv[iarg]);
    if (strArg.find("-i") == 0)
      strInputSDFile = strArg.substr(2);
    else if (strArg.find("-o") == 0)
      strOutputSDFile = strArg.substr(2);
    else if (strArg.find("-ap") == 0)
      bPosIonise = true;
    else if (strArg.find("-an") == 0)
      bNegIonise = true;
    else if (strArg.find("-allH") == 0)
      bImplH = false;
    else if (strArg.find("-tr") == 0)
      bCheckAmides = true;
    else if (strArg.find("-l") == 0)
      bList = true;
    else {
      std::cout << " ** INVALID ARGUMENT" << std::endl;
      return 1;
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;

  bool bWriteLigand = !strOutputSDFile.empty();

  if (!bList) {
    std::cout << std::setw(8) << "RECORD" << std::setw(30) << "NAME"
              << std::setw(8) << "#ATOMS" << std::setw(8) << "#BONDS"
              << std::setw(8) << "MW" << std::setw(8) << "#ROT" << std::setw(8)
              << "#HBD" << std::setw(8) << "#HBA" << std::setw(8) << "#AROM"
              << std::setw(8) << "#GUAN" << std::setw(8) << "#LIPOC"
              << std::endl;
  }

  try {
    // Prepare the SD file sink for saving the autoionised structures for each
    // ligand
    MolecularFileSinkPtr spMdlFileSink;
    if (bWriteLigand)
      spMdlFileSink =
          MolecularFileSinkPtr(new MdlFileSink(strOutputSDFile, ModelPtr()));

    ///////////////////////////////////
    // MAIN LOOP OVER LIGAND RECORDS
    ///////////////////////////////////
    MolecularFileSourcePtr spMdlFileSource(
        new MdlFileSource(strInputSDFile, bPosIonise, bNegIonise, bImplH));
    // DM 4 June 1999 - only read the largest segment
    // spMdlFileSource->SetSegmentFilterMap(ConvertStringToSegmentMap("H"));

    for (int nRec = 1; spMdlFileSource->FileStatusOK();
         spMdlFileSource->NextRecord(), nRec++) {
      Error molStatus = spMdlFileSource->Status();
      if (!molStatus.isOK()) {
        std::cout << molStatus.what() << std::endl;
        continue;
      }

      // DM 26 Jul 1999 - only read the largest segment (guaranteed to be called
      // H)
      spMdlFileSource->SetSegmentFilterMap(ConvertStringToSegmentMap("H"));

      // Create the ligand model
      ModelPtr spModel(new Model(spMdlFileSource));
      std::string strModelName = spModel->GetTitleList()[0];
      if (strModelName.empty())
        strModelName = spModel->GetName();

      if (bList) {
        AtomList atomList = spModel->GetAtomList();
        std::cout << std::endl << "Model = " << spModel->GetName() << std::endl;
        for (AtomListConstIter iter = atomList.begin(); iter != atomList.end();
             iter++) {
          std::cout << (**iter) << std::endl;
        }
        // continue;
      }

      // DM 11 Nov 1999
      double dMolWt = spModel->GetTotalMass();

      // Determine the interaction center counts
      int nAtoms = spModel->GetNumAtoms();

      int nBonds = spModel->GetNumBonds();
      AtomList atomList = spModel->GetAtomList();
      BondList bondList = spModel->GetBondList();

      isHybridState_eq bIsArom(Atom::AROM);
      int nLipoC = GetNumAtomsWithPredicate(atomList, isAtomLipophilic());
      int nAromAtoms =
          GetNumAtomsWithPredicate(atomList, bIsArom); //# aromatic atoms
      int nNHBD = GetNumAtomsWithPredicate(atomList, isAtomHBondDonor());
      int nMetal = GetNumAtomsWithPredicate(atomList, isAtomMetal());
      int nGuan = GetNumAtomsWithPredicate(atomList, isAtomGuanidiniumCarbon());
      int nNHBA = GetNumAtomsWithPredicate(atomList, isAtomHBondAcceptor());

      // Total +ve and -ve charges
      double posChg(0.0);
      double negChg(0.0);

      for (AtomListConstIter iter = atomList.begin(); iter != atomList.end();
           iter++) {
        double chg = (*iter)->GetGroupCharge();
        if (chg > 0.0) {
          posChg += chg;
        } else if (chg < 0.0) {
          negChg += chg;
        }
      }

      // Count aromatic rings
      AtomListList ringLists = spModel->GetRingAtomLists();
      int nAromRings(0); //# aromatic rings
      for (AtomListListIter rIter = ringLists.begin(); rIter != ringLists.end();
           rIter++) {
        if (GetNumAtomsWithPredicate(*rIter, isPiAtom()) == (*rIter).size())
          nAromRings++;
      }

      BondList rotatableBondList =
          GetBondListWithPredicate(bondList, isBondRotatable());
      // DM 24 Sep 2001 - remove bonds to terminal NH3+
      rotatableBondList =
          GetBondListWithPredicate(rotatableBondList, std::not1(isBondToNH3()));
      // DM 19 Sep 2002 - remove bonds to terminal OH
      rotatableBondList =
          GetBondListWithPredicate(rotatableBondList, std::not1(isBondToOH()));
      int nROT = rotatableBondList.size();

      if (bList) {
        std::cout << nROT << " rotatable bonds:" << std::endl;
        for (BondListConstIter bIter = rotatableBondList.begin();
             bIter != rotatableBondList.end(); bIter++) {
          AtomPtr spAtom1 = (*bIter)->GetAtom1Ptr();
          AtomPtr spAtom2 = (*bIter)->GetAtom2Ptr();
          std::cout << spAtom1->GetName() << "(" << spAtom1->GetFFType()
                    << ") - " << spAtom2->GetName() << "("
                    << spAtom2->GetFFType() << ")" << std::endl;
        }
      } else {
        // Output the ligand info
        std::cout << std::setw(8) << nRec << std::setw(30)
                  << strModelName.c_str() << std::setw(8) << nAtoms
                  << std::setw(8) << nBonds << std::setw(8) << dMolWt
                  << std::setw(8) << nROT << std::setw(8) << nNHBD
                  << std::setw(8) << nNHBA << std::setw(8) << nAromRings
                  << std::setw(8) << nGuan << std::setw(8) << nLipoC
                  << std::endl;
      }

      // Dump model to MdlFileSink
      if (bWriteLigand) {
        if (bCheckAmides) {
          CheckAmideBonds(spModel); // DM 7 June 1999 - check and fix amide
                                    // bonds
        }
        // DM 25 Jul 2002 - do not change ligand coords
        // spModel->AlignPrincipalAxes();//Align principal axes with Cartesian
        // axes

        // DM 11 Nov 1999 Set Data fields for mol.wt, #interaction centers etc
        // NB Molecular weight is calculated for the fragment actually read
        // i.e. largest fragment present if multifragment record.
        spModel->SetDataValue("LIG_MW", Variant(dMolWt));
        spModel->SetDataValue("LIG_NATOMS", Variant(nAtoms));
        spModel->SetDataValue("LIG_NLIPOC", Variant(nLipoC));
        spModel->SetDataValue("LIG_NAROMATOMS", Variant(nAromAtoms));
        spModel->SetDataValue("LIG_NAROMRINGS", Variant(nAromRings));
        spModel->SetDataValue("LIG_NHBD", Variant(nNHBD));
        spModel->SetDataValue("LIG_NMETAL", Variant(nMetal));
        spModel->SetDataValue("LIG_NGUAN", Variant(nGuan));
        spModel->SetDataValue("LIG_NHBA", Variant(nNHBA));
        spModel->SetDataValue("LIG_NROT", Variant(nROT));

        // Atom type percentages
        spModel->SetDataValue("LIG_PERC_LIPOC",
                              Variant(100.0 * nLipoC / nAtoms));
        spModel->SetDataValue("LIG_PERC_AROMATOMS",
                              Variant(100.0 * nAromAtoms / nAtoms));
        spModel->SetDataValue("LIG_PERC_HBD", Variant(100.0 * nNHBD / nAtoms));
        spModel->SetDataValue("LIG_PERC_METAL",
                              Variant(100.0 * nMetal / nAtoms));
        spModel->SetDataValue("LIG_PERC_GUAN", Variant(100.0 * nGuan / nAtoms));
        spModel->SetDataValue("LIG_PERC_HBA", Variant(100.0 * nNHBA / nAtoms));

        // Charges
        spModel->SetDataValue("LIG_POS_CHG", Variant(posChg));
        spModel->SetDataValue("LIG_NEG_CHG", Variant(negChg));
        spModel->SetDataValue("LIG_TOT_CHG", Variant(posChg + negChg));

        spMdlFileSink->SetModel(spModel);
        spMdlFileSink->Render();
      }
    }
    // END OF MAIN LOOP OVER LIGAND RECORDS
    ////////////////////////////////////////////////////
  } catch (Error &e) {
    std::cout << e.what() << std::endl;
  } catch (...) {
    std::cout << "Unknown exception" << std::endl;
  }

  _RBTOBJECTCOUNTER_DUMP_(std::cout)

  return 0;
}
