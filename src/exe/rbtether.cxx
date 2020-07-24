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

// Filters a ligand library with a SMARTS query and writes down the
// corresponding tethered atoms in the sd files

#include <iomanip>
#include <sstream>

#include "BiMolWorkSpace.h"
#include "MdlFileSink.h"
#include "MdlFileSource.h"
#include "PRMFactory.h"
#include "ParameterFileSource.h"
#include "Smarts.h"

using namespace rxdock;

namespace rxdock {

void print_atoms(AtomList &atoms, std::ostringstream &ost);

} // namespace rxdock

/////////////////////////////////////////////////////////////////////
// MAIN PROGRAM STARTS HERE
/////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[]) {
  std::cout.setf(std::ios_base::left, std::ios_base::adjustfield);

  // Strip off the path to the executable, leaving just the file name
  std::string strExeName(argv[0]);
  std::string::size_type i = strExeName.rfind("/");
  if (i != std::string::npos)
    strExeName.erase(0, i + 1);

  // Print a standard header
  PrintStdHeader(std::cout, strExeName);

  // Command line arguments and default values
  std::string strLigandMdlFile;
  bool bOutput(false);
  std::string strRunName;
  std::string strReferenceSDFile; // Receptor param file
  std::string strQuery;           // Smarts query
  bool bPosIonise(false);
  bool bNegIonise(false);
  bool bImplH(true);

  // Brief help message
  if (argc == 1) {
    std::cout << std::endl
              << "rbtether - marks the tethered atoms for use in docking"
              << std::endl;
    std::cout << std::endl
              << "Usage:\trbtether -i<InputSDFile> -o<OutputRoot>"
              << "-r<ReferenceSDFile> -s<Smarts Query>"
              << " [-ap] [-an] [-allH]" << std::endl;
    std::cout << std::endl
              << "Options:\t-i<InputSDFile> - input ligand SD file"
              << std::endl;
    std::cout << "\t\t-o<OutputRoot> - root name for output SD file(s)"
              << "(.sd suffix not required)" << std::endl;
    std::cout << "\t\t-r<ReferenceSDFile> - reference SD file "
              << "(with the correct coordinates for the tethered atoms)"
              << std::endl;
    std::cout << "\t\t-s<Smarts Query> - SMARTS query "
              << "(specifying the atoms to be tethered)" << std::endl;
    std::cout << "\t\t-ap - protonate all neutral amines, guanidines, "
              << "imidazoles (default=disabled)" << std::endl;
    std::cout
        << "\t\t-an - deprotonate all carboxylic, sulphur and phosphorous "
        << "acid groups (default=disabled)" << std::endl;
    std::cout << "\t\t-allH - read all hydrogens present (default=polar "
              << "hydrogens only)" << std::endl;
    return 1;
  }

  // Check command line arguments
  std::cout << std::endl << "Command line args:" << std::endl;
  for (int iarg = 1; iarg < argc; iarg++) {
    std::cout << argv[iarg];
    std::string strArg(argv[iarg]);
    if (strArg.find("-i") == 0)
      strLigandMdlFile = strArg.substr(2);
    else if (strArg.find("-o") == 0) {
      strRunName = strArg.substr(2);
      bOutput = true;
    } else if (strArg.find("-r") == 0)
      strReferenceSDFile = strArg.substr(2);
    else if (strArg.find("-s") == 0)
      strQuery = strArg.substr(2);
    else if (strArg.find("-ap") == 0)
      bPosIonise = true;
    else if (strArg.find("-an") == 0)
      bNegIonise = true;
    else if (strArg.find("-allH") == 0)
      bImplH = false;
    else {
      std::cout << " ** INVALID ARGUMENT" << std::endl;
      return 1;
    }
    std::cout << std::endl;
  }

  std::cout << std::endl;

  // DM 20 Apr 1999 - set the auto-ionise flags
  if (bPosIonise)
    std::cout << "Automatically protonating positive ionisable groups "
              << "(amines, imidazoles, guanidines)" << std::endl;
  if (bNegIonise)
    std::cout << "Automatically deprotonating negative ionisable groups "
              << "(carboxylic acids, phosphates, sulphates, sulphonates)"
              << std::endl;
  if (bImplH)
    std::cout << "Reading polar hydrogens only from ligand SD file"
              << std::endl;
  else
    std::cout << "Reading all hydrogens from ligand SD file" << std::endl;

  try {
    // Read the reference SD file
    MolecularFileSourcePtr spReferenceSD(
        new MdlFileSource(strReferenceSDFile, bPosIonise, bNegIonise, bImplH));
    ModelPtr spReferenceMdl(new Model(spReferenceSD));
    std::string strSmiles;
    AtomListList tetherAtoms =
        DT::QueryModel(spReferenceMdl, strQuery, strSmiles);
    std::cout << std::endl
              << "REFERENCE SD FILE: " << strReferenceSDFile << std::endl;
    if (strQuery.empty()) {
      std::cout
          << "No SMARTS query defined - tethering to whole of reference ligand"
          << std::endl;
      strQuery = strSmiles;
    }
    std::cout << "SMILES: " << strSmiles << std::endl;
    std::cout << "SMARTS: " << strQuery << std::endl;
    for (AtomListListIter alli = tetherAtoms.begin(); alli != tetherAtoms.end();
         alli++) {
      std::cout << "Path: ";
      for (AtomListIter ali = alli->begin(); ali != alli->end(); ali++) {
        std::cout << (*ali)->GetName() << "\t";
      }
      std::cout << std::endl;
    }
    if ((argc == 2) && !strReferenceSDFile.empty()) {
      // smile expression found for the reference file. Nothing
      // else to be done
      std::exit(0);
    }
    if (tetherAtoms.empty()) {
      throw Error(_WHERE_, "SMARTS query should have at least one match in "
                           "the reference SD file");
    }
    if (tetherAtoms.size() > 1) {
      // throw Error(_WHERE_,
      std::cout << "More than one SMARTS match found with the reference SD file"
                << std::endl;
      std::cout << "Will only tether to the first matching path..."
                << std::endl;
    }

    // DM 1 Jul 2002 - calculate principal axes and COM of tethered atoms
    // for the reference structure. Use to prealign each ligand with the
    // reference.
    AtomList tetheredAtomList(*(tetherAtoms.begin()));
    PrincipalAxes refAxes = GetPrincipalAxesOfAtoms(tetheredAtomList);

    // Prepare the SD file sink for saving the docked conformations
    // for each ligand
    if (bOutput) {
      // writing down the reference sd file for the query
      std::ostringstream ost;
      print_atoms(tetheredAtomList, ost);
      Variant vTetherAtoms(ost.str());
      MolecularFileSinkPtr spRefMdlFileSink(
          new MdlFileSink(strRunName + "_reference.sd", spReferenceMdl));
      spReferenceMdl->SetDataValue("TETHERED ATOMS", vTetherAtoms);
      spRefMdlFileSink->Render();
    }

    // DM 18 May 1999
    // Variants describing the library version, exe version,
    // parameter file, and current directory
    // Will be stored in the ligand SD files
    Variant vLib(GetProduct() + "/" + GetProgramVersion());
    Variant vExe(strExeName + "/" + GetProgramVersion());
    Variant vDir(GetCurrentWorkingDirectory());
    ///////////////////////////////////
    // MAIN LOOP OVER LIGAND RECORDS
    ///////////////////////////////////
    // DM 20 Apr 1999 - add explicit bPosIonise and bNegIonise flags
    // to MdlFileSource constructor
    MolecularFileSinkPtr spMdlFileSink(
        new MdlFileSink(strRunName + ".sd", ModelPtr()));
    MolecularFileSourcePtr spMdlFileSource(
        new MdlFileSource(strLigandMdlFile, bPosIonise, bNegIonise, bImplH));
    for (int nRec = 1; spMdlFileSource->FileStatusOK();
         spMdlFileSource->NextRecord(), nRec++) {
      std::cout.setf(std::ios_base::left, std::ios_base::adjustfield);
      std::cout << std::endl
                << "**************************************************"
                << std::endl
                << "RECORD #" << nRec << std::endl;
      Error molStatus = spMdlFileSource->Status();
      if (!molStatus.isOK()) {
        std::cout << std::endl
                  << molStatus << std::endl
                  << "************************************************"
                  << std::endl;
        continue;
      }

      // DM 26 Jul 1999 - only read the largest segment
      //(guaranteed to be called H)
      spMdlFileSource->SetSegmentFilterMap(ConvertStringToSegmentMap("H"));

      if (spMdlFileSource->isDataFieldPresent("Name"))
        std::cout << "NAME:   " << spMdlFileSource->GetDataValue("Name")
                  << std::endl;
      if (spMdlFileSource->isDataFieldPresent("REG_Number"))
        std::cout << "REG_Num:" << spMdlFileSource->GetDataValue("REG_Number")
                  << std::endl;

      // Create the ligand model
      ModelPtr spLigand(new Model(spMdlFileSource));
      std::string strMolName = spLigand->GetName();

      tetherAtoms = DT::QueryModel(spLigand, strQuery, strSmiles);
      std::cout << "SMILES: " << strSmiles << std::endl;
      std::cout << "SMARTS: " << strQuery << std::endl;
      if (tetherAtoms.empty()) {
        std::cout << "No SMARTS match, structure not output" << std::endl;
        continue; // ligand filtered out
      }
      AtomList ligAtomList = spLigand->GetAtomList();
      for (AtomListListIter alli = tetherAtoms.begin();
           alli != tetherAtoms.end(); alli++) {
        std::cout << "Path: ";
        for (AtomListIter ali = alli->begin(); ali != alli->end(); ali++) {
          std::cout << (*ali)->GetName() << "\t";
        }
        std::cout << std::endl;
        /////////////////////////////////
        // DM 1 Jul 2002 - prealign each ligand with the reference tether
        // Calculate quat needed to overlay principal axes of tethered fragment
        // with reference tether
        PrincipalAxes prAxes = GetPrincipalAxesOfAtoms(*alli);
        Quat q = GetQuatFromAlignAxes(prAxes, refAxes);
        // Apply the transformation to the whole ligand
        // 1. Translate COM of tethered fragment to the origin
        // 2. Perform the rotation
        // 3. Translate COM of tethered fragment to reference tether COM
        std::for_each(ligAtomList.begin(), ligAtomList.end(),
                      TranslateAtom(-prAxes.com));
        std::for_each(ligAtomList.begin(), ligAtomList.end(),
                      RotateAtomUsingQuat(q));
        std::for_each(ligAtomList.begin(), ligAtomList.end(),
                      TranslateAtom(refAxes.com));
        ///////////////////////////////////

        std::ostringstream ost;
        print_atoms(*alli, ost);
        Variant vTetherAtoms(ost.str());
        // DM 18 May 1999 - store run info in model data
        // Clear any previous Rbt.* data fields
        spLigand->ClearAllDataFields("Rbt.");
        spLigand->ClearAllDataFields("TETHERED ATOMS");
        spLigand->SetDataValue("Rbt.Library", vLib);
        spLigand->SetDataValue("Rbt.Executable", vExe);
        spLigand->SetDataValue("Rbt.Current_Directory", vDir);
        spLigand->SetDataValue("TETHERED ATOMS", vTetherAtoms);
        spMdlFileSink->SetModel(spLigand);
        spMdlFileSink->Render();
      }
    }
    // END OF MAIN LOOP OVER LIGAND RECORDS
    ////////////////////////////////////////////////////
    std::cout << std::endl << "END OF RUN" << std::endl;
  } catch (Error &e) {
    std::cout << e << std::endl;
  } catch (...) {
    std::cout << "Unknown exception" << std::endl;
  }

  _RBTOBJECTCOUNTER_DUMP_(std::cout)

  return 0;
}

void rxdock::print_atoms(AtomList &atoms, std::ostringstream &ost) {
  ost.clear();
  for (int iter = 0; iter < atoms.size(); iter++) {
    ost << atoms[iter]->GetAtomId();
    //    if ((((iter + 1) % 35) == 0) || (iter == (atoms.size() - 1)))
    if (iter == (atoms.size() - 1))
      ost << "\n";
    else if (((iter + 1) % 35) == 0)
      ost << ",\n";
    else
      ost << ",";
  }
}
