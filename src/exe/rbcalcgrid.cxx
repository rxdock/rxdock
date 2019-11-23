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

// Calculates vdW grids for use by RbtVdwGridSF scoring function class

#include "RbtBiMolWorkSpace.h"
#include "RbtPRMFactory.h"
#include "RbtParameterFileSource.h"
#include "RbtRealGrid.h"
#include "RbtSFFactory.h"
#include "RbtTriposAtomType.h"
#include <cstring>
#include <fstream>
#include <iomanip>

const std::string EXEVERSION =
    " ($Id: //depot/dev/client3/rdock/2013.1/src/exe/rbcalcgrid.cxx#3 $)";
const std::string _ROOT_SF = "SCORE";

// Creates list of probe models
// NOTE: MUST BE IN ORDER OF ASCENDING RADII
RbtModelList CreateProbes() {
  RbtModelList probes;
  RbtTriposAtomTypeList atomTypes;
  atomTypes.push_back(RbtTriposAtomType::UNDEFINED);
  atomTypes.push_back(RbtTriposAtomType::Br);
  atomTypes.push_back(RbtTriposAtomType::C_cat);
  atomTypes.push_back(RbtTriposAtomType::C_1);
  atomTypes.push_back(RbtTriposAtomType::C_1_H1);
  atomTypes.push_back(RbtTriposAtomType::C_2);
  atomTypes.push_back(RbtTriposAtomType::C_2_H1);
  atomTypes.push_back(RbtTriposAtomType::C_2_H2);
  atomTypes.push_back(RbtTriposAtomType::C_3);
  atomTypes.push_back(RbtTriposAtomType::C_3_H1);
  atomTypes.push_back(RbtTriposAtomType::C_3_H2);
  atomTypes.push_back(RbtTriposAtomType::C_3_H3);
  atomTypes.push_back(RbtTriposAtomType::C_ar);
  atomTypes.push_back(RbtTriposAtomType::C_ar_H1);
  atomTypes.push_back(RbtTriposAtomType::Cl);
  atomTypes.push_back(RbtTriposAtomType::F);
  atomTypes.push_back(RbtTriposAtomType::H);
  atomTypes.push_back(RbtTriposAtomType::H_P);
  atomTypes.push_back(RbtTriposAtomType::I);
  atomTypes.push_back(RbtTriposAtomType::N_1);
  atomTypes.push_back(RbtTriposAtomType::N_2);
  atomTypes.push_back(RbtTriposAtomType::N_3);
  atomTypes.push_back(RbtTriposAtomType::N_4);
  atomTypes.push_back(RbtTriposAtomType::N_am);
  atomTypes.push_back(RbtTriposAtomType::N_ar);
  atomTypes.push_back(RbtTriposAtomType::N_pl3);
  atomTypes.push_back(RbtTriposAtomType::O_2);
  atomTypes.push_back(RbtTriposAtomType::O_3);
  atomTypes.push_back(RbtTriposAtomType::O_co2);
  atomTypes.push_back(RbtTriposAtomType::P_3);
  atomTypes.push_back(RbtTriposAtomType::S_2);
  atomTypes.push_back(RbtTriposAtomType::S_3);
  atomTypes.push_back(RbtTriposAtomType::S_o);
  atomTypes.push_back(RbtTriposAtomType::S_o2);
  for (RbtTriposAtomTypeListConstIter iter = atomTypes.begin();
       iter != atomTypes.end(); iter++) {
    RbtAtomList atomList;
    RbtBondList bondList;
    RbtAtomPtr spAtom(new RbtAtom(1));
    spAtom->SetTriposType(*iter);
    spAtom->SetAtomicMass(
        12.0); // Mass is irrelevant as we are not rotating models around COM
    atomList.push_back(spAtom);
    probes.push_back(new RbtModel(atomList, bondList));
  }
  return probes;
}

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
  Rbt::PrintStdHeader(std::cout, strExeName + EXEVERSION);

  // Command line arguments and default values
  std::string strSuffix(".grd");
  std::string strReceptorPrmFile;             // Receptor param file
  std::string strSFFile("calcgrid_attr.prm"); // Scoring function file
  double gs(0.5);                             // grid step
  double border(1.0);                         // grid border around docking site

  // Brief help message
  if (argc == 1) {
    std::cout << std::endl
              << "rbcalcgrid - calculates vdw grids for each atom type"
              << std::endl;
    std::cout << std::endl
              << "Usage:\trbcalcgrid -o<OutputRoot> -r<ReceptorPrmFile> "
                 "-p<SFPrmFile> [-g<GridStep>]"
              << std::endl;
    std::cout
        << std::endl
        << "Options:\t-o<OutputSuffix> - suffix for grid (.grd IS required)"
        << std::endl;
    std::cout
        << "\t\t-r<ReceptorPrmFile> - receptor param file (contains active "
           "site params)"
        << std::endl;
    std::cout << "\t\t-p<SFPrmFile> - scoring function param file (either "
                 "calcgrid_vdw1.prm or calcgrid_vdw5.prm)"
              << std::endl;
    std::cout << "\t\t-g<GridStep> - grid step (default=0.5A)" << std::endl;
    std::cout
        << "\t\t-b<Border> - grid border around docking site (default=1.0A)"
        << std::endl;
    return 1;
  }

  // Check command line arguments
  std::cout << std::endl << "Command line args:" << std::endl;
  for (int iarg = 1; iarg < argc; iarg++) {
    std::cout << argv[iarg];
    std::string strArg(argv[iarg]);
    if (strArg.find("-o") == 0)
      strSuffix = strArg.substr(2);
    else if (strArg.find("-r") == 0)
      strReceptorPrmFile = strArg.substr(2);
    else if (strArg.find("-p") == 0)
      strSFFile = strArg.substr(2);
    else if (strArg.find("-g") == 0) {
      std::string strGridStep = strArg.substr(2);
      gs = std::atof(strGridStep.c_str());
    } else if (strArg.find("-b") == 0) {
      std::string strBorder = strArg.substr(2);
      border = std::atof(strBorder.c_str());
    } else {
      std::cout << " ** INVALID ARGUMENT" << std::endl;
      return 1;
    }
    std::cout << std::endl;
  }

  std::cout << std::endl;

  try {
    // Create a bimolecular workspace
    RbtBiMolWorkSpacePtr spWS(new RbtBiMolWorkSpace());
    // Set the workspace name to the root of the receptor .prm filename
    std::vector<std::string> componentList =
        Rbt::ConvertDelimitedStringToList(strReceptorPrmFile, ".");
    std::string wsName = componentList.front();
    spWS->SetName(wsName);

    // Read the receptor parameter file
    RbtParameterFileSourcePtr spRecepPrmSource(new RbtParameterFileSource(
        Rbt::GetRbtFileName("data/receptors", strReceptorPrmFile)));
    std::cout << std::endl
              << "RECEPTOR:" << std::endl
              << spRecepPrmSource->GetFileName() << std::endl
              << spRecepPrmSource->GetTitle() << std::endl;

    // Read the scoring function file
    RbtParameterFileSourcePtr spSFSource(
        new RbtParameterFileSource(Rbt::GetRbtFileName("data/sf", strSFFile)));
    RbtSFFactoryPtr spSFFactory(
        new RbtSFFactory()); // Factory class for scoring functions
    RbtSFAggPtr spSF(spSFFactory->CreateAggFromFile(
        spSFSource, _ROOT_SF)); // Root SF aggregate

    // Register the scoring function with the workspace
    // Dump details to std::cout
    spWS->SetSF(spSF);
    std::cout << std::endl
              << "SCORING FUNCTION DETAILS:" << std::endl
              << *spSF << std::endl;

    // Create the receptor model from the file names in the receptor parameter
    // file
    spRecepPrmSource->SetSection();
    RbtPRMFactory prmFactory(spRecepPrmSource);
    RbtModelPtr spReceptor = prmFactory.CreateReceptor();
    // Trap multiple receptor conformations here: this SF does not support them
    // yet
    bool bEnsemble = (spReceptor->GetNumSavedCoords() > 1);
    if (bEnsemble) {
      std::string message(
          "rbcalcgrid does not support multiple receptor conformations yet");
      throw RbtInvalidRequest(_WHERE_, message);
    }

    // Read docking site from file and register with workspace
    std::string strASFile = spWS->GetName() + ".as";
    std::string strInputFile = Rbt::GetRbtFileName("data/grids", strASFile);
    // DM 26 Sep 2000 - std::ios_base::binary is invalid with IRIX CC
#if defined(__sgi) && !defined(__GNUC__)
    std::ifstream istr(strInputFile.c_str(), std::ios_base::in);
#else
    std::ifstream istr(strInputFile.c_str(),
                       std::ios_base::in | std::ios_base::binary);
#endif
    RbtDockingSitePtr spDS(new RbtDockingSite(istr));
    istr.close();
    spWS->SetDockingSite(spDS);

    // Register receptor with workspace
    spWS->SetReceptor(spReceptor);
    std::cout << std::endl
              << "DOCKING SITE" << std::endl
              << (*spDS) << std::endl;

    // Create a grid covering the docking site, plus user-defined border
    RbtCoord minCoord = spDS->GetMinCoord() - border;
    RbtCoord maxCoord = spDS->GetMaxCoord() + border;
    RbtVector recepExtent = maxCoord - minCoord;
    RbtVector gridStep(gs, gs, gs);
    unsigned int nX = int(recepExtent.x / gridStep.x) + 1;
    unsigned int nY = int(recepExtent.y / gridStep.y) + 1;
    unsigned int nZ = int(recepExtent.z / gridStep.z) + 1;
    std::cout << "Constructing grid of size " << nX << " x " << nY << " x "
              << nZ << std::endl;
    RbtRealGridPtr spGrid(new RbtRealGrid(minCoord, gridStep, nX, nY, nZ));
    float *gridData = spGrid->GetGridData();

    // Create probes
    RbtModelList probes = CreateProbes();

    // Open output file
    std::string strOutputFile(spWS->GetName() + strSuffix);
#if defined(__sgi) && !defined(__GNUC__)
    std::ofstream ostr(strOutputFile.c_str(),
                       std::ios_base::out | std::ios_base::trunc);
#else
    std::ofstream ostr(strOutputFile.c_str(), std::ios_base::out |
                                                  std::ios_base::binary |
                                                  std::ios_base::trunc);
#endif
    // Write header string (RbtVdwGridSF)
    const char *const header = "RbtVdwGridSF";
    int length = strlen(header);
    Rbt::WriteWithThrow(ostr, (const char *)&length, sizeof(length));
    Rbt::WriteWithThrow(ostr, header, length);
    // Write number of grids
    int nGrids = probes.size();
    Rbt::WriteWithThrow(ostr, (const char *)&nGrids, sizeof(nGrids));

    // Store regular pointers to avoid smart pointer dereferencing overheads
    RbtRealGrid *pGrid(spGrid);
    RbtSFAgg *pSF(spSF);
    RbtTriposAtomType triposType;
    // Main loop over each probe model
    for (RbtModelListConstIter mIter = probes.begin(); mIter != probes.end();
         mIter++) {
      RbtModelPtr spLigand(*mIter);
      RbtAtom *pAtom = spLigand->GetAtomList().front();
      RbtTriposAtomType::eType atomType = pAtom->GetTriposType();
      std::string strType = triposType.Type2Str(atomType);
      std::cout << "Atom type=" << strType << std::endl;
      // Register ligand with workspace
      spWS->SetLigand(spLigand);
      pGrid->SetAllValues(0.0);
      // Loop over all grid coords and calculate the score at each position
      for (unsigned int i = 0; i < spGrid->GetN(); i++) {
        pAtom->SetCoords(pGrid->GetCoord(i));
        gridData[i] = pSF->Score();
      }
      // Write the atom type string to the grid file, before the grid itself
      const char *const szType = strType.c_str();
      int l = strlen(szType);
      Rbt::WriteWithThrow(ostr, (const char *)&l, sizeof(l));
      Rbt::WriteWithThrow(ostr, szType, l);
      pGrid->Write(ostr);
    }
    ostr.close();
  } catch (RbtError &e) {
    std::cout << e << std::endl;
  } catch (...) {
    std::cout << "Unknown exception" << std::endl;
  }

  _RBTOBJECTCOUNTER_DUMP_(std::cout)

  return 0;
}
