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

// Calculates vdW grids for use by VdwGridSF scoring function class

#include "rxdock/BiMolWorkSpace.h"
#include "rxdock/PRMFactory.h"
#include "rxdock/ParameterFileSource.h"
#include "rxdock/RealGrid.h"
#include "rxdock/SFFactory.h"
#include "rxdock/TriposAtomType.h"
#include <cstring>
#include <fstream>
#include <iomanip>

using namespace rxdock;

namespace rxdock {

const std::string _ROOT_SF = "SCORE";

// Creates list of probe models
// NOTE: MUST BE IN ORDER OF ASCENDING RADII
ModelList CreateProbes() {
  ModelList probes;
  TriposAtomTypeList atomTypes;
  atomTypes.push_back(TriposAtomType::UNDEFINED);
  atomTypes.push_back(TriposAtomType::Br);
  atomTypes.push_back(TriposAtomType::C_cat);
  atomTypes.push_back(TriposAtomType::C_1);
  atomTypes.push_back(TriposAtomType::C_1_H1);
  atomTypes.push_back(TriposAtomType::C_2);
  atomTypes.push_back(TriposAtomType::C_2_H1);
  atomTypes.push_back(TriposAtomType::C_2_H2);
  atomTypes.push_back(TriposAtomType::C_3);
  atomTypes.push_back(TriposAtomType::C_3_H1);
  atomTypes.push_back(TriposAtomType::C_3_H2);
  atomTypes.push_back(TriposAtomType::C_3_H3);
  atomTypes.push_back(TriposAtomType::C_ar);
  atomTypes.push_back(TriposAtomType::C_ar_H1);
  atomTypes.push_back(TriposAtomType::Cl);
  atomTypes.push_back(TriposAtomType::F);
  atomTypes.push_back(TriposAtomType::H);
  atomTypes.push_back(TriposAtomType::H_P);
  atomTypes.push_back(TriposAtomType::I);
  atomTypes.push_back(TriposAtomType::N_1);
  atomTypes.push_back(TriposAtomType::N_2);
  atomTypes.push_back(TriposAtomType::N_3);
  atomTypes.push_back(TriposAtomType::N_4);
  atomTypes.push_back(TriposAtomType::N_am);
  atomTypes.push_back(TriposAtomType::N_ar);
  atomTypes.push_back(TriposAtomType::N_pl3);
  atomTypes.push_back(TriposAtomType::O_2);
  atomTypes.push_back(TriposAtomType::O_3);
  atomTypes.push_back(TriposAtomType::O_co2);
  atomTypes.push_back(TriposAtomType::P_3);
  atomTypes.push_back(TriposAtomType::S_2);
  atomTypes.push_back(TriposAtomType::S_3);
  atomTypes.push_back(TriposAtomType::S_o);
  atomTypes.push_back(TriposAtomType::S_o2);
  for (TriposAtomTypeListConstIter iter = atomTypes.begin();
       iter != atomTypes.end(); iter++) {
    AtomList atomList;
    BondList bondList;
    AtomPtr spAtom(new Atom(1));
    spAtom->SetTriposType(*iter);
    spAtom->SetAtomicMass(
        12.0); // Mass is irrelevant as we are not rotating models around COM
    atomList.push_back(spAtom);
    probes.push_back(new Model(atomList, bondList));
  }
  return probes;
}

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
  std::string strSuffix(".json");
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
        << "Options:\t-o<OutputSuffix> - suffix for grid (.json IS required)"
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
    BiMolWorkSpacePtr spWS(new BiMolWorkSpace());
    // Set the workspace name to the root of the receptor .prm filename
    std::vector<std::string> componentList =
        ConvertDelimitedStringToList(strReceptorPrmFile, ".");
    std::string wsName = componentList.front();
    spWS->SetName(wsName);

    // Read the receptor parameter file
    ParameterFileSourcePtr spRecepPrmSource(new ParameterFileSource(
        GetDataFileName("data/receptors", strReceptorPrmFile)));
    std::cout << std::endl
              << "RECEPTOR:" << std::endl
              << spRecepPrmSource->GetFileName() << std::endl
              << spRecepPrmSource->GetTitle() << std::endl;

    // Read the scoring function file
    ParameterFileSourcePtr spSFSource(
        new ParameterFileSource(GetDataFileName("data/sf", strSFFile)));
    SFFactoryPtr spSFFactory(
        new SFFactory()); // Factory class for scoring functions
    SFAggPtr spSF(spSFFactory->CreateAggFromFile(
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
    PRMFactory prmFactory(spRecepPrmSource);
    ModelPtr spReceptor = prmFactory.CreateReceptor();
    // Trap multiple receptor conformations here: this SF does not support them
    // yet
    bool bEnsemble = (spReceptor->GetNumSavedCoords() > 1);
    if (bEnsemble) {
      std::string message(
          "rbcalcgrid does not support multiple receptor conformations yet");
      throw InvalidRequest(_WHERE_, message);
    }

    // Read docking site from file and register with workspace
    std::string strDockingSiteFile = spWS->GetName() + ".as";
    std::string strInputFile =
        GetDataFileName("data/grids", strDockingSiteFile);
    std::ifstream inputFile(strInputFile.c_str());
    json siteData;
    inputFile >> siteData;
    inputFile.close();
    DockingSitePtr spDS(new DockingSite(siteData.at("docking-site")));
    spWS->SetDockingSite(spDS);

    // Register receptor with workspace
    spWS->SetReceptor(spReceptor);
    std::cout << std::endl
              << "DOCKING SITE" << std::endl
              << (*spDS) << std::endl;

    // Create a grid covering the docking site, plus user-defined border
    Coord minCoord = spDS->GetMinCoord() - border;
    Coord maxCoord = spDS->GetMaxCoord() + border;
    Vector recepExtent = maxCoord - minCoord;
    Vector gridStep(gs, gs, gs);
    Eigen::Vector3d nXYZ = recepExtent.xyz.array() / gridStep.xyz.array();
    unsigned int nX = static_cast<unsigned int>(nXYZ(0)) + 1;
    unsigned int nY = static_cast<unsigned int>(nXYZ(1)) + 1;
    unsigned int nZ = static_cast<unsigned int>(nXYZ(2)) + 1;
    std::cout << "Constructing grid of size " << nX << " x " << nY << " x "
              << nZ << std::endl;
    RealGridPtr spGrid(new RealGrid(minCoord, gridStep, nX, nY, nZ));
    float *gridData = spGrid->GetGridData();

    // Create probes
    ModelList probes = CreateProbes();

    // Open output file
    std::string strOutputFile(spWS->GetName() + strSuffix);
    std::ofstream ostr(strOutputFile.c_str(),
                       std::ios_base::out | std::ios_base::trunc);

    json gridList;

    // Store regular pointers to avoid smart pointer dereferencing overheads
    RealGrid *pGrid(spGrid);
    SFAgg *pSF(spSF);
    TriposAtomType triposType;
    // Main loop over each probe model
    for (ModelListConstIter mIter = probes.begin(); mIter != probes.end();
         mIter++) {
      ModelPtr spLigand(*mIter);
      Atom *pAtom = spLigand->GetAtomList().front();
      TriposAtomType::eType atomType = pAtom->GetTriposType();
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
      json grid{{"tripos-type", strType}, {"real-grid", *pGrid}};
      gridList.push_back(grid);
    }
    json vdwGrids;
    vdwGrids["vdw-grids"] = gridList;
    ostr << vdwGrids;
    ostr.close();
  } catch (Error &e) {
    std::cout << e.what() << std::endl;
  } catch (...) {
    std::cout << "Unknown exception" << std::endl;
  }

  _RBTOBJECTCOUNTER_DUMP_(std::cout)

  return 0;
}
