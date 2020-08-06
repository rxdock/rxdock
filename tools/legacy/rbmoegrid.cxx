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

#include <cxxopts.hpp>
#include <fstream>
#include <iomanip>

#include "rxdock/BiMolWorkSpace.h"
#include "rxdock/MOEGrid.h"
#include "rxdock/PRMFactory.h"
#include "rxdock/ParameterFileSource.h"
#include "rxdock/RealGrid.h"
#include "rxdock/SFFactory.h"
#include "rxdock/TriposAtomType.h"

using namespace rxdock;

namespace rxdock {

const std::string _ROOT_SF = "SCORE";

// Creates list of probe models
// NOTE: MUST BE IN ORDER OF ASCENDING RADII
ModelList CreateProbes(std::string anAtomTypeStr) {
  ModelList probes;
  TriposAtomTypeList atomTypes;
  TriposAtomType tAtomType;
  atomTypes.push_back(TriposAtomType::UNDEFINED);
  //	keep is as a vector though there is only one element in there
  for (TriposAtomTypeListConstIter iter = atomTypes.begin();
       iter != atomTypes.end(); iter++) {
    AtomList atomList;
    BondList bondList;
    AtomPtr spAtom(new Atom(1));
    spAtom->SetTriposType(tAtomType.Str2Type(anAtomTypeStr));
    if (spAtom->GetTriposType() == TriposAtomType::UNDEFINED)
      throw BadArgument(_WHERE_,
                        "Undefined atom type (try standard Tripos types)");
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

  cxxopts::Options options(
      strExeName, "rbmoegrid - calculates grids for a given atom type");

  // Command line arguments and default values
  cxxopts::OptionAdder adder = options.add_options();
  adder("o,output", "output file name prefix (.grd is suffixed)",
        cxxopts::value<std::string>()->default_value("moegrid"));
  adder("r,receptor-param",
        "receptor param file(s) (contain active site params, can be specified "
        "multiple times)",
        cxxopts::value<std::vector<std::string>>());
  adder("p,sf-param", "scoring function param file",
        cxxopts::value<std::string>()->default_value("calcgrid_vdw.prm"));
  adder("g,grid-step", "grid step (in angstrom)",
        cxxopts::value<double>()->default_value("0.5"));
  adder("b,grid-border", "grid border around docking site (in angstrom)",
        cxxopts::value<double>()->default_value("1.0"));
  adder("t,tripos-atom-type", "Tripos atom type",
        cxxopts::value<std::string>()->default_value("C.3"));
  adder("h,help", "Print help");

  // start processing receptors
  try {
    auto result = options.parse(argc, argv);

    if (result.count("h")) {
      std::cout << options.help() << std::endl;
      return 0;
    }

    // Command line arguments and default values
    std::string strSuffix(".grd"); // file out suffix
    std::string strOutfName =
        result["o"].as<std::string>(); // default output filename
    std::vector<std::string> strReceptorPrmFiles;
    if (result.count("r")) {
      strReceptorPrmFiles =
          result["r"].as<std::vector<std::string>>(); // Receptor param files
                                                      // (you can have more!)
    }
    std::string strSFFile =
        result["p"].as<std::string>();        // Scoring function file
    double gs = result["g"].as<double>();     // grid step
    double border = result["b"].as<double>(); // grid border around docking site
    std::string strAtomType = result["t"].as<std::string>();

    if (strReceptorPrmFiles.empty() || strReceptorPrmFiles[0].empty()) {
      std::cout << "Missing receptor parameter file name(s)" << std::endl;
      return 1;
    }

    // generate output file name
    std::string strOutputFile(strOutfName + strSuffix);
    MOEGrid theMOEGrid;
    theMOEGrid.SetOutputFileName(strOutputFile);
    // read each receptors one-by-one to get a common set of grid extents and
    // origin
    theMOEGrid.CalculateCommonExtents(strReceptorPrmFiles);
    Coord minCoord = theMOEGrid.GetMinExtents() - border;
    Coord maxCoord = theMOEGrid.GetMaxExtents() + border;
    Vector gridStep(gs, gs, gs);
    Vector recepExtent = maxCoord - minCoord;
    Eigen::Vector3d nXYZ = recepExtent.xyz.array() / gridStep.xyz.array();
    unsigned int nX = static_cast<unsigned int>(nXYZ(0)) + 1;
    unsigned int nY = static_cast<unsigned int>(nXYZ(1)) + 1;
    unsigned int nZ = static_cast<unsigned int>(nXYZ(2)) + 1;
    std::cout << "Constructing grid of size " << nX << " x " << nY << " x "
              << nZ << std::endl;
    std::vector<std::string>::iterator strReceptorPrmFilesIter;
    for (strReceptorPrmFilesIter = strReceptorPrmFiles.begin();
         strReceptorPrmFilesIter != strReceptorPrmFiles.end();
         ++strReceptorPrmFilesIter) {
      // Create a bimolecular workspace
      BiMolWorkSpacePtr spWS(new BiMolWorkSpace());
      // Set the workspace name to the root of the receptor .prm filename
      std::vector<std::string> componentList =
          ConvertDelimitedStringToList((*strReceptorPrmFilesIter), ".");
      std::string wsName = componentList.front();
      spWS->SetName(wsName);

      // Read the receptor parameter file
      ParameterFileSourcePtr spRecepPrmSource(new ParameterFileSource(
          GetDataFileName("data/receptors", (*strReceptorPrmFilesIter))));
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

      // Read docking site from file and register with workspace
      std::string strDockingSiteFile = spWS->GetName() + ".as";
      std::string strInputFile =
          GetDataFileName("data/grids", strDockingSiteFile);
      std::ifstream inputFile(strInputFile);
      json siteData;
      inputFile >> siteData;
      inputFile.close();
      DockingSitePtr spDS(new DockingSite(siteData["docking-site"]));
      spWS->SetDockingSite(spDS);

      // Register receptor with workspace
      spWS->SetReceptor(spReceptor);
      std::cout << std::endl
                << "DOCKING SITE" << std::endl
                << (*spDS) << std::endl;

      // Create a grid covering the docking site, plus user-defined border
      // Coord minCoord = spDS->GetMinCoord()-border;
      // Coord maxCoord = spDS->GetMaxCoord()+border;
      // Vector recepExtent = maxCoord-minCoord;
      // Vector gridStep(gs,gs,gs);
      // Eigen::Vector3d nXYZ = recepExtent.xyz.array() / gridStep.xyz.array();
      // UInt nX = static_cast<unsigned int>(nXYZ(0)) + 1;
      // UInt nY = static_cast<unsigned int>(nXYZ(1)) + 1;
      // UInt nZ = static_cast<unsigned int>(nXYZ(2)) + 1;
      // std::cout << "Constructing grid of size " << nX << " x " << nY << " x "
      // << nZ << std::endl;
      RealGridPtr spGrid(new RealGrid(minCoord, gridStep, nX, nY, nZ));

      // Create probes
      ModelList probes = CreateProbes(strAtomType);

      // write the grid into a MOE grid file
      std::cout << "==================================" << std::endl;
      std::cout << "Generating MOE grid " << std::endl;
      // grid extents are in minCoord and maxCoord
      std::vector<double> grid_origin(
          minCoord.xyz.data(), minCoord.xyz.data() + minCoord.xyz.size());
      std::vector<double> grid_extents(
          maxCoord.xyz.data(), maxCoord.xyz.data() + maxCoord.xyz.size());
      // create MOE grid shape vector
      MOEGridShape theShape(grid_origin, grid_extents, gs);
      // MOE grid data vector
      MOEGridData theData;
      // reserve memory for data
      theData.reserve(theShape.GetDataSize());
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
        // Loop over all grid coords and calculate the score at each position
        for (unsigned int i = 0; i < spGrid->GetN(); i++) {
          pAtom->SetCoords(pGrid->GetCoord(i));
          MOEGridPoint p;
          p.SetValue(pSF->Score());
          theData.push_back(p);
        }
      }
      // create grid object from shape and data
      theMOEGrid.SetShape(theShape);
      theMOEGrid.SetData(theData);
      std::cout << "writing grid... ";
      // if it is the very first, then override
      if (strReceptorPrmFilesIter == strReceptorPrmFiles.begin()) {
        theMOEGrid.WriteGrid(std::ios_base::trunc);
        // otherwise append
      } else {
        theMOEGrid.WriteGrid(std::ios_base::app);
      }
      std::cout << "done." << std::endl;
    }
  } catch (const cxxopts::OptionException &e) {
    std::cout << "Error parsing options: " << e.what() << std::endl;
    return 1;
  } catch (Error &e) {
    std::cout << e.what() << std::endl;
  } catch (...) {
    std::cout << "Unknown exception" << std::endl;
  }

  _RBTOBJECTCOUNTER_DUMP_(std::cout)

  return 0;
}
