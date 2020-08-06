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

// Standalone executable for generating docking site .json files for rbdock

#include <algorithm>
#include <cxxopts.hpp>
#include <iomanip>

#include "BiMolWorkSpace.h"
#include "CrdFileSink.h"
#include "DockingSite.h"
#include "PRMFactory.h"
#include "ParameterFileSource.h"
#include "PsfFileSink.h"
#include "SiteMapperFactory.h"

using namespace rxdock;

/////////////////////////////////////////////////////////////////////
// MAIN PROGRAM STARTS HERE
/////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[]) {
  // Handle obsolete arguments, if any
  for (int i = 0; i < argc; i++) {
    std::string opt = argv[i];
    if (opt == "-was" || opt == "-ras") {
      std::cout
          << "Options -was and -ras are no longer supported; use -W and -R "
             "(respectively) instead."
          << std::endl;
      return 2;
    }
  }

  // Strip off the path to the executable, leaving just the file name
  std::string strExeName(argv[0]);
  std::string::size_type i = strExeName.rfind("/");
  if (i != std::string::npos)
    strExeName.erase(0, i + 1);

  // Print a standard header
  PrintStdHeader(std::cout, strExeName);

  cxxopts::Options options(strExeName, "rbcavity - calculate docking cavities");

  // Command line arguments and default values
  cxxopts::OptionAdder adder = options.add_options();
  adder("r,receptor-param", "receptor param file (contains active site params)",
        cxxopts::value<std::string>());
  adder("W,write-docking-cavities",
        "write docking cavities (plus distance grid) to .json file");
  adder("R,read-docking-cavities",
        "read docking cavities (plus distance grid) from .json file");
  adder("d,write-insightii-grids",
        "dump InsightII grids for each cavity for visualisation");
  adder("v,write-psf-crd", "dump target PSF/CRD files for rDock Viewer");
  adder("l,list-atoms-dist",
        "list receptor atoms within specified distance of any cavity (in "
        "angstrom)",
        cxxopts::value<double>()->default_value("5.0"));
  adder("s,print-site", "print SITE descriptors (counts of exposed atoms)");
  adder(
      "b,border",
      "set the border around the cavities for the distance grid (in angstrom)",
      cxxopts::value<double>()->default_value("8.0"));
  adder("m,write-moe-grid", "write active site into a MOE grid");
  adder("h,help", "Print help");

  try {
    auto result = options.parse(argc, argv);

    if (result.count("h")) {
      std::cout << options.help() << std::endl;
      return 0;
    }

    std::string strReceptorPrmFile;
    if (result.count("r")) {
      strReceptorPrmFile = result["r"].as<std::string>();
    }
    bool bReadDS = result.count("R");  // If true, read Active Site from file
    bool bWriteAS = result.count("W"); // If true, write Active Site to file
    bool bDump =
        result.count("d"); // If true, dump cavity grids in Insight format
    bool bViewer =
        result.count("v"); // If true, dump PSF/CRD files for rDock Viewer
    bool bList =
        result.count("l"); // If true, list atoms within 'distance' of cavity
    bool bSite = result.count(
        "s"); // If true, print out "SITE" descriptors (counts of exposed atoms)
    bool bMOEgrid = result.count(
        "m"); // If true, create a MOE grid file for AS visualisation
    bool bBorderArg = result.count("b"); // If true, border was specified in the
                                         // command line
    double border =
        result["b"]
            .as<double>(); // Border to allow around cavities for distance grid
    double dist = result["l"].as<double>();

    // check for parameter file name
    /*if(prmFile) {
            strReceptorPrmFile	= prmFile;
    }*/
    if (strReceptorPrmFile.empty()) {
      std::cout << "Missing receptor parameter file name" << std::endl;
      return 1;
    }
    // writing command line arguments
    std::cout << "Command line arguments:" << std::endl;
    std::cout << "-r " << strReceptorPrmFile << std::endl;
    if (bList)
      std::cout << "-l " << dist << std::endl;
    if (bBorderArg)
      std::cout << "-b " << border << std::endl;
    if (bWriteAS)
      std::cout << "-W" << std::endl;
    if (bReadDS)
      std::cout << "-R" << std::endl;
    if (bMOEgrid)
      std::cout << "-m" << std::endl;
    if (bDump)
      std::cout << "-d" << std::endl;
    if (bSite)
      std::cout << "-s" << std::endl;
    if (bViewer)
      std::cout << "-v" << std::endl;

    std::cout.setf(std::ios_base::left, std::ios_base::adjustfield);

    // Create a bimolecular workspace
    BiMolWorkSpacePtr spWS(new BiMolWorkSpace());
    // Set the workspace name to the root of the receptor .prm filename
    std::vector<std::string> componentList =
        ConvertDelimitedStringToList(strReceptorPrmFile, ".");
    std::string wsName = componentList.front();
    spWS->SetName(wsName);

    // Read the protocol parameter file
    ParameterFileSourcePtr spRecepPrmSource(new ParameterFileSource(
        GetDataFileName("data/receptors", strReceptorPrmFile)));

    // Create the receptor model from the file names in the parameter file
    spRecepPrmSource->SetSection();
    PRMFactory prmFactory(spRecepPrmSource);
    ModelPtr spReceptor = prmFactory.CreateReceptor();

    DockingSitePtr spDockSite;
    std::string strDockingSiteFile = wsName + "-docking-site.json";

    // Either read the docking site from the .as file
    if (bReadDS) {
      std::string strInputFile =
          GetDataFileName("data/grids", strDockingSiteFile);
      std::ifstream inputFile(strInputFile.c_str());
      json siteData;
      inputFile >> siteData;
      inputFile.close();
      spDockSite = DockingSitePtr(new DockingSite(siteData.at("docking-site")));
    }
    // Or map the site using the prescribed mapping algorithm
    else {
      SiteMapperFactoryPtr spMapperFactory(new SiteMapperFactory());
      SiteMapperPtr spMapper =
          spMapperFactory->CreateFromFile(spRecepPrmSource, "MAPPER");
      spMapper->Register(spWS);
      spWS->SetReceptor(spReceptor);
      std::cout << *spMapper << std::endl;

      int nRI = spReceptor->GetNumSavedCoords() - 1;
      if (nRI == 0) {
        spDockSite = DockingSitePtr(new DockingSite((*spMapper)(), border));
      } else {
        CavityList allCavs;
        for (int i = 1; i <= nRI; i++) {
          spReceptor->RevertCoords(i);
          CavityList cavs((*spMapper)());
          std::copy(cavs.begin(), cavs.end(), std::back_inserter(allCavs));
        }
        spDockSite = DockingSitePtr(new DockingSite(allCavs, border));
      }
    }

    std::cout << std::endl
              << "DOCKING SITE" << std::endl
              << (*spDockSite) << std::endl;

    if (bWriteAS) {
      std::ofstream ostr(strDockingSiteFile.c_str(),
                         std::ios_base::out | std::ios_base::trunc);
      json dockingSite;
      dockingSite["docking-site"] = *spDockSite;
      ostr << dockingSite;
      ostr.close();
    }

    // Write PSF/CRD files to keep the rDock Viewer happy (it doesn't read MOL2
    // files yet)
    if (bViewer) {
      MolecularFileSinkPtr spRecepSink =
          new PsfFileSink(wsName + "_for_viewer.psf", spReceptor);
      std::cout << "Writing PSF file: " << spRecepSink->GetFileName()
                << std::endl;
      spRecepSink->Render();
      spRecepSink = new CrdFileSink(wsName + "_for_viewer.crd", spReceptor);
      std::cout << "Writing CRD file: " << spRecepSink->GetFileName()
                << std::endl;
      spRecepSink->Render();
    }

    // Write an ASCII InsightII grid file for each defined cavity
    if (bDump) {
      CavityList cavList = spDockSite->GetCavityList();
      for (unsigned int i = 0; i < cavList.size(); i++) {
        std::ostringstream filename;
        filename << wsName << "_cav" << i + 1 << ".grd";
        std::ofstream dumpFile(filename.str());
        if (dumpFile) {
          cavList[i]->GetGrid()->PrintInsightGrid(dumpFile);
          dumpFile.close();
        }
      }
    }
    // writing active site into MOE grid
    if (bMOEgrid) {
      std::cout << "MOE grid feature not yet implemented, sorry." << std::endl;
    }
    // List all receptor atoms within given distance of any cavity
    if (bList) {
      RealGridPtr spGrid = spDockSite->GetGrid();
      AtomList atomList =
          spDockSite->GetAtomList(spReceptor->GetAtomList(), 0.0, dist);
      std::cout << atomList.size() << " receptor atoms within " << dist
                << " A of any cavity" << std::endl;
      std::cout << std::endl << "DISTANCE,ATOM" << std::endl;
      for (AtomListConstIter iter = atomList.begin(); iter != atomList.end();
           iter++) {
        std::cout << spGrid->GetSmoothedValue((*iter)->GetCoords()) << "\t"
                  << **iter << std::endl;
      }
      std::cout << std::endl;
    }

    // DM 15 Jul 2002 - print out SITE descriptors
    // Use a crude measure of solvent accessibility - count #atoms within 4A of
    // each atom Use an empirical threshold to determine if atom is exposed or
    // not
    if (bSite) {
      double cavDist = 4.0; // Use a fixed definition of cavity atoms - all
                            // those within 4A of docking volume
      double neighbR = 4.0; // Sphere radius for counting nearest neighbours
      double threshold =
          15; // Definition of solvent exposed: neighbours < threshold
      // RealGridPtr spGrid = spDockSite->GetGrid();
      AtomList recepAtomList = spReceptor->GetAtomList();
      AtomList cavAtomList = spDockSite->GetAtomList(recepAtomList, cavDist);
      AtomList exposedAtomList; // The list of exposed cavity atoms
      std::cout << std::endl << "SOLVENT EXPOSED CAVITY ATOMS" << std::endl;
      std::cout << "1) Consider atoms within " << cavDist << "A of docking site"
                << std::endl;
      std::cout << "2) Determine #neighbours within " << neighbR
                << "A of each atom" << std::endl;
      std::cout << "3) If #neighbours < " << threshold << " => exposed"
                << std::endl;
      std::cout
          << "4) Calculate SITE* descriptors over exposed cavity atoms only"
          << std::endl;
      std::cout << std::endl << "ATOM NAME\t#NEIGHBOURS" << std::endl;

      // Get the list of solvent exposed cavity atoms
      for (AtomListConstIter iter = cavAtomList.begin();
           iter != cavAtomList.end(); iter++) {
        int nNeighb = GetNumAtomsWithPredicate(
            recepAtomList, isAtomInsideSphere((*iter)->GetCoords(), neighbR));
        nNeighb--;
        if (nNeighb < threshold) {
          std::cout << (*iter)->GetFullAtomName() << "\t" << nNeighb
                    << std::endl;
          exposedAtomList.push_back(*iter);
        }
      }

      // Total +ve and -ve charges
      double posChg(0.0);
      double negChg(0.0);
      for (AtomListConstIter iter = exposedAtomList.begin();
           iter != exposedAtomList.end(); iter++) {
        double chg = (*iter)->GetGroupCharge();
        if (chg > 0.0) {
          posChg += chg;
        } else if (chg < 0.0) {
          negChg += chg;
        }
      }

      // Atom type counts
      isHybridState_eq bIsArom(Atom::AROM);
      int nAtoms = exposedAtomList.size();
      int nLipoC =
          GetNumAtomsWithPredicate(exposedAtomList, isAtomLipophilic());
      int nArom = GetNumAtomsWithPredicate(exposedAtomList, bIsArom);
      int nNHBD = GetNumAtomsWithPredicate(exposedAtomList, isAtomHBondDonor());
      int nMetal = GetNumAtomsWithPredicate(exposedAtomList, isAtomMetal());
      int nGuan =
          GetNumAtomsWithPredicate(exposedAtomList, isAtomGuanidiniumCarbon());
      int nNHBA =
          GetNumAtomsWithPredicate(exposedAtomList, isAtomHBondAcceptor());

      // Cavity volume
      std::cout << std::endl
                << wsName << ",SITE_VOL," << spDockSite->GetVolume()
                << std::endl;
      // Atom type counts
      std::cout << wsName << ",SITE_NATOMS," << nAtoms << std::endl;
      std::cout << wsName << ",SITE_NLIPOC," << nLipoC << std::endl;
      std::cout << wsName << ",SITE_NAROMATOMS," << nArom << std::endl;
      std::cout << wsName << ",SITE_NHBD," << nNHBD << std::endl;
      std::cout << wsName << ",SITE_NMETAL," << nMetal << std::endl;
      std::cout << wsName << ",SITE_NGUAN," << nGuan << std::endl;
      std::cout << wsName << ",SITE_NHBA," << nNHBA << std::endl;
      // Atom type percentages
      std::cout << wsName << ",SITE_PERC_LIPOC," << 100.0 * nLipoC / nAtoms
                << std::endl;
      std::cout << wsName << ",SITE_PERC_AROMATOMS," << 100.0 * nArom / nAtoms
                << std::endl;
      std::cout << wsName << ",SITE_PERC_HBD," << 100.0 * nNHBD / nAtoms
                << std::endl;
      std::cout << wsName << ",SITE_PERC_METAL," << 100.0 * nMetal / nAtoms
                << std::endl;
      std::cout << wsName << ",SITE_PERC_GUAN," << 100.0 * nGuan / nAtoms
                << std::endl;
      std::cout << wsName << ",SITE_PERC_HBA," << 100.0 * nNHBA / nAtoms
                << std::endl;
      // Charges
      std::cout << wsName << ",SITE_POS_CHG," << posChg << std::endl;
      std::cout << wsName << ",SITE_NEG_CHG," << negChg << std::endl;
      std::cout << wsName << ",SITE_TOT_CHG," << posChg + negChg << std::endl;
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
