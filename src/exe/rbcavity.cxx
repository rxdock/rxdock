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

// Standalone executable for generating docking site .as files for rbdock

#include <algorithm>
#include <iomanip>
#ifdef _WIN32
#else
#include <unistd.h>
#endif

#include "RbtBiMolWorkSpace.h"
#include "RbtCrdFileSink.h"
#include "RbtDockingSite.h"
#include "RbtPRMFactory.h"
#include "RbtParameterFileSource.h"
#include "RbtPsfFileSink.h"
#include "RbtSiteMapperFactory.h"

const std::string EXEVERSION =
    " ($Id: //depot/dev/client3/rdock/2013.1/src/exe/rbcavity.cxx#3 $)";

void PrintUsage(void) {
  std::cout << std::endl
            << "rbcavity - calculate docking cavities" << std::endl;
  std::cout << std::endl
            << "Usage:\trbcavity -r <ReceptorPrmFile> [-W] [-R] [-d] [-v] [-l "
               "<dist>] [-s] [-b <border>] [-m]"
            << std::endl;
  std::cout << std::endl << "Options:" << std::endl;
  std::cout
      << "\t\t-r <PrmFile> - receptor param file (contains active site params)"
      << std::endl;
  std::cout
      << "\t\t-W          - write docking cavities (plus distance grid) to "
         ".as file"
      << std::endl;
  std::cout
      << "\t\t-R          - read docking cavities (plus distance grid) from "
         ".as file"
      << std::endl;
  std::cout << "\t\t-d          - dump InsightII grids for each cavity for "
               "visualisation"
            << std::endl;
  std::cout << "\t\t-v          - dump target PSF/CRD files for rDock Viewer"
            << std::endl;
  std::cout
      << "\t\t-l <dist>   - list receptor atoms with <dist> A of any cavity"
      << std::endl;
  std::cout
      << "\t\t-s          - print SITE descriptors (counts of exposed atoms)"
      << std::endl;
  std::cout << "\t\t-b <border> - set the border around the cavities for the "
               "distance grid (default=8A)"
            << std::endl;
  std::cout << "\t\t-m          - write active site into a MOE grid"
            << std::endl;
}

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
  Rbt::PrintStdHeader(std::cout, strExeName + EXEVERSION);

  // Display brief help message if no args
  if (argc < 2) {
    PrintUsage();
    return 1;
  }

  // Command line arguments and default values
  std::string strReceptorPrmFile;
  bool bReadAS(false);  // If true, read Active Site from file
  bool bWriteAS(false); // If true, write Active Site to file
  bool bDump(false);    // If true, dump cavity grids in Insight format
  bool bViewer(false);  // If true, dump PSF/CRD files for rDock Viewer
  bool bList(false);    // If true, list atoms within 'distance' of cavity
  bool bSite(
      false); // If true, print out "SITE" descriptors (counts of exposed atoms)
  bool bMOEgrid(false); // If true, create a MOE grid file for AS visualisation
  bool bBorderArg(false); // If true, border was specified in the command
                          // line
  double border(8.0);     // Border to allow around cavities for distance grid
  double dist(5.0);

  char c;
  while ((c = getopt(argc, argv, "r:b:RWdvl:ms")) != -1) {
    switch (c) {
    case 'r':
      strReceptorPrmFile = optarg;
      break;
    case 'b':
      bBorderArg = true;
      border = std::atof(optarg);
      break;
    case 'R': // also for -ras
      bReadAS = true;
      break;
    case 'W':
      bWriteAS = true;
      break;
    case 'd':
      bDump = true;
      break;
    case 'v':
      bViewer = true;
      break;
    case 'l':
      bList = true;
      dist = std::atof(optarg);
      break;
    case 'm':
      bMOEgrid = true;
      break;
    case 's':
      bSite = true;
      break;
    case '?':
      break;
    default:
      std::cout << "WARNING: unknown argument: " << c << std::endl;
      ;
      break;
    }
  }
  std::cout << std::endl;

  // check for parameter file name
  /*if(prmFile) {
          strReceptorPrmFile	= prmFile;
  }*/
  if (strReceptorPrmFile.empty()) {
    std::cout << "Missing receptor parameter file name" << std::endl;
    std::exit(1);
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
  if (bReadAS)
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

  try {
    // Create a bimolecular workspace
    RbtBiMolWorkSpacePtr spWS(new RbtBiMolWorkSpace());
    // Set the workspace name to the root of the receptor .prm filename
    std::vector<std::string> componentList =
        Rbt::ConvertDelimitedStringToList(strReceptorPrmFile, ".");
    std::string wsName = componentList.front();
    spWS->SetName(wsName);

    // Read the protocol parameter file
    RbtParameterFileSourcePtr spRecepPrmSource(new RbtParameterFileSource(
        Rbt::GetRbtFileName("data/receptors", strReceptorPrmFile)));

    // Create the receptor model from the file names in the parameter file
    spRecepPrmSource->SetSection();
    RbtPRMFactory prmFactory(spRecepPrmSource);
    RbtModelPtr spReceptor = prmFactory.CreateReceptor();

    RbtDockingSitePtr spDockSite;
    std::string strASFile = wsName + ".as";

    // Either read the docking site from the .as file
    if (bReadAS) {
      std::string strInputFile = Rbt::GetRbtFileName("data/grids", strASFile);
#if defined(__sgi) && !defined(__GNUC__)
      std::ifstream istr(strInputFile.c_str(), std::ios_base::in);
#else
      std::ifstream istr(strInputFile.c_str(),
                         std::ios_base::in | std::ios_base::binary);
#endif
      spDockSite = RbtDockingSitePtr(new RbtDockingSite(istr));
      istr.close();
    }
    // Or map the site using the prescribed mapping algorithm
    else {
      RbtSiteMapperFactoryPtr spMapperFactory(new RbtSiteMapperFactory());
      RbtSiteMapperPtr spMapper =
          spMapperFactory->CreateFromFile(spRecepPrmSource, "MAPPER");
      spMapper->Register(spWS);
      spWS->SetReceptor(spReceptor);
      std::cout << *spMapper << std::endl;

      int nRI = spReceptor->GetNumSavedCoords() - 1;
      if (nRI == 0) {
        spDockSite =
            RbtDockingSitePtr(new RbtDockingSite((*spMapper)(), border));
      } else {
        RbtCavityList allCavs;
        for (int i = 1; i <= nRI; i++) {
          spReceptor->RevertCoords(i);
          RbtCavityList cavs((*spMapper)());
          std::copy(cavs.begin(), cavs.end(), std::back_inserter(allCavs));
        }
        spDockSite = RbtDockingSitePtr(new RbtDockingSite(allCavs, border));
      }
    }

    std::cout << std::endl
              << "DOCKING SITE" << std::endl
              << (*spDockSite) << std::endl;

    if (bWriteAS) {
#if defined(__sgi) && !defined(__GNUC__)
      std::ofstream ostr(strASFile.c_str(),
                         std::ios_base::out | std::ios_base::trunc);
#else
      std::ofstream ostr(strASFile.c_str(), std::ios_base::out |
                                                std::ios_base::binary |
                                                std::ios_base::trunc);
#endif
      spDockSite->Write(ostr);
      ostr.close();
    }

    // Write PSF/CRD files to keep the rDock Viewer happy (it doesn't read MOL2
    // files yet)
    if (bViewer) {
      RbtMolecularFileSinkPtr spRecepSink =
          new RbtPsfFileSink(wsName + "_for_viewer.psf", spReceptor);
      std::cout << "Writing PSF file: " << spRecepSink->GetFileName()
                << std::endl;
      spRecepSink->Render();
      spRecepSink = new RbtCrdFileSink(wsName + "_for_viewer.crd", spReceptor);
      std::cout << "Writing CRD file: " << spRecepSink->GetFileName()
                << std::endl;
      spRecepSink->Render();
    }

    // Write an ASCII InsightII grid file for each defined cavity
    if (bDump) {
      RbtCavityList cavList = spDockSite->GetCavityList();
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
      RbtRealGridPtr spGrid = spDockSite->GetGrid();
      RbtAtomList atomList =
          spDockSite->GetAtomList(spReceptor->GetAtomList(), 0.0, dist);
      std::cout << atomList.size() << " receptor atoms within " << dist
                << " A of any cavity" << std::endl;
      std::cout << std::endl << "DISTANCE,ATOM" << std::endl;
      for (RbtAtomListConstIter iter = atomList.begin(); iter != atomList.end();
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
      // RbtRealGridPtr spGrid = spDockSite->GetGrid();
      RbtAtomList recepAtomList = spReceptor->GetAtomList();
      RbtAtomList cavAtomList = spDockSite->GetAtomList(recepAtomList, cavDist);
      RbtAtomList exposedAtomList; // The list of exposed cavity atoms
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
      for (RbtAtomListConstIter iter = cavAtomList.begin();
           iter != cavAtomList.end(); iter++) {
        int nNeighb = Rbt::GetNumAtoms(
            recepAtomList,
            Rbt::isAtomInsideSphere((*iter)->GetCoords(), neighbR));
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
      for (RbtAtomListConstIter iter = exposedAtomList.begin();
           iter != exposedAtomList.end(); iter++) {
        double chg = (*iter)->GetGroupCharge();
        if (chg > 0.0) {
          posChg += chg;
        } else if (chg < 0.0) {
          negChg += chg;
        }
      }

      // Atom type counts
      Rbt::isHybridState_eq bIsArom(RbtAtom::AROM);
      int nAtoms = exposedAtomList.size();
      int nLipoC = Rbt::GetNumAtoms(exposedAtomList, Rbt::isAtomLipophilic());
      int nArom = Rbt::GetNumAtoms(exposedAtomList, bIsArom);
      int nNHBD = Rbt::GetNumAtoms(exposedAtomList, Rbt::isAtomHBondDonor());
      int nMetal = Rbt::GetNumAtoms(exposedAtomList, Rbt::isAtomMetal());
      int nGuan =
          Rbt::GetNumAtoms(exposedAtomList, Rbt::isAtomGuanidiniumCarbon());
      int nNHBA = Rbt::GetNumAtoms(exposedAtomList, Rbt::isAtomHBondAcceptor());

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
  } catch (RbtError &e) {
    std::cout << e << std::endl;
  } catch (...) {
    std::cout << "Unknown exception" << std::endl;
  }

  _RBTOBJECTCOUNTER_DUMP_(std::cout)

  return 0;
}
