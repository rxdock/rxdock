//===-- CavitySearch.cxx - CavitySearch operation ---------------*- C++ -*-===//
//
// Part of the RxDock project, under the GNU LGPL version 3.
// Visit https://www.rxdock.org/ for more information.
// Copyright (c) 1998--2006 RiboTargets (subsequently Vernalis (R&D) Ltd)
// Copyright (c) 2006--2012 University of York
// Copyright (c) 2012--2014 University of Barcelona
// Copyright (c) 2019--2020 RxTx
// SPDX-License-Identifier: LGPL-3.0-only
//
//===----------------------------------------------------------------------===//
///
/// \file
/// CavitySearch operation.
///
//===----------------------------------------------------------------------===//

#include "rxdock/operation/CavitySearch.h"
#include "rxdock/BiMolWorkSpace.h"
#include "rxdock/CrdFileSink.h"
#include "rxdock/Error.h"
#include "rxdock/PRMFactory.h"
#include "rxdock/PsfFileSink.h"
#include "rxdock/SiteMapperFactory.h"

#include <fmt/format.h>
#include <fmt/ostream.h>

static const std::string _MAPPER = "mapper";

int rxdock::operation::cavitySearch(std::string strReceptorPrmFile,
                                    bool readDockingSite, bool writeDockingSite,
                                    bool writeMOEGrid, bool writeInsightII,
                                    bool writePsfCrd, bool listAtoms,
                                    double listDistance,
                                    bool printSiteDescriptors, double border) {
  try {
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
    spRecepPrmSource->SetSection(PRMFactory::_REC_SECTION);
    PRMFactory prmFactory(spRecepPrmSource);
    ModelPtr spReceptor = prmFactory.CreateReceptor();

    DockingSitePtr spDockSite;
    std::string strDockingSiteFile = wsName + "-docking-site.json";

    // Either read the docking site from the .as file
    if (readDockingSite) {
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
          spMapperFactory->CreateFromFile(spRecepPrmSource, _MAPPER);
      spMapper->Register(spWS);
      spWS->SetReceptor(spReceptor);
      fmt::print("Site mapper: {}\n", *spMapper);

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

    fmt::print("Docking site: {}\n", *spDockSite);

    if (writeDockingSite) {
      std::ofstream ostr(strDockingSiteFile.c_str());
      json dockingSite;
      dockingSite["docking-site"] = *spDockSite;
      ostr << dockingSite;
      ostr.close();
    }

    // writing active site into MOE grid
    if (writeMOEGrid) {
      fmt::print("MOE grid writing not yet implemented, sorry.");
      return EXIT_FAILURE;
    }

    // Write PSF/CRD files to keep the rDock Viewer happy (it doesn't read MOL2
    // files yet)
    if (writePsfCrd) {
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
    if (writeInsightII) {
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

    // List all receptor atoms within given distance of any cavity
    if (listAtoms) {
      RealGridPtr spGrid = spDockSite->GetGrid();
      AtomList atomList =
          spDockSite->GetAtomList(spReceptor->GetAtomList(), 0.0, listDistance);
      fmt::print("Listing {} receptor atoms within {} A of any cavity\n",
                 atomList.size(), listDistance);
      fmt::print("DISTANCE ATOM\n");
      for (AtomListConstIter iter = atomList.begin(); iter != atomList.end();
           iter++) {
        fmt::print("{:5.3f} {}\n",
                   spGrid->GetSmoothedValue((*iter)->GetCoords()), **iter);
      }
    }

    // DM 15 Jul 2002 - print out SITE descriptors
    // Use a crude measure of solvent accessibility - count #atoms within 4A of
    // each atom Use an empirical threshold to determine if atom is exposed or
    // not
    if (printSiteDescriptors) {
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
  } catch (Error &e) {
    fmt::print("{}", e.what());
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
