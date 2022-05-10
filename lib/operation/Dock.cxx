//===-- Dock.cxx - Dock operation -------------------------------*- C++ -*-===//
//
// Part of the RxDock project, under the GNU LGPL version 3.
// Visit https://rxdock.gitlab.io/ for more information.
// Copyright (c) 1998--2006 RiboTargets (subsequently Vernalis (R&D) Ltd)
// Copyright (c) 2006--2012 University of York
// Copyright (c) 2012--2014 University of Barcelona
// Copyright (c) 2019--2020 RxTx
// SPDX-License-Identifier: LGPL-3.0-only
//
//===----------------------------------------------------------------------===//
///
/// \file
/// Dock operation.
///
//===----------------------------------------------------------------------===//

#include "rxdock/operation/Dock.h"
#include "rxdock/BiMolWorkSpace.h"
#include "rxdock/CrdFileSink.h"
#include "rxdock/DockingError.h"
#include "rxdock/Error.h"
#include "rxdock/LigandError.h"
#include "rxdock/MdlFileSink.h"
#include "rxdock/MdlFileSource.h"
#include "rxdock/PRMFactory.h"
#include "rxdock/ParameterFileSource.h"
#include "rxdock/SFFactory.h"
#include "rxdock/TransformFactory.h"

#include <fmt/chrono.h>
#include <fmt/format.h>
#include <fmt/ostream.h>

#include <chrono>

namespace rxdock {
namespace operation {

static const std::string _ROOT_SF = "rxdock.score";
static const std::string _RESTRAINT_SF = "restr";
static const std::string _ROOT_TRANSFORM = "dock";

} // namespace operation
} // namespace rxdock

int rxdock::operation::dock(
    std::string strLigandMdlFile, std::string strOutputMdlFile, bool bOutputCrd,
    std::string strOutputCrdFile, bool bOutputHistory,
    std::string strOutputHistoryFilePrefix, std::string strReceptorPrmFile,
    std::string strParamFile, bool bFilter, std::string strFilterFile,
    bool bDockingRuns, std::size_t nDockingRuns, bool bPosIonise,
    bool bNegIonise, bool bExplH, bool bTarget, double dTargetScore,
    bool bContinue, bool bSeed, std::size_t nSeed) {
  try {
    // Create a bimolecular workspace
    BiMolWorkSpacePtr spWS(new BiMolWorkSpace());
    // Set the workspace name to the root of the receptor .prm filename
    std::vector<std::string> componentList =
        ConvertDelimitedStringToList(strReceptorPrmFile, ".");
    std::string wsName = componentList.front();
    spWS->SetName(wsName);

    // Read the docking protocol parameter file
    ParameterFileSourcePtr spParamSource(
        new ParameterFileSource(GetDataFileName("data/scripts", strParamFile)));
    // Read the receptor parameter file
    ParameterFileSourcePtr spRecepPrmSource(new ParameterFileSource(
        GetDataFileName("data/receptors", strReceptorPrmFile)));
    fmt::print("Docking protocol: {} {}\n", spParamSource->GetFileName(),
               spParamSource->GetTitle());
    fmt::print("Receptor: {} {}\n", spRecepPrmSource->GetFileName(),
               spRecepPrmSource->GetTitle());

    // Create the scoring function from the rxdock.score section of the docking
    // protocol prm file Format is: SECTION rxdock.score
    //    inter    InterSF.prm
    //    intra IntraSF.prm
    // END_SECTION
    //
    // Notes:
    // Section name must be rxdock.score. This is also the name of the root SF
    // aggregate An aggregate is created for each parameter in the section.
    // Parameter name becomes the name of the subaggregate (e.g.
    // rxdock.score.inter) Parameter value is the file name for the subaggregate
    // definition Default directory is $RBT_ROOT/data/sf
    SFFactoryPtr spSFFactory(
        new SFFactory());               // Factory class for scoring functions
    SFAggPtr spSF(new SFAgg(_ROOT_SF)); // Root SF aggregate
    spParamSource->SetSection(_ROOT_SF);
    std::vector<std::string> sfList(spParamSource->GetParameterList());
    // Loop over all parameters in the rxdock.score section
    for (std::vector<std::string>::const_iterator sfIter = sfList.begin();
         sfIter != sfList.end(); sfIter++) {
      // sfFile = file name for scoring function subaggregate
      std::string sfFile(GetDataFileName(
          "data/sf", spParamSource->GetParameterValueAsString(*sfIter)));
      ParameterFileSourcePtr spSFSource(new ParameterFileSource(sfFile));
      // Create and add the subaggregate
      spSF->Add(spSFFactory->CreateAggFromFile(spSFSource, *sfIter));
    }

    // Add the RESTRAINT subaggregate scoring function from any SF definitions
    // in the receptor prm file
    spSF->Add(spSFFactory->CreateAggFromFile(spRecepPrmSource, _RESTRAINT_SF));

    // Create the docking transform aggregate from the transform definitions in
    // the docking prm file
    TransformFactoryPtr spTransformFactory(new TransformFactory());
    spParamSource->SetSection();
    TransformAggPtr spTransform(
        spTransformFactory->CreateAggFromFile(spParamSource, _ROOT_TRANSFORM));

    // Print the scoring function and the transform details
    fmt::print("Scoring function details: {}\n", *spSF);
    fmt::print("Search details: {}\n", *spTransform);

    // Register the scoring function and the transform with the workspace
    spWS->SetSF(spSF);
    spWS->SetTransform(spTransform);

    // DM 18 May 1999
    // Variants describing the library version, parameter file, and current
    // directory will be stored in the ligand SD files
    Variant vLib(GetProduct() + "/" + GetProgramVersion());
    Variant vRecep(spRecepPrmSource->GetFileName());
    Variant vPrm(spParamSource->GetFileName());
    Variant vDir(GetCurrentWorkingDirectory());

    spRecepPrmSource->SetSection();
    // Read docking site from file and register with workspace
    std::string strDockingSiteFile = spWS->GetName() + "-docking-site.json";
    std::string strInputFile =
        GetDataFileName("data/grids", strDockingSiteFile);
    // DM 14 June 2006 - bug fix to one of the longest standing rDock issues
    //(the cryptic "Error reading from input stream" message, if cavity file was
    // missing)
    // std::string message = "Cavity file (" + strDockingSiteFile +
    //                       ") not found in current directory or $RBT_HOME";
    // message += " - run rbcavity first";
    // throw FileReadError(_WHERE_, message);
    std::ifstream inputFile(strInputFile.c_str());
    json siteData;
    inputFile >> siteData;
    inputFile.close();
    DockingSitePtr spDS(new DockingSite(siteData.at("docking-site")));
    spWS->SetDockingSite(spDS);
    fmt::print("Docking site: {}\n", *spDS);

    // Prepare the SD file sink for saving the docked conformations for each
    // ligand DM 3 Dec 1999 - replaced ostrstream with String in determining
    // SD file name SRC 2014 moved here this block to allow WRITE_ERROR TRUE
    MolecularFileSinkPtr spMdlFileSink(
        new MdlFileSink(strOutputMdlFile, ModelPtr()));
    spWS->SetSink(spMdlFileSink);

    PRMFactory prmFactory(spRecepPrmSource, spDS);
    // Create the receptor model from the file names in the receptor parameter
    // file
    ModelPtr spReceptor = prmFactory.CreateReceptor();
    spWS->SetReceptor(spReceptor);

    // Register any solvent
    ModelList solventList = prmFactory.CreateSolvent();
    spWS->SetSolvent(solventList);
    if (spWS->hasSolvent()) {
      int nSolvent = spWS->GetSolvent().size();
      fmt::print("{} solvent molecules registered\n", nSolvent);
    } else {
      fmt::print("No solvent registered\n");
    }

    // Seed the random number generator
    Rand &theRand = GetRandInstance(); // ref to random number generator
    if (bSeed) {
      theRand.Seed(nSeed);
    }

    // BGD 26 Feb 2003 - Create filters to simulate old rbdock
    // behaviour
    std::ostringstream strFilter;
    if (!bFilter) {
      if (bTarget) // -t<TS>
      {
        fmt::print("Lower target intermolecular score = {}\n", dTargetScore);
        if (!bDockingRuns) // -t<TS> only
        {
          strFilter << fmt::format("0 1 - {}score.inter ", GetMetaDataPrefix())
                    << dTargetScore << std::endl;
        } else             // -t<TS> -n<N> need to check if -cont present
                           // for all other cases it doesn't matter
            if (bContinue) // -t<TS> -n<N> -cont
        {
          strFilter << fmt::format("1 if - {}score.NRUNS ", GetMetaDataPrefix())
                    << (nDockingRuns - 1)
                    << fmt::format(" 0.0 -1.0,\n1 - {}score.inter ",
                                   GetMetaDataPrefix())
                    << dTargetScore << std::endl;
        } else // -t<TS> -n<N>
        {
          strFilter << "1 if - " << dTargetScore
                    << fmt::format(" {}score.inter 0.0 ", GetMetaDataPrefix())
                    << fmt::format("if - {}score.NRUNS ", GetMetaDataPrefix())
                    << (nDockingRuns - 1)
                    << fmt::format(" 0.0 -1.0,\n1 - {}score.inter ",
                                   GetMetaDataPrefix())
                    << dTargetScore << std::endl;
        }
      }                      // no target score, no filter
      else if (bDockingRuns) // -n<N>
      {
        strFilter << fmt::format("1 if - {}score.NRUNS ", GetMetaDataPrefix())
                  << (nDockingRuns - 1) << " 0.0 -1.0,\n0";
      } else // no -t no -n
      {
        strFilter << "0 0\n";
      }
    }

    // Create the filter object for controlling early termination of protocol
    FilterPtr spfilter;
    if (bFilter) {
      spfilter = new Filter(strFilterFile);
      if (bDockingRuns) {
        spfilter->SetMaxNRuns(nDockingRuns);
      }
    } else {
      spfilter = new Filter(strFilter.str(), true);
    }

    // Register the Filter with the workspace
    spWS->SetFilter(spfilter);

    // MAIN LOOP OVER LIGAND RECORDS
    // DM 20 Apr 1999 - set the auto-ionise flags
    if (bPosIonise) {
      fmt::print("Automatically protonating positive ionisable groups (amines, "
                 "imidazoles, and guanidines)\n");
    }
    if (bNegIonise) {
      fmt::print(
          "Automatically deprotonating negative ionisable groups (carboxylic "
          "acids, phosphates, sulphates, and sulphonates)\n");
    }
    if (!bExplH) {
      fmt::print("Reading polar hydrogens only from ligand SD file\n");
    } else {
      fmt::print("Reading all hydrogens from ligand SD file\n");
    }
    // DM 20 Apr 1999 - add explicit bPosIonise and bNegIonise flags to
    // MdlFileSource constructor
    MolecularFileSourcePtr spMdlFileSource(
        new MdlFileSource(strLigandMdlFile, bPosIonise, bNegIonise, !bExplH));
    std::chrono::duration<double> totalDuration(0.0);
    std::size_t nFailedLigands = 0;
    std::size_t nUnnamedLigands = 0;
    std::chrono::system_clock::time_point loopBegin =
        std::chrono::system_clock::now();
    std::size_t nRec;
    for (nRec = 0; spMdlFileSource->FileStatusOK();
         spMdlFileSource->NextRecord(), nRec++) {
      fmt::print("SDfile record #{}\n", nRec + 1);
      Error molStatus = spMdlFileSource->Status();
      if (!molStatus.isOK()) {
        fmt::print("{}\n", molStatus.what());
        continue;
      }

      auto startTime = std::chrono::high_resolution_clock::now();
      bool bLigandError = false;

      // DM 26 Jul 1999 - only read the largest segment (guaranteed to be called
      // H) BGD 07 Oct 2002 - catching errors created by the ligands, so rbdock
      // continues with the next one, instead of completely stopping
      try {
        spMdlFileSource->SetSegmentFilterMap(ConvertStringToSegmentMap("H"));

        if (spMdlFileSource->isDataFieldPresent("Name")) {
          Variant molName = spMdlFileSource->GetDataValue("Name");
          if (molName.isEmpty())
            nUnnamedLigands++;
          fmt::print("Name: {}\n", molName.GetString());
        }
        if (spMdlFileSource->isDataFieldPresent("REG_Number"))
          fmt::print("REG_Number: {}\n",
                     spMdlFileSource->GetDataValue("REG_Number").GetString());
        fmt::print("RNG seed: ");
        if (bSeed) {
          fmt::print("{}\n", nSeed);
        } else {
          fmt::print("std::random_device\n");
        }

        // Create and register the ligand model
        ModelPtr spLigand = prmFactory.CreateLigand(spMdlFileSource);
        std::string strMolName = spLigand->GetName();
        spWS->SetLigand(spLigand);
        // Update any model coords from embedded chromosomes in the ligand file
        spWS->UpdateModelCoordsFromChromRecords(spMdlFileSource);

        // DM 18 May 1999 - store run info in model data
        // Clear any previous rxdock.program.* data fields
        spLigand->ClearAllDataFields(GetMetaDataPrefix() + "program.");
        spLigand->SetDataValue(GetMetaDataPrefix() + "program.library", vLib);
        spLigand->SetDataValue(GetMetaDataPrefix() + "program.receptor",
                               vRecep);
        spLigand->SetDataValue(GetMetaDataPrefix() + "program.parameter_file",
                               vPrm);
        spLigand->SetDataValue(
            GetMetaDataPrefix() + "program.current_directory", vDir);

        // DM 10 Dec 1999 - if in target mode, loop until target score is
        // reached
        bool bTargetMet = false;

        ////////////////////////////////////////////////////
        // MAIN LOOP OVER EACH SIMULATED ANNEALING RUN
        // Create a history file sink, just in case it's needed by any
        // of the transforms
        std::size_t iRun = 0;
        std::size_t nErrors = 0;
        // need to check this here. The termination
        // filter is only run once at least
        // one docking run has been done.
        if (nDockingRuns < 1)
          bTargetMet = true;
        while (!bTargetMet) {
          // Catching errors with this specific run
          if (nErrors > 10) {
            fmt::print(
                "Target not met, but giving up on ligand after {} errors",
                nErrors);
            bLigandError = true;
            break;
          }
          try {
            if (bOutputHistory) {
              std::ostringstream histr;
              histr << strOutputHistoryFilePrefix << "_" << strMolName
                    << nRec + 1 << "_his_" << iRun + 1 << ".sd";
              MolecularFileSinkPtr spHistoryFileSink(
                  new MdlFileSink(histr.str(), spLigand));
              spWS->SetHistorySink(spHistoryFileSink);
            }
            spWS->Run(); // Dock!
            bool bterm = spfilter->Terminate();
            bool bwrite = spfilter->Write();
            if (bterm)
              bTargetMet = true;
            if (bwrite) {
              spWS->Save();
            }
            iRun++;
          } catch (DockingError &e) {
            fmt::print("{}\n", e.what());
            nErrors++;
          }
        }
        // END OF MAIN LOOP OVER EACH SIMULATED ANNEALING RUN
        ////////////////////////////////////////////////////

        // here we use iRun since it got incremented in the last iteration to
        // the number of runs done
        fmt::print("Numer of docking runs done:   {} ({} errors)\n", iRun,
                   nErrors);
      }
      // END OF TRY
      catch (LigandError &e) {
        fmt::print("{}\n", e.what());
        bLigandError = true;
      }

      if (!bLigandError) {
        auto endTime = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> recordDuration = endTime - startTime;
        fmt::print("Ligand docking duration:      {} second(s)\n",
                   recordDuration.count());
        totalDuration += recordDuration;
        // report average every 10th record starting from the 1st
        // nRec ligand is done here so the number of docked ligands is nRec + 1
        if (nRec % 10 == 0) {
          fmt::print("\nAverage duration per ligand:  {} second(s)\n",
                     totalDuration.count() / static_cast<double>(nRec + 1));
          std::size_t estNumRecords = spMdlFileSource->GetEstimatedNumRecords();
          if (estNumRecords > 0) {
            std::chrono::duration<double> estimatedTimeRemaining =
                estNumRecords * (totalDuration / static_cast<double>(nRec + 1));
            std::chrono::system_clock::time_point loopEnd =
                loopBegin + std::chrono::duration_cast<std::chrono::seconds>(
                                estimatedTimeRemaining);
            std::time_t loopEndTime =
                std::chrono::system_clock::to_time_t(loopEnd);
            fmt::print(
                "Approximately {} record(s) remaining, will finish {:%c}\n",
                estNumRecords - (nRec + 1), fmt::localtime(loopEndTime));
          }
        }
      } else {
        nFailedLigands++;
      }
    }
    // END OF MAIN LOOP OVER LIGAND RECORDS
    ////////////////////////////////////////////////////

    // FileStatusOK becomes false when nRec, counting from zero, becomes equal
    // to number of ligands in the file
    fmt::print("Total number of ligands: {}", nRec);
    if (nFailedLigands > 0)
      fmt::print(", of which {} failed to dock\n", nFailedLigands);
    else
      fmt::print(", all ligands docked without errors\n");

    if (nRec - nFailedLigands > 0) {
      auto hTotal =
          std::chrono::duration_cast<std::chrono::hours>(totalDuration);
      totalDuration -= hTotal;
      auto mTotal =
          std::chrono::duration_cast<std::chrono::minutes>(totalDuration);
      totalDuration -= mTotal;

      fmt::print("\nDocking duration for {} ligand(s): ",
                 nRec - nFailedLigands);

      if (hTotal.count() > 0) {
        fmt::print("{} hour(s), ", hTotal.count());
      }
      if (hTotal.count() > 0 || mTotal.count() > 0) {
        fmt::print("{} minute(s), ", mTotal.count());
      }
      fmt::print("{} second(s)\n", totalDuration.count());
    }

    if (nUnnamedLigands > 0) {
      fmt::print("Warning: {} ligand(s) are unnamed. Post-processing tools "
                 "might have an issue correctly identifying different poses of "
                 "the same ligand.\n",
                 nUnnamedLigands);
    }

    if (bOutputCrd) {
      MolecularFileSinkPtr spRecepSink(
          new CrdFileSink(strOutputCrdFile, spReceptor));
      spRecepSink->Render();
    }
    std::cout << std::endl;
    PrintBibliographyItem(std::cout, "RiboDock2004");
    PrintBibliographyItem(std::cout, "rDock2014");
#if !defined(__sun) && !defined(_MSC_VER)
    PrintBibliographyItem(std::cout, "PCG2014");
#endif
    std::cout << std::endl;
    std::cout << "Thank you for using " << GetProgramName() << " "
              << GetProgramVersion() << "." << std::endl;
  } catch (Error &e) {
    fmt::print("{}\n", e.what());
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
