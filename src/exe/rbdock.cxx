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

// Main docking application
#include <cerrno>
#include <chrono>
#include <cxxopts.hpp>
#include <iomanip>

#include "BiMolWorkSpace.h"
#include "CrdFileSink.h"
#include "DockingError.h"
#include "FileError.h"
#include "Filter.h"
#include "LigandError.h"
#include "MdlFileSink.h"
#include "MdlFileSource.h"
#include "ModelError.h"
#include "PRMFactory.h"
#include "ParameterFileSource.h"
#include "Rand.h"
#include "SFFactory.h"
#include "SFRequest.h"
#include "TransformFactory.h"

using namespace rxdock;

namespace rxdock {

// Section name in docking prm file containing scoring function definition
const std::string _ROOT_SF = "SCORE";
const std::string _RESTRAINT_SF = "RESTR";
const std::string _ROOT_TRANSFORM = "DOCK";

} // namespace rxdock

/////////////////////////////////////////////////////////////////////
// MAIN PROGRAM STARTS HERE
/////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[]) {
  // Handle obsolete arguments, if any
  int tIndex = 0;
  bool CFound = false;
  for (int i = 0; i < argc; i++) {
    std::string opt = argv[i];
    if (opt == "-ap" || opt == "-an" || opt == "-allH" || opt == "-cont") {
      std::cout
          << "Options -ap, -an, -allH, and -cont are no longer supported; use "
             "-P, -D, -H, and -C (respectively) instead."
          << std::endl;
      return 2;
    }
    if (opt == "-t") {
      tIndex = i;
    }
    if (opt == "-C") {
      CFound = true;
    }
  }
  if (tIndex > 0 && !CFound) {
    std::cout << "Option -t without -C is no longer supported; use -f instead."
              << std::endl;
    return 2;
  }

  std::cout.setf(std::ios_base::left, std::ios_base::adjustfield);

  // Strip off the path to the executable, leaving just the file name
  std::string strExeName(argv[0]);
  std::string::size_type i = strExeName.rfind("/");
  if (i != std::string::npos)
    strExeName.erase(0, i + 1);

  // Print a standard header
  PrintStdHeader(std::cout, strExeName);

  cxxopts::Options options(strExeName, "rbdock - docking engine");

  // Command line arguments and default values
  cxxopts::OptionAdder adder = options.add_options();
  adder("i,input", "input ligand SD file", cxxopts::value<std::string>());
  adder("o,output", "output file name(s) prefix",
        cxxopts::value<std::string>());
  adder("r,receptor-param", "receptor parameter file",
        cxxopts::value<std::string>());
  adder("p,docking-param", "docking protocol parameter file",
        cxxopts::value<std::string>());
  adder("n,number", "number of runs per ligand (0 = unlimited)",
        cxxopts::value<unsigned int>()->default_value("0"));
  adder("P,protonate", "protonate all neutral amines, guanidines, imidazoles");
  adder("D,deprotonate",
        "deprotonate all carboxylic, sulphur and phosphorous acid groups");
  adder("H,all-hydrogens",
        "read all hydrogens present instead of only polar hydrogens");
  adder("t,threshold", "score threshold", cxxopts::value<double>());
  adder("C,continue",
        "continue if score threshold is met instead of terminating ligand");
  adder("f,filter", "filter file name", cxxopts::value<std::string>());
  adder("T,trace",
        "controls output level for debugging (0 = minimal, >0 = more verbose)",
        cxxopts::value<int>()->default_value("0"));
  adder("s,seed", "random number seed to use instead of std::random_device",
        cxxopts::value<unsigned int>());
  adder("h,help", "Print help");

  try {
    auto result = options.parse(argc, argv);

    if (result.count("h")) {
      std::cout << options.help() << std::endl;
      return 0;
    }

    // Command line arguments and default values
    std::string strLigandMdlFile;
    if (result.count("i")) {
      strLigandMdlFile = result["i"].as<std::string>();
    }
    bool bOutput = result.count("o");
    std::string strRunName;
    if (bOutput) {
      strRunName = result["o"].as<std::string>();
    }
    std::string strReceptorPrmFile;
    if (result.count("r")) {
      strReceptorPrmFile = result["r"].as<std::string>(); // Receptor param file
    }
    std::string strParamFile;
    if (result.count("p")) {
      strParamFile = result["p"].as<std::string>(); // Docking run param file
    }
    bool bFilter = result.count("f");
    std::string strFilterFile;
    if (bFilter) {
      strFilterFile = result["f"].as<std::string>(); // Filter file
    }
    bool bDockingRuns = result.count("n"); // is argument -n present?
    unsigned int nDockingRuns =
        result["n"].as<unsigned int>(); // Defaults to zero, so can detect later
                                        // whether user explictly typed -n

    bool bPosIonise = result.count("P");
    bool bNegIonise = result.count("D");
    bool bExplH = result.count("H"); // if true, read only polar hydrogens from
                                     // SD file, else read all H's present

    // Params for target score
    bool bTarget = result.count("t");
    double dTargetScore;
    if (bTarget) {
      dTargetScore = result["t"].as<double>();
    }
    bool bContinue =
        result.count("C"); // DM 25 May 2001 - if true, stop once target is met

    bool bSeed = result.count(
        "s"); // Random number seed (default = from std::random_device)
    unsigned int nSeed;
    if (bSeed) {
      nSeed = result["s"].as<unsigned int>();
    }
    bool bTrace = result.count("T");
    int iTrace = result["T"].as<int>(); // Trace level, for debugging

    // input ligand file, receptor and parameter is compulsory
    if (strLigandMdlFile.empty() || strReceptorPrmFile.empty() ||
        strParamFile.empty()) { // if any of them is missing
      std::cout << "Missing required parameter(s): -i, -r, or -p" << std::endl;
      return 1;
    }

    // print out arguments
    std::cout << std::endl << "Command line args:" << std::endl;
    std::cout << " -i " << strLigandMdlFile << std::endl;
    std::cout << " -r " << strReceptorPrmFile << std::endl;
    std::cout << " -p " << strParamFile << std::endl;
    // output is not that important but good to have
    if (!strRunName.empty()) {
      bOutput = true;
      std::cout << " -o " << strRunName << std::endl;
    } else {
      std::cout << "WARNING: output file name is missing." << std::endl;
    }
    // docking runs
    if (nDockingRuns >= 1) { // User typed -n explicitly
      std::cout << " -n " << nDockingRuns << std::endl;
    } else {
      nDockingRuns = 1; // User didn't type -n explicitly, so fall back to the
                        // default of n=1
      std::cout << " -n " << nDockingRuns << " (default) " << std::endl;
    }
    if (bSeed) // random seed (if provided)
      std::cout << " -s " << nSeed << std::endl;
    if (bTrace) // random seed (if provided)
      std::cout << " -T " << iTrace << std::endl;
    if (bPosIonise) // protonate
      std::cout << " -ap " << std::endl;
    if (bNegIonise) // deprotonate
      std::cout << " -an " << std::endl;
    if (bExplH) // all-H
      std::cout << " -allH " << std::endl;
    if (bContinue) // stop after target
      std::cout << " -cont " << std::endl;
    if (bTarget)
      std::cout << " -t " << dTargetScore << std::endl;

    // BGD 26 Feb 2003 - Create filters to simulate old rbdock
    // behaviour
    std::ostringstream strFilter;
    if (!bFilter) {
      if (bTarget) // -t<TS>
      {
        if (!bDockingRuns) // -t<TS> only
        {
          strFilter << "0 1 - SCORE.INTER " << dTargetScore << std::endl;
        } else             // -t<TS> -n<N> need to check if -cont present
                           // for all other cases it doesn't matter
            if (bContinue) // -t<TS> -n<N> -cont
        {
          strFilter << "1 if - SCORE.NRUNS " << (nDockingRuns - 1)
                    << " 0.0 -1.0,\n1 - SCORE.INTER " << dTargetScore
                    << std::endl;
        } else // -t<TS> -n<N>
        {
          strFilter << "1 if - " << dTargetScore << " SCORE.INTER 0.0 "
                    << "if - SCORE.NRUNS " << (nDockingRuns - 1)
                    << " 0.0 -1.0,\n1 - SCORE.INTER " << dTargetScore
                    << std::endl;
        }
      }                      // no target score, no filter
      else if (bDockingRuns) // -n<N>
      {
        strFilter << "1 if - SCORE.NRUNS " << (nDockingRuns - 1)
                  << " 0.0 -1.0,\n0";
      } else // no -t no -n
      {
        strFilter << "0 0\n";
      }
    }

    // DM 20 Apr 1999 - set the auto-ionise flags
    if (bPosIonise)
      std::cout
          << "Automatically protonating positive ionisable groups (amines, "
             "imidazoles, guanidines)"
          << std::endl;
    if (bNegIonise)
      std::cout << "Automatically deprotonating negative ionisable groups "
                   "(carboxylic "
                   "acids, phosphates, sulphates, sulphonates)"
                << std::endl;
    if (!bExplH)
      std::cout << "Reading polar hydrogens only from ligand SD file"
                << std::endl;
    else
      std::cout << "Reading all hydrogens from ligand SD file" << std::endl;

    if (bTarget) {
      std::cout << std::endl
                << "Lower target intermolecular score = " << dTargetScore
                << std::endl;
    }

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
    std::cout << std::endl
              << "DOCKING PROTOCOL:" << std::endl
              << spParamSource->GetFileName() << std::endl
              << spParamSource->GetTitle() << std::endl;
    std::cout << std::endl
              << "RECEPTOR:" << std::endl
              << spRecepPrmSource->GetFileName() << std::endl
              << spRecepPrmSource->GetTitle() << std::endl;

    // Create the scoring function from the SCORE section of the docking
    // protocol prm file Format is: SECTION SCORE
    //    INTER    InterSF.prm
    //    INTRA IntraSF.prm
    // END_SECTION
    //
    // Notes:
    // Section name must be SCORE. This is also the name of the root SF
    // aggregate An aggregate is created for each parameter in the section.
    // Parameter name becomes the name of the subaggregate (e.g. SCORE.INTER)
    // Parameter value is the file name for the subaggregate definition
    // Default directory is $RBT_ROOT/data/sf
    SFFactoryPtr spSFFactory(
        new SFFactory());               // Factory class for scoring functions
    SFAggPtr spSF(new SFAgg(_ROOT_SF)); // Root SF aggregate
    spParamSource->SetSection(_ROOT_SF);
    std::vector<std::string> sfList(spParamSource->GetParameterList());
    // Loop over all parameters in the SCORE section
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

    // Override the TRACE levels for the scoring function and transform
    // Dump details to std::cout
    // Register the scoring function and the transform with the workspace
    if (bTrace) {
      RequestPtr spTraceReq(new SFSetParamRequest("TRACE", iTrace));
      spSF->HandleRequest(spTraceReq);
      spTransform->HandleRequest(spTraceReq);
    }
    if (iTrace > 0) {
      std::cout << std::endl
                << "SCORING FUNCTION DETAILS:" << std::endl
                << *spSF << std::endl;
      std::cout << std::endl
                << "SEARCH DETAILS:" << std::endl
                << *spTransform << std::endl;
    }
    spWS->SetSF(spSF);
    spWS->SetTransform(spTransform);

    // DM 18 May 1999
    // Variants describing the library version, exe version, parameter file, and
    // current directory Will be stored in the ligand SD files
    Variant vLib(GetProduct() + "/" + GetProgramVersion());
    Variant vExe(strExeName + "/" + GetProgramVersion());
    Variant vRecep(spRecepPrmSource->GetFileName());
    Variant vPrm(spParamSource->GetFileName());
    Variant vDir(GetCurrentWorkingDirectory());

    spRecepPrmSource->SetSection();
    // Read docking site from file and register with workspace
    std::string strASFile = spWS->GetName() + ".as";
    std::string strInputFile = GetDataFileName("data/grids", strASFile);
    // DM 26 Sep 2000 - std::ios_base::binary is invalid with IRIX CC
#if defined(__sgi) && !defined(__GNUC__)
    std::ifstream istr(strInputFile.c_str(), std::ios_base::in);
#else
    std::ifstream istr(strInputFile.c_str(),
                       std::ios_base::in | std::ios_base::binary);
#endif
    // DM 14 June 2006 - bug fix to one of the longest standing rDock issues
    //(the cryptic "Error reading from input stream" message, if cavity file was
    // missing)
    if (!istr) {
      std::string message = "Cavity file (" + strASFile +
                            ") not found in current directory or $RBT_HOME";
      message += " - run rbcavity first";
      throw FileReadError(_WHERE_, message);
    }
    DockingSitePtr spDS(new DockingSite(istr));
    istr.close();
    spWS->SetDockingSite(spDS);
    std::cout << std::endl
              << "DOCKING SITE" << std::endl
              << (*spDS) << std::endl;

    // Prepare the SD file sink for saving the docked conformations for each
    // ligand DM 3 Dec 1999 - replaced ostrstream with String in determining
    // SD file name SRC 2014 moved here this block to allow WRITE_ERROR TRUE
    if (bOutput) {
      MolecularFileSinkPtr spMdlFileSink(
          new MdlFileSink(strRunName + ".sd", ModelPtr()));
      spWS->SetSink(spMdlFileSink);
    }

    PRMFactory prmFactory(spRecepPrmSource, spDS);
    prmFactory.SetTrace(iTrace);
    // Create the receptor model from the file names in the receptor parameter
    // file
    ModelPtr spReceptor = prmFactory.CreateReceptor();
    spWS->SetReceptor(spReceptor);

    // Register any solvent
    ModelList solventList = prmFactory.CreateSolvent();
    spWS->SetSolvent(solventList);
    if (spWS->hasSolvent()) {
      int nSolvent = spWS->GetSolvent().size();
      std::cout << std::endl
                << nSolvent << " solvent molecules registered" << std::endl;
    } else {
      std::cout << std::endl << "No solvent" << std::endl;
    }

    // SRC 2014 removed sector bOutput from here to some blocks above, for
    // WRITEERRORS TRUE

    // Seed the random number generator
    Rand &theRand = GetRandInstance(); // ref to random number generator
    if (bSeed) {
      theRand.Seed(nSeed);
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
    if (bTrace) {
      RequestPtr spTraceReq(new SFSetParamRequest("TRACE", iTrace));
      spfilter->HandleRequest(spTraceReq);
    }

    // Register the Filter with the workspace
    spWS->SetFilter(spfilter);

    // MAIN LOOP OVER LIGAND RECORDS
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
    for (nRec = 1; spMdlFileSource->FileStatusOK();
         spMdlFileSource->NextRecord(), nRec++) {
      std::cout.setf(std::ios_base::left, std::ios_base::adjustfield);
      std::cout << std::endl
                << "**************************************************"
                << std::endl
                << "RECORD #" << nRec << std::endl;
      Error molStatus = spMdlFileSource->Status();
      if (!molStatus.isOK()) {
        std::cout << std::endl
                  << molStatus.what() << std::endl
                  << "************************************************"
                  << std::endl;
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
          std::cout << "NAME:   " << molName << std::endl;
        }
        if (spMdlFileSource->isDataFieldPresent("REG_Number"))
          std::cout << "REG_Num:" << spMdlFileSource->GetDataValue("REG_Number")
                    << std::endl;
        std::cout << std::setw(30) << "RNG seed:";
        if (bSeed) {
          std::cout << nSeed;
        } else {
          std::cout << "std::random_device";
        }
        std::cout << std::endl;

        // Create and register the ligand model
        ModelPtr spLigand = prmFactory.CreateLigand(spMdlFileSource);
        std::string strMolName = spLigand->GetName();
        spWS->SetLigand(spLigand);
        // Update any model coords from embedded chromosomes in the ligand file
        spWS->UpdateModelCoordsFromChromRecords(spMdlFileSource, iTrace);

        // DM 18 May 1999 - store run info in model data
        // Clear any previous Rbt.* data fields
        spLigand->ClearAllDataFields("Rbt.");
        spLigand->SetDataValue("Rbt.Library", vLib);
        spLigand->SetDataValue("Rbt.Executable", vExe);
        spLigand->SetDataValue("Rbt.Receptor", vRecep);
        spLigand->SetDataValue("Rbt.Parameter_File", vPrm);
        spLigand->SetDataValue("Rbt.Current_Directory", vDir);

        // DM 10 Dec 1999 - if in target mode, loop until target score is
        // reached
        bool bTargetMet = false;

        ////////////////////////////////////////////////////
        // MAIN LOOP OVER EACH SIMULATED ANNEALING RUN
        // Create a history file sink, just in case it's needed by any
        // of the transforms
        int iRun = 1;
        std::size_t nErrors = 0;
        // need to check this here. The termination
        // filter is only run once at least
        // one docking run has been done.
        if (nDockingRuns < 1)
          bTargetMet = true;
        while (!bTargetMet) {
          // Catching errors with this specific run
          if (nErrors > 10) {
            std::cout << "Target not met, but giving up on ligand after "
                      << nErrors << " errors" << std::endl;
            bLigandError = true;
            break;
          }
          try {
            if (bOutput) {
              std::ostringstream histr;
              histr << strRunName << "_" << strMolName << nRec << "_his_"
                    << iRun << ".sd";
              MolecularFileSinkPtr spHistoryFileSink(
                  new MdlFileSink(histr.str(), spLigand));
              spWS->SetHistorySink(spHistoryFileSink);
            }
            spWS->Run(); // Dock!
            bool bterm = spfilter->Terminate();
            bool bwrite = spfilter->Write();
            if (bterm)
              bTargetMet = true;
            if (bOutput && bwrite) {
              spWS->Save();
            }
            iRun++;
          } catch (DockingError &e) {
            std::cout << e.what() << std::endl;
            nErrors++;
          }
        }
        // END OF MAIN LOOP OVER EACH SIMULATED ANNEALING RUN
        ////////////////////////////////////////////////////

        // here we use iRun - 1 since iRun got incremented in the last iteration
        std::cout << "Numer of docking runs done:   " << iRun - 1 << " ("
                  << nErrors << " errors)" << std::endl;
      }
      // END OF TRY
      catch (LigandError &e) {
        std::cout << e.what() << std::endl;
        bLigandError = true;
      }

      if (!bLigandError) {
        auto endTime = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> recordDuration = endTime - startTime;
        std::cout << "Ligand docking duration:      " << recordDuration.count()
                  << " second(s)" << std::endl;
        totalDuration += recordDuration;
        // report average every 10th record starting from the 1st
        if (nRec % 10 == 1) {
          std::cout << std::endl
                    << "Average duration per ligand:  "
                    << totalDuration.count() / static_cast<double>(nRec)
                    << " second(s)" << std::endl;
          std::size_t estNumRecords = spMdlFileSource->GetEstimatedNumRecords();
          if (estNumRecords > 0) {
            std::chrono::duration<double> estimatedTimeRemaining =
                estNumRecords * (totalDuration / static_cast<double>(nRec));
            std::chrono::system_clock::time_point loopEnd =
                loopBegin + std::chrono::duration_cast<std::chrono::seconds>(
                                estimatedTimeRemaining);
            std::time_t loopEndTime =
                std::chrono::system_clock::to_time_t(loopEnd);
            std::cout << "Approximately " << estNumRecords - nRec
                      << " record(s) remaining, will finish "
                      << std::put_time(std::localtime(&loopEndTime), "%c")
                      << std::endl;
          }
        }
      } else {
        nFailedLigands++;
      }
    }
    // END OF MAIN LOOP OVER LIGAND RECORDS
    ////////////////////////////////////////////////////

    std::cout << std::endl
              << "**************************************************"
              << std::endl
              << "Total number of ligands: " << nRec - 1;
    if (nFailedLigands > 0)
      std::cout << ", of which" << nFailedLigands << " failed to dock"
                << std::endl;
    else
      std::cout << ", all ligands docked without errors" << std::endl;

    // here we use nRec - 1 for the number of ligands since nRec got incremented
    // in the iteration in which spMdlFileSource->FileStatusOK() returned false
    // and for loop ended
    if (nRec - 1 - nFailedLigands > 0) {
      auto hTotal =
          std::chrono::duration_cast<std::chrono::hours>(totalDuration);
      totalDuration -= hTotal;
      auto mTotal =
          std::chrono::duration_cast<std::chrono::minutes>(totalDuration);
      totalDuration -= mTotal;

      std::cout << "Docking duration for " << nRec - 1 - nFailedLigands
                << " ligand(s): ";

      if (hTotal.count() > 0) {
        std::cout << hTotal.count() << " hour(s), ";
      }
      if (hTotal.count() > 0 || mTotal.count() > 0) {
        std::cout << mTotal.count() << " minute(s), ";
      }
      std::cout << totalDuration.count() << " second(s)" << std::endl;
    }

    if (nUnnamedLigands > 0) {
      std::cout
          << std::endl
          << "WARNING: " << nUnnamedLigands
          << " ligand(s) are unnamed. Post-processing tools might have an "
             "issue correctly identifying different poses of the same ligand."
          << std::endl;
    }

    std::cout << std::endl << "END OF RUN" << std::endl;
    //    if (bOutput && flexRec) {
    //      MolecularFileSinkPtr spRecepSink(new
    //      CrdFileSink(strRunName+".crd",spReceptor));
    //      spRecepSink->Render();
    //    }
    std::cout << std::endl;
    PrintBibliographyItem(std::cout, "RiboDock2004");
    PrintBibliographyItem(std::cout, "rDock2014");
#if !defined(__sun) && !defined(_MSC_VER)
    PrintBibliographyItem(std::cout, "PCG2014");
#endif
    std::cout << std::endl;
    std::cout << "Thank you for using " << GetProgramName() << " "
              << GetProgramVersion() << "." << std::endl;
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
