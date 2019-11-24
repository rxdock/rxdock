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
#include <iomanip>
#ifdef _WIN32
#else
#include <unistd.h>
#endif

#include "RbtBiMolWorkSpace.h"
#include "RbtCrdFileSink.h"
#include "RbtDockingError.h"
#include "RbtFileError.h"
#include "RbtFilter.h"
#include "RbtLigandError.h"
#include "RbtMdlFileSink.h"
#include "RbtMdlFileSource.h"
#include "RbtModelError.h"
#include "RbtPRMFactory.h"
#include "RbtParameterFileSource.h"
#include "RbtRand.h"
#include "RbtSFFactory.h"
#include "RbtSFRequest.h"
#include "RbtTransformFactory.h"

const std::string EXEVERSION =
    " ($Id: //depot/dev/client3/rdock/2013.1/src/exe/rbdock.cxx#4 $)";
// Section name in docking prm file containing scoring function definition
const std::string _ROOT_SF = "SCORE";
const std::string _RESTRAINT_SF = "RESTR";
const std::string _ROOT_TRANSFORM = "DOCK";

void PrintUsage(void) {
  std::cout << std::endl << "Usage:" << std::endl;
  std::cout << "rbdock -i <sdFile> -o <outputRoot> -r <recepPrmFile> -p "
               "<protoPrmFile> [-n <nRuns>] [-P] [-D] [-H]"
            << std::endl;
  std::cout
      << "       [-t <targetScore|targetFilterFile>] [-c] [-T <traceLevel>] "
         "[-s <rndSeed>]"
      << std::endl;
  std::cout << std::endl
            << "Options:\t-i <sdFile> - input ligand SD file" << std::endl;
  std::cout << "\t\t-o <outputRoot> - root name for output file(s)"
            << std::endl;
  std::cout << "\t\t-r <recepPrmFile> - receptor parameter file " << std::endl;
  std::cout << "\t\t-p <protoPrmFile> - docking protocol parameter file"
            << std::endl;
  std::cout << "\t\t-n <nRuns> - number of runs/ligand (default=1)"
            << std::endl;
  std::cout << "\t\t-P - protonate all neutral amines, guanidines, imidazoles "
               "(default=disabled)"
            << std::endl;
  std::cout
      << "\t\t-D - deprotonate all carboxylic, sulphur and phosphorous acid "
         "groups (default=disabled)"
      << std::endl;
  std::cout
      << "\t\t-H - read all hydrogens present (default=polar hydrogens only)"
      << std::endl;
  std::cout << "\t\t-t - score threshold OR filter file name" << std::endl;
  std::cout << "\t\t-C - continue if score threshold is met (use with -t "
               "<targetScore>, default=terminate ligand)"
            << std::endl;
  std::cout << "\t\t-T <traceLevel> - controls output level for debugging (0 = "
               "minimal, >0 = more verbose)"
            << std::endl;
  std::cout << "\t\t-s <rndSeed> - random number seed (default=from "
               "std::random_device)"
            << std::endl;
}

/////////////////////////////////////////////////////////////////////
// MAIN PROGRAM STARTS HERE
/////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[]) {
  // Handle obsolete arguments, if any
  for (int i = 0; i < argc; i++) {
    std::string opt = argv[i];
    if (opt == "-ap" || opt == "-an" || opt == "-allH" || opt == "-cont") {
      std::cout
          << "Options -ap, -an, -allH, and -cont are no longer supported; use "
             "-P, -D, -H, and -C (respectively) instead."
          << std::endl;
      return 2;
    }
  }

  std::cout.setf(std::ios_base::left, std::ios_base::adjustfield);

  // Strip off the path to the executable, leaving just the file name
  std::string strExeName(argv[0]);
  std::string::size_type i = strExeName.rfind("/");
  if (i != std::string::npos)
    strExeName.erase(0, i + 1);

  // Print a standard header
  Rbt::PrintStdHeader(std::cout, strExeName + EXEVERSION);

  // Brief help message
  if (argc < 2) {
    PrintUsage();
    return 1;
  }

  // Command line arguments and default values
  std::string strLigandMdlFile;
  bool bOutput(false);
  std::string strRunName;
  std::string strReceptorPrmFile; // Receptor param file
  std::string strParamFile;       // Docking run param file
  std::string strFilterFile;      // Filter file
  int nDockingRuns(
      0); // Init to zero, so can detect later whether user explictly typed -n

  // Params for target score
  bool bTarget(false);
  bool bStop(true);         // DM 25 May 2001 - if true, stop once target is met
  bool bDockingRuns(false); // is argument -n present?
  double dTargetScore(0.0);
  bool bFilter(false);

  bool bPosIonise(false);
  bool bNegIonise(false);
  bool bImplH(true); // if true, read only polar hydrogens from SD file, else
                     // read all H's present
  bool bSeed(false); // Random number seed (default = from std::random_device)
  int nSeed(0);
  bool bTrace(false);
  int iTrace(0); // Trace level, for debugging

  double val(0.0);
  char c;
  while ((c = getopt(argc, argv, "i:o:r:p:n:PDHtCT:s:")) != -1) {
    switch (c) {
    case 'i':
      strLigandMdlFile = optarg;
      break;
    case 'o':
      strRunName = optarg;
      break;
    case 'r':
      strReceptorPrmFile = optarg;
      break;
    case 'p':
      strParamFile = optarg;
      break;
    case 'n':
      nDockingRuns = std::atoi(optarg);
      break;
    case 'P': // protonate
      bPosIonise = true;
      break;
    case 'D': // deprotonate
      bNegIonise = true;
      break;
    case 'H': // all-H model
      bImplH = false;
      break;
    case 'C':
      bStop = false;
      break;
    case 't':
      // If str can be translated to an integer, I assume is a
      // threshold. Otherwise, I assume is the filter file name
      char *error;
      errno = 0;
      val = std::strtod(optarg, &error);
      if (!errno && !*error) // Is it a number?
      {
        dTargetScore = val;
        bTarget = true;
      } else // Assume is the filter file name
      {
        bFilter = true;
        strFilterFile = optarg;
      }
      break;
    case 's':
      bSeed = true;
      nSeed = std::atoi(optarg);
      break;
    case 'T':
      bTrace = true;
      iTrace = std::atoi(optarg);
      break;
    default:
      break;
    }
  }
  std::cout << std::endl;

  // print out arguments
  // input ligand file, receptor and parameter is compulsory
  std::cout << std::endl << "Command line args:" << std::endl;
  if (strLigandMdlFile.empty() || strReceptorPrmFile.empty() ||
      strParamFile.empty()) { // if any of them is missing
    std::cout << "Missing required parameter(s)" << std::endl;
    std::exit(1);
  }
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
  if (nDockingRuns >=
      1) { // User typed -n explicitly, so set bDockingRuns to true
    bDockingRuns = true;
    std::cout << " -n " << nDockingRuns << std::endl;
  } else {
    nDockingRuns =
        1; // User didn't type -n explicitly, so fall back to the default of n=1
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
  if (!bImplH) // all-H
    std::cout << " -allH " << std::endl;
  if (!bStop) // stop after target
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
      } else          // -t<TS> -n<N> need to check if -cont present
                      // for all other cases it doesn't matter
          if (!bStop) // -t<TS> -n<N> -cont
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
    std::cout << "Automatically protonating positive ionisable groups (amines, "
                 "imidazoles, guanidines)"
              << std::endl;
  if (bNegIonise)
    std::cout
        << "Automatically deprotonating negative ionisable groups (carboxylic "
           "acids, phosphates, sulphates, sulphonates)"
        << std::endl;
  if (bImplH)
    std::cout << "Reading polar hydrogens only from ligand SD file"
              << std::endl;
  else
    std::cout << "Reading all hydrogens from ligand SD file" << std::endl;

  if (bTarget) {
    std::cout << std::endl
              << "Lower target intermolecular score = " << dTargetScore
              << std::endl;
  }

  try {
    // Create a bimolecular workspace
    RbtBiMolWorkSpacePtr spWS(new RbtBiMolWorkSpace());
    // Set the workspace name to the root of the receptor .prm filename
    std::vector<std::string> componentList =
        Rbt::ConvertDelimitedStringToList(strReceptorPrmFile, ".");
    std::string wsName = componentList.front();
    spWS->SetName(wsName);

    // Read the docking protocol parameter file
    RbtParameterFileSourcePtr spParamSource(new RbtParameterFileSource(
        Rbt::GetRbtFileName("data/scripts", strParamFile)));
    // Read the receptor parameter file
    RbtParameterFileSourcePtr spRecepPrmSource(new RbtParameterFileSource(
        Rbt::GetRbtFileName("data/receptors", strReceptorPrmFile)));
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
    //    INTER    RbtInterSF.prm
    //    INTRA RbtIntraSF.prm
    // END_SECTION
    //
    // Notes:
    // Section name must be SCORE. This is also the name of the root SF
    // aggregate An aggregate is created for each parameter in the section.
    // Parameter name becomes the name of the subaggregate (e.g. SCORE.INTER)
    // Parameter value is the file name for the subaggregate definition
    // Default directory is $RBT_ROOT/data/sf
    RbtSFFactoryPtr spSFFactory(
        new RbtSFFactory()); // Factory class for scoring functions
    RbtSFAggPtr spSF(new RbtSFAgg(_ROOT_SF)); // Root SF aggregate
    spParamSource->SetSection(_ROOT_SF);
    std::vector<std::string> sfList(spParamSource->GetParameterList());
    // Loop over all parameters in the SCORE section
    for (std::vector<std::string>::const_iterator sfIter = sfList.begin();
         sfIter != sfList.end(); sfIter++) {
      // sfFile = file name for scoring function subaggregate
      std::string sfFile(Rbt::GetRbtFileName(
          "data/sf", spParamSource->GetParameterValueAsString(*sfIter)));
      RbtParameterFileSourcePtr spSFSource(new RbtParameterFileSource(sfFile));
      // Create and add the subaggregate
      spSF->Add(spSFFactory->CreateAggFromFile(spSFSource, *sfIter));
    }

    // Add the RESTRAINT subaggregate scoring function from any SF definitions
    // in the receptor prm file
    spSF->Add(spSFFactory->CreateAggFromFile(spRecepPrmSource, _RESTRAINT_SF));

    // Create the docking transform aggregate from the transform definitions in
    // the docking prm file
    RbtTransformFactoryPtr spTransformFactory(new RbtTransformFactory());
    spParamSource->SetSection();
    RbtTransformAggPtr spTransform(
        spTransformFactory->CreateAggFromFile(spParamSource, _ROOT_TRANSFORM));

    // Override the TRACE levels for the scoring function and transform
    // Dump details to std::cout
    // Register the scoring function and the transform with the workspace
    if (bTrace) {
      RbtRequestPtr spTraceReq(new RbtSFSetParamRequest("TRACE", iTrace));
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
    RbtVariant vLib(Rbt::GetProduct() + " (" + Rbt::GetVersion() + ", Build" +
                    Rbt::GetBuild() + ")");
    RbtVariant vExe(strExeName + EXEVERSION);
    RbtVariant vRecep(spRecepPrmSource->GetFileName());
    RbtVariant vPrm(spParamSource->GetFileName());
    RbtVariant vDir(Rbt::GetCurrentWorkingDirectory());

    spRecepPrmSource->SetSection();
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
    // DM 14 June 2006 - bug fix to one of the longest standing rDock issues
    //(the cryptic "Error reading from input stream" message, if cavity file was
    // missing)
    if (!istr) {
      std::string message = "Cavity file (" + strASFile +
                            ") not found in current directory or $RBT_HOME";
      message += " - run rbcavity first";
      throw RbtFileReadError(_WHERE_, message);
    }
    RbtDockingSitePtr spDS(new RbtDockingSite(istr));
    istr.close();
    spWS->SetDockingSite(spDS);
    std::cout << std::endl
              << "DOCKING SITE" << std::endl
              << (*spDS) << std::endl;

    // Prepare the SD file sink for saving the docked conformations for each
    // ligand DM 3 Dec 1999 - replaced ostrstream with RbtString in determining
    // SD file name SRC 2014 moved here this block to allow WRITE_ERROR TRUE
    if (bOutput) {
      RbtMolecularFileSinkPtr spMdlFileSink(
          new RbtMdlFileSink(strRunName + ".sd", RbtModelPtr()));
      spWS->SetSink(spMdlFileSink);
    }

    RbtPRMFactory prmFactory(spRecepPrmSource, spDS);
    prmFactory.SetTrace(iTrace);
    // Create the receptor model from the file names in the receptor parameter
    // file
    RbtModelPtr spReceptor = prmFactory.CreateReceptor();
    spWS->SetReceptor(spReceptor);

    // Register any solvent
    RbtModelList solventList = prmFactory.CreateSolvent();
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
    RbtRand &theRand = Rbt::GetRbtRand(); // ref to random number generator
    if (bSeed) {
      theRand.Seed(nSeed);
    }

    // Create the filter object for controlling early termination of protocol
    RbtFilterPtr spfilter;
    if (bFilter) {
      spfilter = new RbtFilter(strFilterFile);
      if (bDockingRuns) {
        spfilter->SetMaxNRuns(nDockingRuns);
      }
    } else {
      spfilter = new RbtFilter(strFilter.str(), true);
    }
    if (bTrace) {
      RbtRequestPtr spTraceReq(new RbtSFSetParamRequest("TRACE", iTrace));
      spfilter->HandleRequest(spTraceReq);
    }

    // Register the Filter with the workspace
    spWS->SetFilter(spfilter);

    // MAIN LOOP OVER LIGAND RECORDS
    // DM 20 Apr 1999 - add explicit bPosIonise and bNegIonise flags to
    // MdlFileSource constructor
    RbtMolecularFileSourcePtr spMdlFileSource(
        new RbtMdlFileSource(strLigandMdlFile, bPosIonise, bNegIonise, bImplH));
    for (int nRec = 1; spMdlFileSource->FileStatusOK();
         spMdlFileSource->NextRecord(), nRec++) {
      std::cout.setf(std::ios_base::left, std::ios_base::adjustfield);
      std::cout << std::endl
                << "**************************************************"
                << std::endl
                << "RECORD #" << nRec << std::endl;
      RbtError molStatus = spMdlFileSource->Status();
      if (!molStatus.isOK()) {
        std::cout << std::endl
                  << molStatus << std::endl
                  << "************************************************"
                  << std::endl;
        continue;
      }

      // DM 26 Jul 1999 - only read the largest segment (guaranteed to be called
      // H) BGD 07 Oct 2002 - catching errors created by the ligands, so rbdock
      // continues with the next one, instead of completely stopping
      try {
        spMdlFileSource->SetSegmentFilterMap(
            Rbt::ConvertStringToSegmentMap("H"));

        if (spMdlFileSource->isDataFieldPresent("Name"))
          std::cout << "NAME:   " << spMdlFileSource->GetDataValue("Name")
                    << std::endl;
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
        RbtModelPtr spLigand = prmFactory.CreateLigand(spMdlFileSource);
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
        // need to check this here. The termination
        // filter is only run once at least
        // one docking run has been done.
        if (nDockingRuns < 1)
          bTargetMet = true;
        while (!bTargetMet) {
          // Catching errors with this specific run
          try {
            if (bOutput) {
              std::ostringstream histr;
              histr << strRunName << "_" << strMolName << nRec << "_his_"
                    << iRun << ".sd";
              RbtMolecularFileSinkPtr spHistoryFileSink(
                  new RbtMdlFileSink(histr.str(), spLigand));
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
          } catch (RbtDockingError &e) {
            std::cout << e << std::endl;
          }
        }
        // END OF MAIN LOOP OVER EACH SIMULATED ANNEALING RUN
        ////////////////////////////////////////////////////
      }
      // END OF TRY
      catch (RbtLigandError &e) {
        std::cout << e << std::endl;
      }
    }
    // END OF MAIN LOOP OVER LIGAND RECORDS
    ////////////////////////////////////////////////////
    std::cout << std::endl << "END OF RUN" << std::endl;
    //    if (bOutput && flexRec) {
    //      RbtMolecularFileSinkPtr spRecepSink(new
    //      RbtCrdFileSink(strRunName+".crd",spReceptor));
    //      spRecepSink->Render();
    //    }
    std::cout << std::endl;
    Rbt::PrintBibliographyItem(std::cout, "RiboDock2004");
    Rbt::PrintBibliographyItem(std::cout, "rDock2014");
#if !defined(__sun) && !defined(_MSC_VER)
    Rbt::PrintBibliographyItem(std::cout, "PCG2014");
#endif
    std::cout << std::endl;
    std::cout << "Thank you for using " << Rbt::GetProgramName() << " "
              << Rbt::GetVersion() << "." << std::endl;
  } catch (RbtError &e) {
    std::cout << e << std::endl;
  } catch (...) {
    std::cout << "Unknown exception" << std::endl;
  }

  _RBTOBJECTCOUNTER_DUMP_(std::cout)

  return 0;
}
