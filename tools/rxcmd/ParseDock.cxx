//===-- ParseDock.cxx - Parse CLI params for Dock ---------------*- C++ -*-===//
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
/// Parse command-line interface parameters for Dock operation.
///
//===----------------------------------------------------------------------===//

#include "ParseDock.h"

#include "rxdock/operation/Dock.h"

#include <cxxopts.hpp>
#include <fmt/format.h>

#include <stdexcept> // TODO

int rxdock::parseDock(int argc, char *argv[]) {
  cxxopts::Options options("rxcmd dock", "Receptor-ligand docking");
  options.positional_help("");

  // Command line arguments and default values
  cxxopts::OptionAdder adder = options.add_options();
  adder("i,input", "Input ligand structure-data file (SDfile) name",
        cxxopts::value<std::string>());
  adder("o,output", "Output docked ligand structure-data file (SDfile) name",
        cxxopts::value<std::string>());
  adder("output-crd", "Output receptor CRD file name",
        cxxopts::value<std::string>());
  adder("output-history", "History structure-data file (SDfile) name prefix",
        cxxopts::value<std::string>());
  adder("r,receptor-param", "Receptor parameters file",
        cxxopts::value<std::string>());
  adder("p,docking-param", "Docking protocol parameters file",
        cxxopts::value<std::string>());
  adder("n,number", "Number of runs per ligand (0 = unlimited)",
        cxxopts::value<std::size_t>()->default_value("50"));
  adder("P,protonate",
        "Protonate all neutral amines, guanidines, and imidazoles");
  adder("D,deprotonate",
        "Deprotonate all carboxylic, sulphur, and phosphorous acid groups");
  adder("H,all-hydrogens",
        "Read all hydrogens present instead of only polar hydrogens");
  adder("t,threshold", "Score threshold", cxxopts::value<double>());
  adder("c,continue",
        "Continue if score threshold is met instead of terminating ligand");
  adder("f,filter", "Filter file name", cxxopts::value<std::string>());
  adder("s,seed", "Random number seed to use instead of std::random_device",
        cxxopts::value<std::size_t>());
  adder("positional",
        "Positional arguments: unused, but useful to have to catch errors",
        cxxopts::value<std::vector<std::string>>());
  adder("h,help", "Print help");

  try {
    options.parse_positional({"positional"});
    auto result = options.parse(argc, argv);

    if (result.count("positional")) {
      fmt::print(
          "Positional arguments are unsupported but were used: {}.\n",
          fmt::join(result["positional"].as<std::vector<std::string>>(), ", "));
      return EXIT_FAILURE;
    }

    if (result.count("h")) {
      fmt::print(options.help());
      return EXIT_SUCCESS;
    }

    // Command line arguments and default values
    std::string strLigandMdlFile;
    if (result.count("i")) {
      strLigandMdlFile = result["i"].as<std::string>();
    } else {
      fmt::print("Input ligand structure-data file name is missing.\n");
      return EXIT_FAILURE;
    }
    // FIXME check if file exists

    std::string strOutputMdlFile;
    if (result.count("o")) {
      strOutputMdlFile = result["o"].as<std::string>();
    } else {
      fmt::print("Output ligand structure-data file name is missing.\n");
      return EXIT_FAILURE;
    }
    // FIXME check if file exists

    bool bOutputCrd = result.count("output-crd");
    std::string strOutputCrdFile;
    if (bOutputCrd) {
      strOutputCrdFile = result["output-crd"].as<std::string>();
    }

    bool bOutputHistory = result.count("output-history");
    std::string strOutputHistoryFilePrefix;
    if (bOutputHistory) {
      strOutputHistoryFilePrefix = result["output-history"].as<std::string>();
    }

    std::string strReceptorPrmFile;
    if (result.count("r")) {
      strReceptorPrmFile = result["r"].as<std::string>();
    } else {
      fmt::print("Receptor parameters file name is missing.\n");
      return EXIT_FAILURE;
    }
    // FIXME check if file exists

    std::string strParamFile;
    if (result.count("p")) {
      strParamFile = result["p"].as<std::string>();
    } else {
      fmt::print("Docking protocol parameters file name is missing.\n");
      return EXIT_FAILURE;
    }
    // FIXME check if file exists

    bool bFilter = result.count("f");
    std::string strFilterFile;
    if (bFilter) {
      strFilterFile = result["f"].as<std::string>();
    }

    bool bDockingRuns = result.count("n");
    std::size_t nDockingRuns = result["n"].as<std::size_t>();

    bool bPosIonise = result.count("P");
    bool bNegIonise = result.count("D");
    bool bExplH = result.count("H");

    bool bTarget = result.count("t");
    double dTargetScore = 0.0;
    if (bTarget) {
      dTargetScore = result["t"].as<double>();
    }
    bool bContinue = result.count("c");

    bool bSeed = result.count("s");
    std::size_t nSeed = 0;
    if (bSeed) {
      nSeed = result["s"].as<std::size_t>();
    }

    operation::dock(strLigandMdlFile, strOutputMdlFile, bOutputCrd,
                    strOutputCrdFile, bOutputHistory,
                    strOutputHistoryFilePrefix, strReceptorPrmFile,
                    strParamFile, bFilter, strFilterFile, bDockingRuns,
                    nDockingRuns, bPosIonise, bNegIonise, bExplH, bTarget,
                    dTargetScore, bContinue, bSeed, nSeed);

  } catch (const cxxopts::OptionException &e) {
    fmt::print("Error parsing options: {}\n", e.what());
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
