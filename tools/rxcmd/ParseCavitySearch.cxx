//===-- ParseCavitySearch.cxx - Parse CLI params for CavSearch --*- C++ -*-===//
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
/// Parse command-line interface parameters for CavitySearch operation.
///
//===----------------------------------------------------------------------===//

#include "ParseCavitySearch.h"

#include "rxdock/operation/CavitySearch.h"

#include <cxxopts.hpp>
#include <fmt/format.h>

#include <stdexcept> // TODO

int rxdock::parseCavitySearch(int argc, char *argv[]) {
  cxxopts::Options options("rxcmd cavity-search", "Receptor cavity search");
  options.positional_help("");

  // Command line arguments and default values
  cxxopts::OptionAdder adder = options.add_options();
  adder("r,receptor-param",
        "Receptor parameters file (contains active site parameters)",
        cxxopts::value<std::string>());
  adder("R,read-docking-site",
        "Read docking site (cavities and distance grid) from JSON file");
  adder("W,write-docking-site",
        "Write docking site (cavities and distance grid) to JSON file");
  adder("M,write-moe-grid", "Write active site into a MOE grid");
  adder("I,write-insightii-grids",
        "Write InsightII grids for each cavity for visualisation");
  adder("C,write-psf-crd", "Write target PSF/CRD files for rDock Viewer");
  adder("l,list-atoms-dist",
        "List receptor atoms within specified distance of any cavity (in "
        "angstrom)",
        cxxopts::value<double>()->default_value("5.0"));
  adder("s,print-site", "Print SITE descriptors (counts of exposed atoms)");
  adder(
      "b,border",
      "Set the border around the cavities for the distance grid (in angstrom)",
      cxxopts::value<double>()->default_value("8.0"));
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
    std::string strReceptorPrmFile;
    if (result.count("r")) {
      strReceptorPrmFile = result["r"].as<std::string>();
    } else {
      fmt::print("Receptor parameters file name is missing.\n");
      return EXIT_FAILURE;
    }
    // FIXME check if file exists

    bool readDockingSite = result.count("R");
    bool writeDockingSite = result.count("W");
    bool writeMOEGrid = result.count("M");
    bool writeInsightII = result.count("I");
    bool writePsfCrd = result.count("C");
    bool listAtoms = result.count("l");
    double listDistance = result["l"].as<double>();
    bool printSiteDescriptors = result.count("s");
    double border = result["b"].as<double>();

    return operation::cavitySearch(strReceptorPrmFile, readDockingSite,
                                   writeDockingSite, writeMOEGrid,
                                   writeInsightII, writePsfCrd, listAtoms,
                                   listDistance, printSiteDescriptors, border);

  } catch (const cxxopts::OptionException &e) {
    fmt::print("Error parsing options: {}\n", e.what());
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
