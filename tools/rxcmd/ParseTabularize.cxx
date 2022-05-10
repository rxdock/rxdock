//===-- ParseTabularize.cxx - Parse CLI params for Tabularize ---*- C++ -*-===//
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
/// Parse command-line interface parameters for Tabularize operation.
///
//===----------------------------------------------------------------------===//

#include "ParseTabularize.h"

#include "rxdock/operation/Tabularize.h"

#include <cxxopts.hpp>
#include <fmt/format.h>

#include <stdexcept> // TODO

int rxdock::parseTabularize(int argc, char *argv[]) {
  cxxopts::Options options(
      "rxcmd tabularize",
      "Convert a structure-data file into a tabular format");
  options.positional_help("");

  // Command line arguments and default values
  cxxopts::OptionAdder adder = options.add_options();
  adder("i,input", "Input structure-data file (SDfile) name",
        cxxopts::value<std::string>());
  adder("o,output", "Output CSV file name", cxxopts::value<std::string>());
  adder("a,atoms", "Number of atoms to output per molecule",
        cxxopts::value<std::size_t>()->default_value("50"));
  adder("d,data-fields", "Number of data fields to output per molecule",
        cxxopts::value<std::size_t>()->default_value("50"));
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
    std::string inputSDFile;
    if (result.count("i")) {
      inputSDFile = result["i"].as<std::string>();
    } else {
      fmt::print("Input molecule SDfile name is missing.\n");
      return EXIT_FAILURE;
    }
    // FIXME check if file exists

    std::string outputCSVFile;
    if (result.count("o")) {
      outputCSVFile = result["o"].as<std::string>();
    } else {
      fmt::print("Output CSV file name is missing.\n");
      return EXIT_FAILURE;
    }
    // FIXME check if file exists

    std::size_t nAtoms = result["a"].as<std::size_t>();
    std::size_t nDataFields = result["d"].as<std::size_t>();

    return operation::tabularize(inputSDFile, outputCSVFile, nAtoms,
                                 nDataFields);

  } catch (const cxxopts::OptionException &e) {
    fmt::print("Error parsing options: {}\n", e.what());
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
