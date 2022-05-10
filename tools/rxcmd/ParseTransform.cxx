//===-- ParseTransform.cxx - Parse CLI params for Transform -----*- C++ -*-===//
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
/// Parse command-line interface parameters for Transform operation.
///
//===----------------------------------------------------------------------===//

#include "ParseTransform.h"

#include "rxdock/operation/Transform.h"

#include <cxxopts.hpp>
#include <fmt/format.h>

#include <stdexcept> // TODO

int rxdock::parseTransform(int argc, char *argv[]) {
  cxxopts::Options options("rxcmd transform", "Transform structure-data files");
  options.positional_help("");

  // Command line arguments and default values
  cxxopts::OptionAdder adder = options.add_options();
  adder("i,input", "Input structure-data file (SDfile) name",
        cxxopts::value<std::string>());
  adder("o,output", "Output structure-data file (SDfile) name",
        cxxopts::value<std::string>());
  adder("n,name", "Only work on molecules matching name",
        cxxopts::value<std::string>());
  adder("set-name",
        "Set name to given value (use <field-name> to get value of any data "
        "item)",
        cxxopts::value<std::string>());
  adder("d,data", "Only work on molecules matching field-name=value",
        cxxopts::value<std::string>());
  adder("set-data",
        "Set value of data item to given value (use <field-name> to get value "
        "of any data item)",
        cxxopts::value<std::string>());
  adder("s,sort", "Sort output data by field", cxxopts::value<std::string>());
  adder("l,limit", "Number of records to output",
        cxxopts::value<std::size_t>());
  adder("p,limit-per-molecule", "Number of records to output per molecule",
        cxxopts::value<std::size_t>());
  adder("b,buffer",
        "Number of molecules to read in a render cycle (ignored if sorting is "
        "used)",
        cxxopts::value<std::size_t>()->default_value("1000"));
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

    std::string outputSDFile;
    if (result.count("o")) {
      outputSDFile = result["o"].as<std::string>();
    } else {
      fmt::print("Output molecule SDfile file name is missing.\n");
      return EXIT_FAILURE;
    }
    // FIXME check if file exists

    bool checkName = false;
    std::string name;
    if (result.count("n")) {
      checkName = true;
      name = result["n"].as<std::string>();
    }

    bool changeName = false;
    std::string newName;
    if (result.count("set-name")) {
      changeName = true;
      newName = result["set-name"].as<std::string>();
    }

    // FIXME specified multiple times either with a=b,c=d or with -d foo=bar -d
    // baz=quux
    bool checkData = true;
    std::string dataFieldName;
    std::string dataValue; // FIXME should be std::variant
    if (result.count("d")) {
      checkData = true;
      dataFieldName = result["d"].as<std::string>();
      dataValue = result["d"].as<std::string>();
    }

    bool doSorting = false;
    std::string sortDataField;
    if (result.count("s")) {
      doSorting = true;
      sortDataField = result["s"].as<std::string>();
      if (sortDataField.empty()) {
        fmt::print("Data field name can not be blank.\n");
        return EXIT_FAILURE;
      }
    }

    bool limitedRecords = false;
    std::size_t maxNumRecords = 0;
    if (result.count("l")) {
      limitedRecords = true;
      maxNumRecords = result["l"].as<std::size_t>();
    }

    bool limitedRecPerMolecule = false;
    std::size_t maxNumRecPerMolecule = 0;
    if (result.count("p")) {
      limitedRecPerMolecule = true;
      maxNumRecPerMolecule = result["p"].as<std::size_t>();
    }

    std::size_t bufferSize = result["b"].as<std::size_t>();

    return operation::transform(
        inputSDFile, outputSDFile, limitedRecPerMolecule, maxNumRecPerMolecule,
        checkName, name, checkData, dataFieldName, dataValue, changeName,
        newName, bufferSize, doSorting, sortDataField, limitedRecords,
        maxNumRecords);

  } catch (const cxxopts::OptionException &e) {
    fmt::print("Error parsing options: {}\n", e.what());
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
