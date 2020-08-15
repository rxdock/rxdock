//===-- rxcmd.cpp - RxDock command-line interface ---------------*- C++ -*-===//
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
/// The starting point of the RxDock command-line interface.
///
//===----------------------------------------------------------------------===//

#include "ParseCavitySearch.h"
#include "ParseDock.h"
#include "ParseTabularize.h"
#include "ParseTransform.h"

#include "rxdock/Rbt.h"
#include "rxdock/support/Quote.h"

#include <cxxopts.hpp>
#include <fmt/chrono.h>
#include <fmt/format.h>
#include <loguru.hpp>

#include <algorithm>
#include <chrono>
#include <functional>
#include <iomanip>
#include <random>
#include <stdexcept>

using namespace rxdock;

namespace rxdock {

///
/// \brief Null function for RxDock commands that aren't yet implemented.
///
int parseNull(int argc, char *argv[]) {
  fmt::print("Command is not implemented yet.\n");
  return EXIT_FAILURE;
}

} // namespace rxdock

///
/// \brief Entry point for the RxDock command-line interface.
///
int main(int argc, char *argv[]) {
  // Set up logging verbosity
  loguru::g_stderr_verbosity = -1; // show ERRORs (-2) and WARNINGs (-1)
  // Parse -v verbosity
  int argcPreLoguru = argc;
  loguru::init(argc, argv);
  fmt::print("Current log verbosity is {}. Supported log verbosities are: \n"
             "-2 (ERROR), -1 (adds WARNING), 0 (adds INFO), 1 (adds DEBUG), "
             "and 2 (adds TRACE)\n",
             loguru::g_stderr_verbosity);
  // Let the user know there is an option to use -v if it wasn't specified
  if (argc == argcPreLoguru) {
    fmt::print("Log verbosity can be set using parameter -v and the desired "
               "value before the command name, e.g.\n"
               "rxcmd -v 1 dock -i input.sd -o output.sd -n 100 ...\n");
  }

  fmt::print(
      "\nWARNING: THE INTERFACE AND THE OUTPUT OF THIS COMMAND IS NOT YET "
      "STABLE.\n");
  fmt::print("Found an issue? Please report it here: "
             "https://gitlab.com/rxdock/rxdock/issues\n");

  fmt::print("\n(({} version {}\n\n", GetProgramName(), GetProgramVersion());

  fmt::print("{}\n", GetCopyright());

  fmt::print("\n{} command-line interface -- rxcmd -- version {}\n\n",
             GetProgramName(), GetProgramVersion());

  fmt::print("Executable:   {}\n", argv[0]);
  fmt::print("Data prefix:  {}\n", GetRoot());
  fmt::print("Working dir:  {}\n", GetCurrentWorkingDirectory());

  std::vector<std::string> args;
  for (int i = 0; i < argc; i++) {
    args.push_back(static_cast<std::string>(argv[i]));
  }

  fmt::print("Command line: {}\n", fmt::join(args, " "));

  std::map<std::string, std::function<int(int, char **)>> commandFunctions = {
      {"dock", parseDock},                  // rbdock
      {"cavity-search", parseCavitySearch}, // rbcavity
      {"grid", parseNull},   // rbcalcgrid, make_grid.csh, rbconvgrid, rbmoegrid
      {"tether", parseNull}, // sdtether
      {"rmsd", parseNull},   // rbrms, sdrmsd
      {"subset", parseNull}, // representative subset for predictions
      {"predict", parseNull},  // rbhtfinder
      {"describe", parseNull}, // rblist (add molecular descriptors and
                               // fingerprints to structure-data files)
      {"transform",
       parseTransform}, // sdfilter, sdsort, sdsplit, sdmodify, sdfield
      {"tabularize", parseTabularize}, // sdreport-like
      {"filter", parseNull},           // merge into transform
      {"sort", parseNull},             // merge into transform
      {"split", parseNull},            // merge into transform
      {"modify", parseNull},           // merge into transform
      {"field", parseNull},            // merge into transform
      {"report",
       parseNull}}; // sdreport for CSV/TSV + user-friendly report to rST/MD

  std::vector<std::string> availableCommands;
  for (const auto &command : commandFunctions) {
    availableCommands.push_back(command.first);
  }

  std::vector<std::string> matchingCommands;
  if (argc > 1) {
    std::string argCommand = argv[1];
    for (const auto &command : commandFunctions) {
      if (command.first.rfind(argCommand, 0) == 0) {
        matchingCommands.push_back(command.first);
      }
    }
    if (matchingCommands.empty()) {
      fmt::print("Unknown command {}, available commands: {}.\n", argCommand,
                 fmt::join(availableCommands, ", "));
      return EXIT_FAILURE;
    } else if (matchingCommands.size() > 1) {
      fmt::print("Ambiguous command {}, matching commands: {}.\n", argCommand,
                 fmt::join(matchingCommands, ", "));
      return EXIT_FAILURE;
    }
  } else {
    fmt::print("Unspecified command, available commands: {}.\n",
               fmt::join(availableCommands, ", "));
    return EXIT_FAILURE;
  }
  std::string matchingCommand = matchingCommands[0];

  // construct argv for subcommand: skip 0 since that's the main command and
  // start from 1 which is a matching command
  args[1] = matchingCommand;

  std::vector<char *> commandArgs;
  /* std::transform(args.begin() + 1, args.end(),
     std::back_inserter(commandArgs),
                 [](std::string s) -> char * {
                   const char *c = s.c_str();
                   return const_cast<char *>(c);
                 }); */
  for (auto it = args.begin() + 1; it < args.end(); it++) {
    const char *arg = (*it).c_str();
    commandArgs.push_back(const_cast<char *>(arg));
  }

  std::function<int(int, char **)> selectedCommandFunction =
      commandFunctions[matchingCommand];
  int exitCode = selectedCommandFunction(argc - 1, commandArgs.data());

  if (exitCode != EXIT_SUCCESS) {
    fmt::print("\nRxDock {} command failed, please look at the error messages "
               "above.\n",
               matchingCommand);
    return exitCode;
  }

  fmt::print("\nRxDock {} command finished successfully.\n", matchingCommand);

  std::string quote = support::getRandomQuote();
  if (!quote.empty()) {
    fmt::print("\n{}\n", quote);
  }

  return EXIT_SUCCESS;
}
