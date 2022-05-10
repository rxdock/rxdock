//===-- Quote.cxx - Quote reading, selection, and formatting ----*- C++ -*-===//
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
/// Quote reading, selection, and formatting.
///
//===----------------------------------------------------------------------===//

#include "rxdock/support/Quote.h"
#include "rxdock/support/Number.h"

#include "rxdock/Rand.h"

#include <nlohmann/json.hpp>

#include <fstream>
#include <sstream>

std::string rxdock::support::getRandomQuote() {
  std::ostringstream qss;

  std::string quotesFileName = GetDataFileName("data/support", "quotes.json");
  std::ifstream quotesFile(quotesFileName.c_str());
  if (quotesFile.is_open()) {
    json quotes;
    quotesFile >> quotes;
    quotesFile.close();

    Rand &rand = GetRandInstance();
    std::size_t quoteNumber = rand.GetRandomInt(quotes.at("quotes").size());

    qss << "Cultural uplift with " << GetProgramName()
        << ", ἐπεισόδιον/episodium " << getRomanNumeral(quoteNumber + 1)
        << ": ";

    std::string quoteText;
    quotes.at("quotes").at(quoteNumber).at("text").get_to(quoteText);
    qss << "\"" << quoteText << "\" ";

    std::string quoteSource;
    quotes.at("quotes").at(quoteNumber).at("source").get_to(quoteSource);
    qss << "(" << quoteSource << ")";
  }
  return qss.str();
}
