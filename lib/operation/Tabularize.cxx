//===-- Tabularize.cxx - Tabularize operation -------------------*- C++ -*-===//
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
/// Tabularize operation.
///
//===----------------------------------------------------------------------===//

#include "rxdock/operation/Tabularize.h"
#include "rxdock/CSVFileSink.h"
#include "rxdock/MdlFileSource.h"
#include "rxdock/ModelError.h"

#include <fmt/format.h>

int rxdock::operation::tabularize(std::string inputSDFile,
                                  std::string outputCSVFile, std::size_t nAtoms,
                                  std::size_t nDataFields) {
  try {
    MolecularFileSourcePtr fileSource(new MdlFileSource(inputSDFile));
    MolecularFileSinkPtr fileSink(
        new CSVFileSink(outputCSVFile, ModelPtr(), nAtoms, nDataFields));

    std::map<std::string, std::size_t> numRecPerMolecule;
    for (std::size_t nRead = 1; fileSource->FileStatusOK();
         fileSource->NextRecord(), nRead++) {
      try {
        ModelPtr model(new Model(fileSource));
        fileSink->SetModel(model);
        fileSink->Render();
      } catch (ModelError &e) {
        fmt::print("{}", e.what());
      }
    }

  } catch (Error &e) {
    fmt::print("{}", e.what());
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
