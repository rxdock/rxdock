//===-- Tabularize.h - Tabularize operation ---------------------*- C++ -*-===//
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
/// Tabularize operation.
///
//===----------------------------------------------------------------------===//

#ifndef RXDOCK_OPERATION_TABULARIZE_H
#define RXDOCK_OPERATION_TABULARIZE_H

#include "rxdock/support/Export.h"

#include <string>

namespace rxdock {
namespace operation {

///
/// \brief Represents structure-data file (SDF) as a CSV table.
///
RBTDLL_EXPORT int tabularize(std::string inputSDFile, std::string outputCSVFile,
                             std::size_t nAtoms, std::size_t nDataFields);

} // namespace operation
} // namespace rxdock

#endif // RXDOCK_OPERATION_TABULARIZE_H
