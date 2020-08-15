//===-- CavitySearch.h - CavitySearch operation -----------------*- C++ -*-===//
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
/// CavitySearch operation.
///
//===----------------------------------------------------------------------===//

#ifndef RXDOCK_OPERATION_CAVITYSEARCH_H
#define RXDOCK_OPERATION_CAVITYSEARCH_H

#include "rxdock/support/Export.h"

#include <string>

namespace rxdock {
namespace operation {

///
/// \brief Finds cavities on the receptor.
///
RBTDLL_EXPORT int cavitySearch(std::string strReceptorPrmFile,
                               bool readDockingSite, bool writeDockingSite,
                               bool writeMOEGrid, bool writeInsightII,
                               bool writePsfCrd, bool listAtoms,
                               double listDistance, bool printSiteDescriptors,
                               double border);

} // namespace operation
} // namespace rxdock

#endif // RXDOCK_OPERATION_CAVITYSEARCH_H
