//===-- ParseCavitySearch.h - Parse CLI params for CavSearch ----*- C++ -*-===//
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
/// Parse command-line interface parameters for CavitySearch operation.
///
//===----------------------------------------------------------------------===//

#ifndef RXDOCK_TOOLS_RXCMD_PARSECAVITYSEARCH_H
#define RXDOCK_TOOLS_RXCMD_PARSECAVITYSEARCH_H

namespace rxdock {

int parseCavitySearch(int argc, char *argv[]);

} // namespace rxdock

#endif // RXDOCK_TOOLS_RXCMD_PARSECAVITYSEARCH_H
