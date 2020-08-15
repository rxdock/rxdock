//===-- Dock.h - Dock operation ---------------------------------*- C++ -*-===//
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
/// Dock operation.
///
//===----------------------------------------------------------------------===//

#ifndef RXDOCK_OPERATION_DOCK_H
#define RXDOCK_OPERATION_DOCK_H

#include "rxdock/support/Export.h"

#include <string>

namespace rxdock {
namespace operation {

///
/// \brief Docks ligand(s) to the receptor.
///
RBTDLL_EXPORT int
dock(std::string strLigandMdlFile, std::string strOutputMdlFile,
     bool bOutputCrd, std::string strOutputCrdFile, bool bOutputHistory,
     std::string strOutputHistoryFilePrefix, std::string strReceptorPrmFile,
     std::string strParamFile, bool bFilter, std::string strFilterFile,
     bool bDockingRuns, std::size_t nDockingRuns, bool bPosIonise,
     bool bNegIonise, bool bExplH, bool bTarget, double dTargetScore,
     bool bContinue, bool bSeed, std::size_t nSeed);

} // namespace operation
} // namespace rxdock

#endif // RXDOCK_OPERATION_DOCK_H
