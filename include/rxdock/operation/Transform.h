//===-- Transform.h - Transform operation -----------------------*- C++ -*-===//
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
/// Transform operation.
///
//===----------------------------------------------------------------------===//

#ifndef RXDOCK_OPERATION_TRANSFORM_H
#define RXDOCK_OPERATION_TRANSFORM_H

#include "rxdock/support/Export.h"

#include <string>

namespace rxdock {
namespace operation {

///
/// \brief Transforms structure-data file (SDF) according to your requirements.
///
RBTDLL_EXPORT int transform(std::string inputSDFile, std::string outputSDFile,
                            bool limitedRecPerMolecule,
                            std::size_t maxNumRecPerMolecule, bool checkName,
                            std::string name, bool checkData,
                            std::string dataFieldName, std::string dataValue,
                            bool changeName, std::string newName,
                            std::size_t bufferSize, bool doSorting,
                            std::string sortDataField, bool limitedRecords,
                            std::size_t maxNumRecords);

} // namespace operation
} // namespace rxdock

#endif // RXDOCK_OPERATION_TRANSFORM_H
