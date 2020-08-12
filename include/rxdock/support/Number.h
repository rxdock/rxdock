//===-- Number.h - Number conversion ----------------------------*- C++ -*-===//
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
/// Number conversion.
///
//===----------------------------------------------------------------------===//

#ifndef RXDOCK_SUPPORT_NUMBER_H
#define RXDOCK_SUPPORT_NUMBER_H

#include "rxdock/support/Export.h"

#include <string>

namespace rxdock {
namespace support {

///
/// \brief Converts a numeric value to the corresponding Roman numeral.
/// \param value a numeric value.
/// \return the corresponding Roman numeral
///
RBTDLL_EXPORT std::string getRomanNumeral(std::size_t value);

} // namespace support
} // namespace rxdock

#endif // RXDOCK_SUPPORT_NUMBER_H
