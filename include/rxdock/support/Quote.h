//===-- Quote.h - Quote reading, selection, and formatting ------*- C++ -*-===//
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
/// Quote reading, selection, and formatting.
///
//===----------------------------------------------------------------------===//

#ifndef RXDOCK_SUPPORT_QUOTE_H
#define RXDOCK_SUPPORT_QUOTE_H

#include "rxdock/support/Export.h"

#include <string>

namespace rxdock {
namespace support {

///
/// \brief Selects a random quote from the database and formats it for printing.
///
/// Will return an empty string if the quotes database can not be opened for
/// reading.
/// \return a random quote formatted for printing
///
RBTDLL_EXPORT std::string getRandomQuote();

} // namespace support
} // namespace rxdock

#endif // RXDOCK_SUPPORT_QUOTE_H
