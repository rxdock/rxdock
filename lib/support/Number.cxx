//===-- Number.cxx - Number conversion --------------------------*- C++ -*-===//
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
/// Number conversion.
///
//===----------------------------------------------------------------------===//

#include "rxdock/support/Number.h"

static std::string getRomanDigit(std::size_t num, char one, char five,
                                 char ten) {
  if (num < 4)
    return std::string(num, one);
  if (num < 6)
    return std::string(5 - num, one) + five;
  if (num < 9)
    return five + std::string(num - 5, one);
  return std::string(10 - num, one) + ten;
}

std::string rxdock::support::getRomanNumeral(std::size_t value) {
  if (value == 0)
    return "nulla";

  std::string RomanNumeral;
  if (value > 1000) {
    RomanNumeral += std::string(value / 1000, 'M');
    value %= 1000;
  }
  if (value > 100) {
    RomanNumeral += getRomanDigit(value / 100, 'C', 'D', 'M');
    value %= 100;
  }
  if (value > 10) {
    RomanNumeral += getRomanDigit(value / 10, 'X', 'L', 'C');
    value %= 10;
  }
  RomanNumeral += getRomanDigit(value, 'I', 'V', 'X');
  return RomanNumeral;
}
