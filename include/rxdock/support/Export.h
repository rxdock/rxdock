//===-- Export.h - Symbol export macros -------------------------*- C++ -*-===//
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
/// Symbol export macros.
///
//===----------------------------------------------------------------------===//

#ifndef RXDOCK_SUPPORT_EXPORT_H
#define RXDOCK_SUPPORT_EXPORT_H

#if defined _WIN32 || defined __CYGWIN__
#ifdef __GNUC__
#define RBTDLL_EXPORT __attribute__((dllexport))
#define RBTDLL_IMPORT __attribute__((dllimport))
#else
#define RBTDLL_EXPORT __declspec(dllexport)
#define RBTDLL_IMPORT __declspec(dllimport)
#endif
#define RBTDLL_LOCAL
#else
#if __GNUC__ >= 4
#define RBTDLL_EXPORT __attribute__((visibility("default")))
#define RBTDLL_IMPORT __attribute__((visibility("default")))
#define RBTDLL_LOCAL __attribute__((visibility("hidden")))
#else
#define RBTDLL_EXPORT
#define RBTDLL_IMPORT
#define RBTDLL_LOCAL
#endif
#endif

// Same rules apply to Loguru logging library
#define LOGURU_EXPORT RBTDLL_EXPORT

#endif // RXDOCK_SUPPORT_EXPORT_H
