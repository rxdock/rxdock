//===-- RbtCSVFileSink.h - RbtModel output to a CSV file --------*- C++ -*-===//
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
/// Outputs RbtModel(s) to a CSV file.
///
//===----------------------------------------------------------------------===//

#ifndef _RBTCSVFILESINK_H_
#define _RBTCSVFILESINK_H_

#include "RbtBaseMolecularFileSink.h"
#include "RbtElementFileSource.h"

class RbtCSVFileSink : public RbtBaseMolecularFileSink {
public:
  RBTDLL_EXPORT RbtCSVFileSink(const std::string &fileName, RbtModelPtr spModel,
                               std::size_t nAtoms, std::size_t nDataFields);

  virtual ~RbtCSVFileSink();

  virtual void Render();

private:
  std::string RenderAtomList(const RbtAtomList &atomList);
  std::string RenderData(const RbtStringVariantMap &dataMap);

  RbtCSVFileSink();                       // Disable default constructor
  RbtCSVFileSink(const RbtCSVFileSink &); // Copy constructor disabled by
                                          // default
  RbtCSVFileSink &
  operator=(const RbtCSVFileSink &); // Copy assignment disabled by default

  std::size_t m_nAtoms;
  std::size_t m_nDataFields;
  RbtElementFileSourcePtr m_spElementData;
};

// Useful typedefs
typedef SmartPtr<RbtCSVFileSink> RbtCSVFileSinkPtr; // Smart pointer

#endif //_RBTCSVFILESINK_H_
