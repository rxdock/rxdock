//===-- CSVFileSink.h - Model output to a CSV file --------------*- C++ -*-===//
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
/// Outputs Model(s) to a CSV file.
///
//===----------------------------------------------------------------------===//

#ifndef _RBTCSVFILESINK_H_
#define _RBTCSVFILESINK_H_

#include "rxdock/BaseMolecularFileSink.h"
#include "rxdock/ElementFileSource.h"

namespace rxdock {

class CSVFileSink : public BaseMolecularFileSink {
public:
  RBTDLL_EXPORT CSVFileSink(const std::string &fileName, ModelPtr spModel,
                            std::size_t nAtoms, std::size_t nDataFields);

  virtual ~CSVFileSink();

  virtual void Render();

private:
  std::string RenderAtomList(const AtomList &atomList);
  std::string RenderData(const StringVariantMap &dataMap);

  CSVFileSink();                    // Disable default constructor
  CSVFileSink(const CSVFileSink &); // Copy constructor disabled by
                                    // default
  CSVFileSink &
  operator=(const CSVFileSink &); // Copy assignment disabled by default

  std::size_t m_nAtoms;
  std::size_t m_nDataFields;
  ElementFileSourcePtr m_spElementData;
};

// Useful typedefs
typedef SmartPtr<CSVFileSink> CSVFileSinkPtr; // Smart pointer

} // namespace rxdock

#endif //_RBTCSVFILESINK_H_
