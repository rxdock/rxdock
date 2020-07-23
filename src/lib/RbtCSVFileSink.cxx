//===-- RbtCSVFileSink.cxx - RbtModel output to a CSV file ------*- C++ -*-===//
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

#include "RbtCSVFileSink.h"

#include <iomanip>

using namespace rxdock;

RbtCSVFileSink::RbtCSVFileSink(const std::string &fileName, RbtModelPtr spModel,
                               std::size_t nAtoms, std::size_t nDataFields)
    : RbtBaseMolecularFileSink(fileName, spModel), m_nAtoms(nAtoms),
      m_nDataFields(nDataFields) {
  SetAppend(true);
  m_spElementData = RbtElementFileSourcePtr(
      new RbtElementFileSource(GetRbtFileName("data", "RbtElements.dat")));
  _RBTOBJECTCOUNTER_CONSTR_("RbtCSVFileSink");
}

RbtCSVFileSink::~RbtCSVFileSink() {
  _RBTOBJECTCOUNTER_DESTR_("RbtCSVFileSink");
}

void RbtCSVFileSink::Render() {
  try {
    RbtModelPtr spModel(GetModel());
    RbtAtomList modelAtomList(spModel->GetAtomList());
    RbtBondList modelBondList(spModel->GetBondList());
    RbtModelList solventList = GetSolvent();
    // Concatenate all solvent atoms and bonds into a single list
    // DM 7 June 2006 - only render the enabled solvent models
    RbtAtomList solventAtomList;
    RbtBondList solventBondList;
    for (RbtModelListConstIter iter = solventList.begin();
         iter != solventList.end(); ++iter) {
      if ((*iter)->GetEnabled()) {
        RbtAtomList atomList = (*iter)->GetAtomList();
        RbtBondList bondList = (*iter)->GetBondList();
        std::copy(atomList.begin(), atomList.end(),
                  std::back_inserter(solventAtomList));
        std::copy(bondList.begin(), bondList.end(),
                  std::back_inserter(solventBondList));
      }
    }

    std::string line;

    // Write title
    std::vector<std::string> titleList(spModel->GetTitleList());
    if (!titleList.empty())
      line += titleList[0] + ",";
    else
      line += spModel->GetName() + ",";

    // Write atoms
    std::ostringstream oss;
    oss << modelAtomList.size() + solventAtomList.size() << ",";
    line += oss.str();
    line += RenderAtomList(modelAtomList);
    line += RenderAtomList(solventAtomList);

    // Write data fields
    line += RenderData(spModel->GetDataMap());

    // Finish the line
    AddLine(line);

    Write();
  } catch (RbtError &error) {
    // Close();
    throw; // Rethrow the RbtError
  }
}

std::string RbtCSVFileSink::RenderAtomList(const RbtAtomList &atomList) {
  std::string atomLine;
  std::size_t i = 0;
  for (RbtAtomListConstIter aIter = atomList.begin();
       aIter != atomList.end() && i < m_nAtoms; aIter++, i++) {
    RbtAtomPtr spAtom(*aIter);
    RbtElementData elData =
        m_spElementData->GetElementData(spAtom->GetAtomicNo());
    int nFormalCharge = spAtom->GetFormalCharge();
    if (nFormalCharge != 0)
      nFormalCharge = 4 - nFormalCharge;
    std::ostringstream oss;
    oss << elData.element.c_str() << ","; // Element name
    oss << spAtom->GetX() << "," << spAtom->GetY() << "," << spAtom->GetZ()
        << ",";                  // X,Y,Z coord
    oss << nFormalCharge << ","; // charge
    atomLine += oss.str();
  }
  for (; i < m_nAtoms; i++) {
    std::ostringstream oss;
    oss << ",";   // Element name
    oss << ",,,"; // X,Y,Z coord
    oss << ",";   // charge
    atomLine += oss.str();
  }
  return atomLine;
}

std::string RbtCSVFileSink::RenderData(const RbtStringVariantMap &dataMap) {
  std::string dataLine;

  std::ostringstream oss;
  oss << dataMap.size() << ",";
  std::size_t i = 0;
  if (!dataMap.empty()) {
    for (RbtStringVariantMapConstIter iter = dataMap.begin();
         iter != dataMap.end() && i < m_nDataFields; iter++, i++) {
      oss.clear();
      oss << "\"" + (*iter).first + "="; // Field name =
      std::vector<std::string> sl = (*iter).second.StringList();
      for (std::vector<std::string>::const_iterator slIter = sl.begin();
           slIter != sl.end(); ++slIter) {
        oss << *slIter; // Field values
        if (slIter != sl.end() - 1) {
          oss << ",";
        }
      }
      oss << "\",";
      dataLine += oss.str();
    }
  }
  for (; i < m_nDataFields; i++) {
    std::ostringstream oss;
    oss << ","; // Field
    dataLine += oss.str();
  }
  return dataLine;
}
