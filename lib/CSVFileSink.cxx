//===-- CSVFileSink.cxx - Model output to a CSV file ------------*- C++ -*-===//
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

#include "rxdock/CSVFileSink.h"

#include <iomanip>

using namespace rxdock;

CSVFileSink::CSVFileSink(const std::string &fileName, ModelPtr spModel,
                         std::size_t nAtoms, std::size_t nDataFields)
    : BaseMolecularFileSink(fileName, spModel), m_nAtoms(nAtoms),
      m_nDataFields(nDataFields) {
  SetAppend(true);
  m_spElementData = ElementFileSourcePtr(
      new ElementFileSource(GetDataFileName("data", "elements.json")));
  _RBTOBJECTCOUNTER_CONSTR_("CSVFileSink");
}

CSVFileSink::~CSVFileSink() { _RBTOBJECTCOUNTER_DESTR_("CSVFileSink"); }

void CSVFileSink::Render() {
  try {
    ModelPtr spModel(GetModel());
    AtomList modelAtomList(spModel->GetAtomList());
    BondList modelBondList(spModel->GetBondList());
    ModelList solventList = GetSolvent();
    // Concatenate all solvent atoms and bonds into a single list
    // DM 7 June 2006 - only render the enabled solvent models
    AtomList solventAtomList;
    BondList solventBondList;
    for (ModelListConstIter iter = solventList.begin();
         iter != solventList.end(); ++iter) {
      if ((*iter)->GetEnabled()) {
        AtomList atomList = (*iter)->GetAtomList();
        BondList bondList = (*iter)->GetBondList();
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
  } catch (Error &error) {
    // Close();
    throw; // Rethrow the Error
  }
}

std::string CSVFileSink::RenderAtomList(const AtomList &atomList) {
  std::string atomLine;
  std::size_t i = 0;
  for (AtomListConstIter aIter = atomList.begin();
       aIter != atomList.end() && i < m_nAtoms; aIter++, i++) {
    AtomPtr spAtom(*aIter);
    ElementData elData = m_spElementData->GetElementData(spAtom->GetAtomicNo());
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

std::string CSVFileSink::RenderData(const StringVariantMap &dataMap) {
  std::string dataLine;

  std::ostringstream oss;
  oss << dataMap.size() << ",";
  std::size_t i = 0;
  if (!dataMap.empty()) {
    for (StringVariantMapConstIter iter = dataMap.begin();
         iter != dataMap.end() && i < m_nDataFields; iter++, i++) {
      oss.clear();
      oss << "\"" + (*iter).first + "="; // Field name =
      std::vector<std::string> sl = (*iter).second.GetStringList();
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
