//===-- Transform.cxx - Transform operation ---------------------*- C++ -*-===//
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

#include "rxdock/operation/Transform.h"
#include "rxdock/MdlFileSink.h"
#include "rxdock/MdlFileSource.h"
#include "rxdock/ModelError.h"

#include <algorithm>
#include <chrono>
#include <functional>

#include <fmt/chrono.h>
#include <fmt/format.h>

namespace rxdock {
namespace operation {

static void changeModelName(ModelPtr model, std::string newName) {
  std::string formattedName;
  std::ostringstream oss;
  bool inField = false;
  for (auto c : newName) {
    if (c == '<') {
      if (inField) {
        throw std::invalid_argument("Incorrect field name specification <" +
                                    oss.str() + "<");
      }
      inField = true;
      formattedName += oss.str();
      oss.str("");
    } else if (c == '>') {
      if (!inField) {
        throw std::invalid_argument("Incorrect field name specification >" +
                                    oss.str() + ">");
      }
      inField = false;
      formattedName += model->GetDataValue(oss.str()).GetString();
      oss.str("");
    } else
      oss.put(c);
  }
  if (inField) {
    throw std::invalid_argument("Incorrect field name specification <" +
                                oss.str());
  }
  if (oss.tellp() > 0) {
    formattedName += oss.str();
  }
  model->SetTitle(0, formattedName);
}

} // namespace operation
} // namespace rxdock

int rxdock::operation::transform(
    std::string inputSDFile, std::string outputSDFile,
    bool limitedRecPerMolecule, std::size_t maxNumRecPerMolecule,
    bool checkName, std::string name, bool checkData, std::string dataFieldName,
    std::string dataValue, bool changeName, std::string newName,
    std::size_t bufferSize, bool doSorting, std::string sortDataField,
    bool limitedRecords, std::size_t maxNumRecords) {
  try {
    MolecularFileSourcePtr fileSource(new MdlFileSource(inputSDFile));
    MolecularFileSinkPtr fileSink(new MdlFileSink(outputSDFile, ModelPtr()));

    auto startTime = std::chrono::system_clock::now();

    auto dataFieldCmp = [&sortDataField](ModelPtr a, ModelPtr b) {
      Variant va = a->GetDataValue(sortDataField);
      Variant vb = b->GetDataValue(sortDataField);
      // FIXME will work for SCORE, what about the rest?
      return va.GetDouble() < vb.GetDouble();
    };

    std::vector<ModelPtr> ligands;
    std::size_t nLigandsRendered = 0;
    std::size_t nLigandsFailed = 0;
    std::map<std::string, std::vector<ModelPtr>> recsPerMolecule;
    for (std::size_t nRead = 1; fileSource->FileStatusOK();
         fileSource->NextRecord(), nRead++) {
      try {
        ModelPtr model(new Model(fileSource));
        // FIXME reconsider generalizing this assumption
        std::string moleculeName = model->GetTitleList()[0];
        // FIXME: multiple data, correct variant, consider order
        if ((checkName && moleculeName != name) ||
            (checkData &&
             model->GetDataValue(dataFieldName).GetString() != dataValue)) {
          continue;
        }

        if (changeName) {
          changeModelName(model, newName);
        }

        /* if (changeData) */

        if (recsPerMolecule.find(moleculeName) == recsPerMolecule.end()) {
          recsPerMolecule[moleculeName] = std::vector<ModelPtr>{};
          recsPerMolecule[moleculeName].reserve(maxNumRecPerMolecule);
          fmt::print(
              "Picking {} (score {}, {} out of {})\n", model->GetTitleList()[0],
              model->GetDataValue(sortDataField).GetDouble(),
              recsPerMolecule[moleculeName].size(), maxNumRecPerMolecule);
          recsPerMolecule[moleculeName].push_back(model);
        } else if (recsPerMolecule.find(moleculeName) !=
                       recsPerMolecule.end() &&
                   recsPerMolecule[moleculeName].size() <
                       maxNumRecPerMolecule) {
          fmt::print(
              "Picking {} (score {}, {} out of {})\n", model->GetTitleList()[0],
              model->GetDataValue(sortDataField).GetDouble(),
              recsPerMolecule[moleculeName].size(), maxNumRecPerMolecule);
          recsPerMolecule[moleculeName].push_back(model);
          if (doSorting) {
            std::sort(recsPerMolecule[moleculeName].begin(),
                      recsPerMolecule[moleculeName].end(), dataFieldCmp);
          }
        } else {
          Variant worstValue =
              recsPerMolecule[moleculeName]
                             [recsPerMolecule[moleculeName].size() - 1]
                                 ->GetDataValue(sortDataField);
          Variant vModel = model->GetDataValue(sortDataField);
          if (vModel.GetDouble() < worstValue.GetDouble()) {
            auto it = std::lower_bound(recsPerMolecule[moleculeName].begin(),
                                       recsPerMolecule[moleculeName].end(),
                                       model, dataFieldCmp);
            fmt::print("Picking {} (score {}, {} out of {})\n",
                       model->GetTitleList()[0],
                       model->GetDataValue(sortDataField).GetDouble(),
                       recsPerMolecule[moleculeName].size(),
                       maxNumRecPerMolecule);
            recsPerMolecule[moleculeName].insert(it, model);
            recsPerMolecule[moleculeName].pop_back();
          } else {
            fmt::print("Not good enough {} (score {}, worst {})\n",
                       model->GetTitleList()[0], vModel.GetDouble(),
                       worstValue.GetDouble());
          }
        }

        /*if (limitedRecPerMolecule &&
            recsPerMolecule[moleculeName].size() > maxNumRecPerMolecule) {
          recsPerMolecule[moleculeName].pop_back();
        }*/

        /* FIXME simplify, remove later check */
        /*if (doSorting && !ligands.empty() &&
            (limitedRecords &&
             nLigandsRendered + ligands.size() >= maxNumRecords)) {
          double worstValue =
              ligands[ligands.size() - 1]->GetDataValue(sortDataField);
          Variant vModel = model->GetDataValue(sortDataField);
          if (vModel.GetDouble() < worstValue) {
            ligands.pop_back();

            auto it = std::lower_bound(ligands.begin(), ligands.end(), model,
                                       dataFieldCmp);
            fmt::print("Ligand printer goes brrr for {} ({} out of {})\n",
                       model->GetTitleList()[0],
                       recsPerMolecule[moleculeName].size(),
                       maxNumRecPerMolecule);
            ligands.insert(it, model);
          }
        } else {
          fmt::print("Ligand printer goes brrr for {} ({} out of {})\n",
                     model->GetTitleList()[0],
                     recsPerMolecule[moleculeName].size(),
                     maxNumRecPerMolecule);
          ligands.push_back(model);
        }
      }

      if (ligands.size() % bufferSize == 0) {
        if (!doSorting) {
          fmt::print("Will write {} molecules now.\n", nLigandsRendered);
          for (auto &model : ligands) {
            fileSink->SetModel(model);
            fileSink->Render();
          }
          nLigandsRendered += ligands.size();
          ligands.clear();
          fmt::print(
              "Wrote {} molecules, no molecules remain to be written.\n",
              nLigandsRendered);
        }
      }

      if (nRead % bufferSize == 0) {
        std::size_t estNumRecords = fileSource->GetEstimatedNumRecords();
        fmt::print(
            "Read and processed {} molecules out of approximately {}.\n",
            nRead, estNumRecords);

        auto endTime = std::chrono::system_clock::now();
        std::chrono::duration<double> totalDuration = endTime - startTime;
        std::chrono::duration<double> estimatedTimeRemaining =
            estNumRecords * (totalDuration / nRead);
        std::chrono::system_clock::time_point predictedEndTime =
            startTime + std::chrono::duration_cast<std::chrono::seconds>(
                            estimatedTimeRemaining);
        std::time_t predictedEndTimeT =
            std::chrono::system_clock::to_time_t(predictedEndTime);

        auto mTotal =
            std::chrono::duration_cast<std::chrono::minutes>(totalDuration);
        totalDuration -= mTotal;

        if (!doSorting)
          fmt::print("Reading, processing, and writing ");
        else
          fmt::print("Reading and processing ");
        fmt::print("of {} molecules took ", nRead);
        if (mTotal.count() > 0)
          fmt::print("{} minute(s) and ", mTotal.count());
        fmt::print("{} second(s).\n", totalDuration.count());

        fmt::print("Approximately {} molecule(s) remaining, reading and "
                   "processing will finish {:%c}.\n",
                   estNumRecords - nRead,
                   *std::localtime(&predictedEndTimeT));
      }*/
      } catch (ModelError &e) {
        fmt::print("{}", e.what());
        nLigandsFailed++;
      }

      /* FIXME uncomment eventually */
      /* if (limitedRecords && nLigandsRendered + ligands.size() >=
        maxNumRecords) break; */
    }

    if (recsPerMolecule.empty()) {
      fmt::print("No molecules matched, output file will not be written.\n");
      return EXIT_FAILURE;
    }

    for (auto const &ligandsName : recsPerMolecule) {
      for (auto const &ligand : ligandsName.second) {
        ligands.push_back(ligand);
      }
    }

    if (doSorting) {
      std::sort(ligands.begin(), ligands.end(), dataFieldCmp);
    }

    if (nLigandsRendered == 0 && ligands.empty()) {
      // FIXME add matched what
      fmt::print("No molecules matched, output file will not be written.\n");
      return EXIT_FAILURE;
    }

    fmt::print("Will write {} molecules more.\n", nLigandsRendered);
    fmt::print("Will not write {} molecules due to errors (see above).\n",
               nLigandsFailed);
    for (auto &model : ligands) {
      fileSink->SetModel(model);
      fileSink->Render();
      nLigandsRendered++;
      if (nLigandsRendered >= maxNumRecords) {
        break;
      }
    }
    // nLigandsRendered += ligands.size();
    ligands.clear();
    fmt::print("Wrote {} molecules, no molecules remain to be written.\n",
               nLigandsRendered);
  } catch (Error &e) {
    fmt::print("{}", e.what());
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
