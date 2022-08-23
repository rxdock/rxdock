/***********************************************************************
 * The rDock program was developed from 1998 - 2006 by the software team
 * at RiboTargets (subsequently Vernalis (R&D) Ltd).
 * In 2006, the software was licensed to the University of York for
 * maintenance and distribution.
 * In 2012, Vernalis and the University of York agreed to release the
 * program as Open Source software.
 * This version is licensed under GNU-LGPL version 3.0 with support from
 * the University of Barcelona.
 * http://rdock.sourceforge.net/
 ***********************************************************************/

#include "rxdock/VdwIntraSF.h"
#include "rxdock/SFRequest.h"

#include <loguru.hpp>

using namespace rxdock;

// Static data members
const std::string VdwIntraSF::_CT = "VdwIntraSF";

// NB - Virtual base class constructor (BaseSF) gets called first,
// implicit constructor for BaseInterSF is called second
VdwIntraSF::VdwIntraSF(const std::string &strName) : BaseSF(_CT, strName) {
  LOG_F(2, "VdwIntraSF parameterised constructor");
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

VdwIntraSF::~VdwIntraSF() {
  LOG_F(2, "VdwIntraSF parameterised constructor");
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

// Request Handling method
// Handles the Partition request
void VdwIntraSF::HandleRequest(RequestPtr spRequest) {
  VariantList params = spRequest->GetParameters();

  switch (spRequest->GetID()) {

    // ID_REQ_SF_PARTITION requests come in two forms:
    // 1-param: param[0] = distance => Partition all scoring functions
    // 2-param: param[0] = SF Name,
    //         param[1] = distance => Partition a named scoring function
  case ID_REQ_SF_PARTITION:
    if (params.size() == 1) {
      LOG_F(1, "VdwIntraSF::HandleRequest: Partitioning {} at distance={}",
            GetFullName(), params[0].GetString());
      Partition(m_ligAtomList, m_vdwIntns, m_prtIntns, params[0]);
    } else if ((params.size() == 2) &&
               (params[0].GetString() == GetFullName())) {
      LOG_F(1, "VdwIntraSF::HandleRequest: Partitioning {} at distance={}",
            GetFullName(), params[1].GetString());
      Partition(m_ligAtomList, m_vdwIntns, m_prtIntns, params[1]);
    }
    break;

    // Pass all other requests to base handler
  default:
    LOG_F(1, "VdwIntraSF::HandleRequest: Partitioning {} at distance={}",
          GetFullName(), params[1].GetString());
    BaseObject::HandleRequest(spRequest);
    break;
  }
}

void VdwIntraSF::SetupScore() {
  LOG_F(2, "VdwIntraSF::SetupScore");
  m_ligAtomList.clear();
  for (AtomRListListIter iter = m_vdwIntns.begin(); iter != m_vdwIntns.end();
       iter++)
    (*iter).clear();
  m_vdwIntns.clear();
  for (AtomRListListIter iter = m_prtIntns.begin(); iter != m_prtIntns.end();
       iter++)
    (*iter).clear();
  m_prtIntns.clear();

  ModelPtr spModel = GetLigand();
  if (spModel.Null())
    return;

  AtomList tmpList = spModel->GetAtomList();
  // Strip off the smart pointers
  std::copy(tmpList.begin(), tmpList.end(), std::back_inserter(m_ligAtomList));

  // Build map of intra-ligand flexible interactions
  m_vdwIntns = AtomRListList(m_ligAtomList.size(), AtomRList());
  m_prtIntns = AtomRListList(m_ligAtomList.size(), AtomRList());
  BuildIntraMap(m_ligAtomList, m_vdwIntns);
  // Partition with zero distance is needed to copy all the vdW interactions
  // into the partitioned list (this is the list that is scored)
  Partition(m_ligAtomList, m_vdwIntns, m_prtIntns, 0.0);
}

double VdwIntraSF::RawScore() const {
  double score = 0.0; // Total score
  // Loop over all ligand atoms
  for (AtomRListConstIter iter = m_ligAtomList.begin();
       iter != m_ligAtomList.end(); iter++) {
    int id = (*iter)->GetAtomId() - 1;
    // XB changed call from "VdwScore" to "VdwScoreIntra" and created new
    // function
    // in "VdwSF.cxx" to avoid using reweighting terms for intra
    // Double s = VdwScoreIntra(*iter,m_prtIntns[id]);
    double s = VdwScore(*iter, m_prtIntns[id]);
    score += s;
  }
  return score;
}

// DM 25 Oct 2000 - track changes to parameter values in local data members
// ParameterUpdated is invoked by ParamHandler::SetParameter
void VdwIntraSF::ParameterUpdated(const std::string &strName) {
  VdwSF::OwnParameterUpdated(strName);
  BaseSF::ParameterUpdated(strName);
}

void rxdock::to_json(json &j, const VdwIntraSF &vdwIntraSF) {
  json atomRListList;
  for (const auto &aListIter : vdwIntraSF.m_vdwIntns) {
    json atomRList;
    for (const auto &aIter : aListIter) {
      json atom = *aIter;
      atomRList.push_back(atom);
    }
    atomRListList.push_back(atomRList);
  }

  json atomRListList2;
  for (const auto &aListIter : vdwIntraSF.m_prtIntns) {
    json atomRList;
    for (const auto &aIter : aListIter) {
      json atom = *aIter;
      atomRList.push_back(atom);
    }
    atomRListList2.push_back(atomRList);
  }

  json atomRList;
  for (const auto &aIter : vdwIntraSF.m_ligAtomList) {
    json atom = *aIter;
    atomRList.push_back(atom);
  }

  j = json{{"vdw-intns", atomRListList},
           {"prt-intns", atomRListList2},
           {"lig-at-li", atomRList}};
}

void rxdock::from_json(const json &j, VdwIntraSF &vdwIntraSF) {
  for (auto &atomList : j.at("vdw-intns")) {
    AtomRList spAtoms;
    for (auto &atom : atomList) {
      AtomPtr spAtom = AtomPtr(new Atom(atom));
      spAtoms.push_back(spAtom);
    }
    vdwIntraSF.m_vdwIntns.push_back(spAtoms);
  }

  for (auto &atomList : j.at("prt-intns")) {
    AtomRList spAtoms;
    for (auto &atom : atomList) {
      AtomPtr spAtom = AtomPtr(new Atom(atom));
      spAtoms.push_back(spAtom);
    }
    vdwIntraSF.m_prtIntns.push_back(spAtoms);
  }

  for (auto &atom : j.at("lig-at-li")) {
    AtomPtr spAtom = AtomPtr(new Atom(atom));
    vdwIntraSF.m_ligAtomList.push_back(spAtom);
  }
}