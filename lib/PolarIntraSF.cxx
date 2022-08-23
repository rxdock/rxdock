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

#include "rxdock/PolarIntraSF.h"
#include "rxdock/SFRequest.h"

#include <loguru.hpp>

using namespace rxdock;

// Static data members
const std::string PolarIntraSF::_CT = "PolarIntraSF";
const std::string PolarIntraSF::_ATTR = "attractive";

// NB - Virtual base class constructor (BaseSF) gets called first,
// implicit constructor for BaseInterSF is called second
PolarIntraSF::PolarIntraSF(const std::string &strName)
    : BaseSF(_CT, strName), m_bAttr(true) {
  LOG_F(2, "PolarIntraSF parameterised constructor");
  // Add parameters
  AddParameter(_ATTR, m_bAttr);
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

PolarIntraSF::~PolarIntraSF() {
  LOG_F(2, "PolarIntraSF destructor");
  ClearModel();
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

void PolarIntraSF::SetupScore() {
  ClearModel();
  ModelPtr spModel = GetLigand();
  if (spModel.Null())
    return;

  AtomList atomList(spModel->GetAtomList());
  m_posList = CreateDonorInteractionCenters(atomList);
  m_negList = CreateAcceptorInteractionCenters(atomList);

  int nAtoms = atomList.size();
  m_intns = InteractionListMap(nAtoms, InteractionCenterList());
  m_prtIntns = InteractionListMap(nAtoms, InteractionCenterList());
  if (m_bAttr) {
    BuildIntraMap(m_posList, m_negList, m_intns);
  } else {
    BuildIntraMap(m_posList, m_intns);
    BuildIntraMap(m_negList, m_intns);
  }
  // Partition with zero distance is needed to copy all the polar interactions
  // into the partitioned list (this is the list that is scored)
  Partition(m_posList, m_negList, m_intns, m_prtIntns, 0.0);
}

double PolarIntraSF::RawScore() const {
  return IntraScore(m_posList, m_negList, m_prtIntns, m_bAttr);
}

void PolarIntraSF::ClearModel() {
  // Clear the interaction maps
  for (InteractionListMapIter iter = m_intns.begin(); iter != m_intns.end();
       iter++) {
    (*iter).clear();
  }
  m_intns.clear();
  for (InteractionListMapIter iter = m_prtIntns.begin();
       iter != m_prtIntns.end(); iter++) {
    (*iter).clear();
  }
  m_prtIntns.clear();

  // Delete the ligand interaction centers
  for (InteractionCenterListIter iter = m_posList.begin();
       iter != m_posList.end(); iter++) {
    delete *iter;
  }
  m_posList.clear();
  for (InteractionCenterListIter iter = m_negList.begin();
       iter != m_negList.end(); iter++) {
    delete *iter;
  }
  m_negList.clear();
}

// DM 25 Oct 2000 - track changes to parameter values in local data members
// ParameterUpdated is invoked by ParamHandler::SetParameter
void PolarIntraSF::ParameterUpdated(const std::string &strName) {
  // DM 25 Oct 2000 - heavily used params
  if (strName == _ATTR) {
    m_bAttr = GetParameter(_ATTR);
  } else {
    PolarSF::OwnParameterUpdated(strName);
    BaseSF::ParameterUpdated(strName);
  }
}

// Request Handling method
// Handles the Partition request
void PolarIntraSF::HandleRequest(RequestPtr spRequest) {
  VariantList params = spRequest->GetParameters();

  switch (spRequest->GetID()) {
    // ID_REQ_SF_PARTITION requests come in two forms:
    // 1-param: param[0] = distance => Partition all scoring functions
    // 2-param: param[0] = SF Name,
    //         param[1] = distance => Partition a named scoring function
  case ID_REQ_SF_PARTITION:
    if (params.size() == 1) {
      LOG_F(1, "PolarIntraSF::HandleRequest: Partitioning {} at distance={}",
            GetFullName(), params[0].GetString());
      Partition(m_posList, m_negList, m_intns, m_prtIntns, params[0]);
    } else if ((params.size() == 2) &&
               (params[0].GetString() == GetFullName())) {
      LOG_F(1, "PolarIntraSF::HandleRequest: Partitioning {} at distance={}",
            GetFullName(), params[1].GetString());
      Partition(m_posList, m_negList, m_intns, m_prtIntns, params[1]);
    }
    break;

    // Pass all other requests to base handler
  default:
    LOG_F(1, "PolarIntraSF::HandleRequest: {} passing request to base handler",
          GetFullName());
    BaseObject::HandleRequest(spRequest);
    break;
  }
}

void rxdock::to_json(json &j, const PolarIntraSF &polarIntraSF) {
  json intersectionCenterList3;
  for (const auto &aIter : polarIntraSF.m_posList) {
    json interactionCenter = *aIter;
    intersectionCenterList3.push_back(interactionCenter);
  }

  json intersectionCenterList4;
  for (const auto &aIter : polarIntraSF.m_negList) {
    json interactionCenter = *aIter;
    intersectionCenterList4.push_back(interactionCenter);
  }

  json intersectionListMap;
  for (const auto &aListIter : polarIntraSF.m_intns) {
    json intersectionCenterList;
    for (const auto &aIter : aListIter) {
      json interactionCenter = *aIter;
      intersectionCenterList.push_back(interactionCenter);
    }
    intersectionListMap.push_back(intersectionCenterList);
  }

  json intersectionListMap2;
  for (const auto &aListIter : polarIntraSF.m_prtIntns) {
    json intersectionCenterList;
    for (const auto &aIter : aListIter) {
      json interactionCenter = *aIter;
      intersectionCenterList.push_back(interactionCenter);
    }
    intersectionListMap2.push_back(intersectionCenterList);
  }

  j = json{{"pos-list", intersectionCenterList3},
           {"neg-list", intersectionCenterList4},
           {"intns", intersectionListMap},
           {"prt-intns", intersectionListMap2},
           {"attr", polarIntraSF.m_bAttr}};
}
/*
void rxdock::from_json(const json &j, PolarIntraSF &polarIntraSF) {
  for (auto &interactionCenter : j.at("pos-list")) {
    InteractionCenter *spInteractionCenter = AtomPtr(new
Atom(interactionCenter)); model.m_atomList.push_back(spAtom);
  }

  for (auto &atom : j.at("atom-list")) {
    AtomPtr spAtom = AtomPtr(new Atom(atom));
    model.m_atomList.push_back(spAtom);
  }

  for (auto &atomList : j.at("ring-list")) {
    AtomList spAtoms;
    for (auto &atom : atomList) {
      AtomPtr spAtom = AtomPtr(new Atom(atom));
      spAtoms.push_back(spAtom);
    }
    model.m_ringList.push_back(spAtoms);
  }
  for (auto &atomList : j.at("ring-list")) {
    AtomList spAtoms;
    for (auto &atom : atomList) {
      AtomPtr spAtom = AtomPtr(new Atom(atom));
      spAtoms.push_back(spAtom);
    }
    model.m_ringList.push_back(spAtoms);
  }
  j.at("attr").get_to(polarIntraSF.m_bAttr);*/