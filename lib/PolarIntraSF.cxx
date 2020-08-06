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
std::string PolarIntraSF::_CT("PolarIntraSF");
std::string PolarIntraSF::_ATTR("ATTR");

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
