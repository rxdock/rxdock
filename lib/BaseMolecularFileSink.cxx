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

#include "rxdock/BaseMolecularFileSink.h"

using namespace rxdock;

////////////////////////////////////////
// Constructors/destructors
// BaseMolecularFileSink::BaseMolecularFileSink(const char* fileName,
// ModelPtr spModel) :
//  BaseFileSink(fileName), m_spModel(spModel), m_bMultiConf(false)
//{
//  Reset();
//  _RBTOBJECTCOUNTER_CONSTR_("BaseMolecularFileSink");
//}

// DM 21 Apr 1999 - Default is to initialise m_bUseModelSegmentNames to true in
// the constructor so that the first model is rendered using the model segment
// names Default in calls to SetModel is to set m_bUseModelSegmentNames to false
// so that subsequent Render()s use numeric segment IDs
BaseMolecularFileSink::BaseMolecularFileSink(const std::string &fileName,
                                             ModelPtr spModel,
                                             bool bUseModelSegmentNames)
    : BaseFileSink(fileName), m_bUseModelSegmentNames(bUseModelSegmentNames),
      m_spModel(spModel), m_bMultiConf(false) {
  Reset();
  _RBTOBJECTCOUNTER_CONSTR_("BaseMolecularFileSink");
}

// Default destructor
BaseMolecularFileSink::~BaseMolecularFileSink() {
  _RBTOBJECTCOUNTER_DESTR_("BaseMolecularFileSink");
}

////////////////////////////////////////
// Public methods
////////////////

// Get/set the Model that is linked to the sink
ModelPtr BaseMolecularFileSink::GetModel() const { return m_spModel; }

void BaseMolecularFileSink::SetModel(ModelPtr spModel,
                                     bool bUseModelSegmentNames) {
  m_spModel = spModel;
  m_bUseModelSegmentNames = bUseModelSegmentNames;
}

// Get/set the multiconformer status
bool BaseMolecularFileSink::GetMultiConf() const { return m_bMultiConf; }

void BaseMolecularFileSink::SetMultiConf(bool bMultiConf) {
  // If the MultiConf status of the sink is changing, we need to ensure
  // that the cache gets written first
  if (m_bMultiConf != bMultiConf) {
    Write();
    Reset();
  }
  m_bMultiConf = bMultiConf;
}

// Force writing the file following several multi-conf Renders
void BaseMolecularFileSink::WriteMultiConf() {
  Write();
  Reset();
}

ModelList BaseMolecularFileSink::GetSolvent() const { return m_solventList; }

void BaseMolecularFileSink::SetSolvent(ModelList solventList) {
  m_solventList = solventList;
}

////////////////////////////////////////
// Protected methods
///////////////////
// Reset the atom, residue and segment numbering
void BaseMolecularFileSink::Reset() {
  if ((!m_bMultiConf) || isCacheEmpty()) {
    m_nAtomId = 0;
    m_nSubunitId = 0;
    m_nSegmentId = 0;
  }
}

////////////////////////////////////////
// Private methods
/////////////////

void rxdock::to_json(json &j, const BaseMolecularFileSink &baseBiMolFSink) {
  json modelList;
  for (const auto &aIter : baseBiMolFSink.m_solventList) {
    json model = *aIter;
    modelList.push_back(model);
  }

  j = json{{"atom-id", baseBiMolFSink.m_nAtomId},
           {"subunit-id", baseBiMolFSink.m_nSubunitId},
           {"segment-id", baseBiMolFSink.m_nSegmentId},
           {"model-seg-names", baseBiMolFSink.m_bUseModelSegmentNames},

           {"model", *baseBiMolFSink.m_spModel},
           {"solvent-list", modelList},
           {"multi-conf", baseBiMolFSink.m_bMultiConf}};
}

void rxdock::from_json(const json &j, BaseMolecularFileSink &baseBiMolFSink) {
  j.at("atom-id").get_to(baseBiMolFSink.m_nAtomId);
  j.at("subunit-id").get_to(baseBiMolFSink.m_nSubunitId);
  j.at("segment-id").get_to(baseBiMolFSink.m_nSegmentId);
  j.at("model-seg-names").get_to(baseBiMolFSink.m_bUseModelSegmentNames);
  j.at("model").get_to(*baseBiMolFSink.m_spModel);

  for (auto &model : j.at("solvent-list")) {
    ModelPtr spModel = ModelPtr(new Model(model));
    baseBiMolFSink.m_solventList.push_back(spModel);
  }

  j.at("multi-conf").get_to(baseBiMolFSink.m_bMultiConf);
}