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

#include "WorkSpace.h"
#include "BaseSF.h"
#include "BaseTransform.h"

using namespace rxdock;

// Static data members
std::string WorkSpace::_CT("WorkSpace");
std::string WorkSpace::_NAME("NAME");

////////////////////////////////////////
// Constructors/destructors

// Create an empty model container of the right size
WorkSpace::WorkSpace(unsigned int nModels)
    : m_models(nModels), m_SF(nullptr), m_transform(nullptr) {
  AddParameter(_NAME, _CT);
#ifdef _DEBUG
  std::cout << "WorkSpace::WorkSpace(): Created model list of size "
            << m_models.size() << std::endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_CONSTR_("WorkSpace");
}

WorkSpace::~WorkSpace() { _RBTOBJECTCOUNTER_DESTR_("WorkSpace"); }

////////////////////////////////////////
// Public methods
////////////////

// Get/set workspace name
std::string WorkSpace::GetName() const { return GetParameter(_NAME); }
void WorkSpace::SetName(const std::string &strName) {
  SetParameter(_NAME, strName);
}

// Model handling

// Returns number of models in workspace
unsigned int WorkSpace::GetNumModels() const { return m_models.size(); }

// Returns vector of all models in workspace
ModelList WorkSpace::GetModels() const { return m_models; }

// Returns a specific (numbered) model
// Throws BadArgument if iModel out of range
ModelPtr WorkSpace::GetModel(unsigned int iModel) const {
  if (iModel >= m_models.size()) {
    throw BadArgument(_WHERE_, "iModel out of range");
  }
  return m_models[iModel];
}

// Replace an existing model, notify all observers
// Throws BadArgument if iModel out of range
void WorkSpace::SetModel(unsigned int iModel, ModelPtr spModel) {
  if (iModel >= m_models.size()) {
    throw BadArgument(_WHERE_, "iModel out of range");
  }
  m_models[iModel] = spModel;
  Notify(); // Notify observers of change in state
}

// Returns vector of models, starting from index iModel
ModelList WorkSpace::GetModels(unsigned int iModel) const {
  ModelList retVal;
  if (iModel >= m_models.size()) {
    throw BadArgument(_WHERE_, "iModel out of range");
  }
  std::copy(m_models.begin() + iModel, m_models.end(),
            std::back_inserter(retVal));
  return retVal;
}

// Append a number of models to the workspace, increasing the total number of
// models
void WorkSpace::AddModels(ModelList modelList) {
  std::copy(modelList.begin(), modelList.end(), std::back_inserter(m_models));
  Notify(); // Notify observers of change in state
}

// Replace a number of existing models
// iModel is the index of the first model to replace
// Throws BadArgument if iModel out of range or modelList too large
void WorkSpace::SetModels(unsigned int iModel, ModelList modelList) {
  if (iModel > (m_models.size() - modelList.size())) {
    throw BadArgument(_WHERE_, "iModel out of range");
  } else if (modelList.size() > m_models.size()) {
    throw BadArgument(_WHERE_, "modelList too large");
  }
  std::copy(modelList.begin(), modelList.end(), m_models.begin() + iModel);
  Notify(); // Notify observers of change in state
}

// Removes a number of models from the workspace
// Removes from index iModel to end of model list
void WorkSpace::RemoveModels(unsigned int iModel) {
  if (iModel >= m_models.size()) {
    throw BadArgument(_WHERE_, "iModel out of range");
  }
  m_models.erase(m_models.begin() + iModel, m_models.end());
  Notify();
}

// Model I/O
// Get/set the molecular file sink (for outputting ligands)
MolecularFileSinkPtr WorkSpace::GetSink() const { return m_spSink; }
void WorkSpace::SetSink(MolecularFileSinkPtr spSink) { m_spSink = spSink; }

// Saves models to file sink
// Base workspace does nothing
// It is up to subclasses to decide what action to take
void WorkSpace::Save(bool bSaveScores) {}

// Get/set the history file sink
MolecularFileSinkPtr WorkSpace::GetHistorySink() const { return m_spHisSink; }
void WorkSpace::SetHistorySink(MolecularFileSinkPtr spHisSink) {
  m_spHisSink = spHisSink;
}
// Saves models to file sink
// Base workspace does nothing
// It is up to subclasses to decide what action to take
void WorkSpace::SaveHistory(bool bSaveScores) {}

// Scoring function handling
// SetSF automatically registers the scoring function with the workspace
BaseSF *WorkSpace::GetSF() const { return m_SF; }
void WorkSpace::SetSF(BaseSF *pSF) {
  if (m_SF) {
    m_SF->Unregister();
  }
  m_SF = pSF;
  if (m_SF) {
    m_SF->Register(this);
    Notify();
  }
}

// Get/Set the transform
// SetTransform automatically registers the transform with the workspace
BaseTransform *WorkSpace::GetTransform() const { return m_transform; }
void WorkSpace::SetTransform(BaseTransform *pTransform) {
  if (m_transform) {
    m_transform->Unregister();
  }
  m_transform = pTransform;
  if (m_transform) {
    m_transform->Register(this);
    Notify();
  }
}

// Run the simulation!
void WorkSpace::Run() {
  if (m_transform) {
    m_transform->Go();
  }
}

// Population handling
// DM 30 Oct 2001 - provide a home for persistent populations
void WorkSpace::SetPopulation(PopulationPtr population) {
  m_population = population;
}

PopulationPtr WorkSpace::GetPopulation() const { return m_population; }

void WorkSpace::ClearPopulation() { m_population.SetNull(); }

// Docking site handling
// DM 09 Apr 2002 - workspace now manages the docking site
DockingSitePtr WorkSpace::GetDockingSite() const { return m_spDockSite; }

void WorkSpace::SetDockingSite(DockingSitePtr spDockSite) {
  m_spDockSite = spDockSite;
}

// Filter handling
// BGD 27 Jan 2003 - SetFilter automotically registers the Filter
// with the workspace
FilterPtr WorkSpace::GetFilter() const { return m_spFilter; }

void WorkSpace::SetFilter(FilterPtr spFilter) {
  if (m_spFilter) {
    m_spFilter->Unregister();
  }
  m_spFilter = spFilter;
  if (m_spFilter) {
    m_spFilter->Register(this);
  }
}
