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

// Manages all models used in simulation, coordinates updates of all registered
// Observer objects (scoring functions and transforms)
// In this version, the number of models in the workspace is fixed in the
// constructor, and can't be changed. Models can still be null however.
//
// Design Pattern: Workspace objects are the Subject in the Observer pattern
//                (Design Patterns, Gamma et al, Addison Wesley, p293).

#ifndef _RBTWORKSPACE_H_
#define _RBTWORKSPACE_H_

#include "rxdock/BaseMolecularFileSink.h"
#include "rxdock/Config.h"
#include "rxdock/DockingSite.h"
#include "rxdock/Filter.h"
#include "rxdock/Model.h"
#include "rxdock/ParamHandler.h"
#include "rxdock/Population.h"
#include "rxdock/Subject.h"

namespace rxdock {

class BaseSF;        // Forward definition
class BaseTransform; // Forward definition

class WorkSpace : public Subject, public ParamHandler {
public:
  // Class type string
  static std::string _CT;
  // Parameter names
  static std::string _NAME;

  ////////////////////////////////////////
  // Constructors/destructors

  WorkSpace(unsigned int nModels = 2);
  virtual ~WorkSpace();

  ////////////////////////////////////////
  // Public methods
  ////////////////

  // Get/set workspace name
  RBTDLL_EXPORT std::string GetName() const;
  RBTDLL_EXPORT void SetName(const std::string &);

  // Model handling

  // Returns number of models in workspace
  RBTDLL_EXPORT unsigned int GetNumModels() const;
  // Returns vector of all models in workspace
  RBTDLL_EXPORT ModelList GetModels() const;
  // Returns a specific (numbered) model
  // Throws BadArgument if iModel out of range
  RBTDLL_EXPORT ModelPtr GetModel(unsigned int iModel) const;

  // Replace an existing model
  // Throws BadArgument if iModel out of range
  void SetModel(unsigned int iModel, ModelPtr spModel);

  // Returns vector of models, starting from index iModel
  ModelList GetModels(unsigned int iModel) const;
  // Append a number of models to the workspace, increasing the total number of
  // models
  void AddModels(ModelList modelList);
  // Replace a number of existing models
  // iModel is the index of the first model to replace
  // Throws BadArgument if iModel out of range or modelList too large
  void SetModels(unsigned int iModel, ModelList modelList);
  // Removes a number of models from the workspace
  // Removes from index iModel to end of model list
  void RemoveModels(unsigned int iModel);

  // Model I/O
  // Get/set the molecular file sink (for outputting ligands)
  MolecularFileSinkPtr GetSink() const;
  RBTDLL_EXPORT void SetSink(MolecularFileSinkPtr);
  // Saves models to file sink
  // If bSaveScores is true, write component scores also
  // Base workspace does nothing
  // It is up to subclasses to decide what action to take
  virtual void Save(bool bSaveScores = true);
  // Get/set the history file sink
  MolecularFileSinkPtr GetHistorySink() const;
  RBTDLL_EXPORT void SetHistorySink(MolecularFileSinkPtr);
  // Saves models to file sink
  // If bSaveScores is true, write component scores also
  // Base workspace does nothing
  // It is up to subclasses to decide what action to take
  virtual void SaveHistory(bool bSaveScores = true);

  // Get/Set the scoring function
  // SetSF automatically registers the scoring function with the workspace
  RBTDLL_EXPORT BaseSF *GetSF() const;
  RBTDLL_EXPORT void SetSF(BaseSF *);

  // Get/Set the transform
  // SetTransform automatically registers the transform with the workspace
  BaseTransform *GetTransform() const;
  RBTDLL_EXPORT void SetTransform(BaseTransform *);

  // Run the simulation!
  virtual void Run();

  // Population handling
  // DM 30 Oct 2001 - provide a home for persistent populations
  void SetPopulation(PopulationPtr population);
  PopulationPtr GetPopulation() const;
  void ClearPopulation();

  // Docking site handling
  // DM 09 Apr 2002 - workspace now manages the docking site
  RBTDLL_EXPORT DockingSitePtr GetDockingSite() const;
  RBTDLL_EXPORT void SetDockingSite(DockingSitePtr spDockSite);

  // Filter handling
  // BGD 27 Jan 2003 - SetFilter automotically registers the Filter
  // with the workspace
  FilterPtr GetFilter() const;
  RBTDLL_EXPORT void SetFilter(FilterPtr spFilter);

protected:
  ////////////////////////////////////////
  // Protected methods
  ///////////////////

private:
  ////////////////////////////////////////
  // Private methods
  /////////////////

  WorkSpace(const WorkSpace &); // Copy constructor disabled by default

  WorkSpace &
  operator=(const WorkSpace &); // Copy assignment disabled by default

protected:
  ////////////////////////////////////////
  // Protected data
  ////////////////

private:
  ////////////////////////////////////////
  // Private data
  //////////////
  ModelList m_models;
  MolecularFileSinkPtr m_spSink;
  MolecularFileSinkPtr m_spHisSink;
  BaseSF *m_SF;
  BaseTransform *m_transform;
  PopulationPtr m_population;
  DockingSitePtr m_spDockSite;
  FilterPtr m_spFilter;
};

// Useful typedefs
typedef SmartPtr<WorkSpace> WorkSpacePtr; // Smart pointer

} // namespace rxdock

#endif //_RBTWORKSPACE_H_
