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
// Rbt::Observer objects (scoring functions and transforms)
// In this version, the number of models in the workspace is fixed in the
// constructor, and can't be changed. Models can still be null however.
//
// Design Pattern: Workspace objects are the Subject in the Observer pattern
//                (Design Patterns, Gamma et al, Addison Wesley, p293).

#ifndef _RBTWORKSPACE_H_
#define _RBTWORKSPACE_H_

#include "RbtBaseMolecularFileSink.h"
#include "RbtConfig.h"
#include "RbtDockingSite.h"
#include "RbtFilter.h"
#include "RbtModel.h"
#include "RbtParamHandler.h"
#include "RbtPopulation.h"
#include "RbtSubject.h"

class RbtBaseSF;        // Forward definition
class RbtBaseTransform; // Forward definition

class RbtWorkSpace : public RbtSubject, public RbtParamHandler {
public:
  // Class type string
  static std::string _CT;
  // Parameter names
  static std::string _NAME;

  ////////////////////////////////////////
  // Constructors/destructors

  RbtWorkSpace(unsigned int nModels = 2);
  virtual ~RbtWorkSpace();

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
  RBTDLL_EXPORT RbtModelList GetModels() const;
  // Returns a specific (numbered) model
  // Throws RbtBadArgument if iModel out of range
  RBTDLL_EXPORT RbtModelPtr GetModel(unsigned int iModel) const;

  // Replace an existing model
  // Throws RbtBadArgument if iModel out of range
  void SetModel(unsigned int iModel, RbtModelPtr spModel);

  // Returns vector of models, starting from index iModel
  RbtModelList GetModels(unsigned int iModel) const;
  // Append a number of models to the workspace, increasing the total number of
  // models
  void AddModels(RbtModelList modelList);
  // Replace a number of existing models
  // iModel is the index of the first model to replace
  // Throws RbtBadArgument if iModel out of range or modelList too large
  void SetModels(unsigned int iModel, RbtModelList modelList);
  // Removes a number of models from the workspace
  // Removes from index iModel to end of model list
  void RemoveModels(unsigned int iModel);

  // Model I/O
  // Get/set the molecular file sink (for outputting ligands)
  RbtMolecularFileSinkPtr GetSink() const;
  RBTDLL_EXPORT void SetSink(RbtMolecularFileSinkPtr);
  // Saves models to file sink
  // If bSaveScores is true, write component scores also
  // Base workspace does nothing
  // It is up to subclasses to decide what action to take
  virtual void Save(bool bSaveScores = true);
  // Get/set the history file sink
  RbtMolecularFileSinkPtr GetHistorySink() const;
  RBTDLL_EXPORT void SetHistorySink(RbtMolecularFileSinkPtr);
  // Saves models to file sink
  // If bSaveScores is true, write component scores also
  // Base workspace does nothing
  // It is up to subclasses to decide what action to take
  virtual void SaveHistory(bool bSaveScores = true);

  // Get/Set the scoring function
  // SetSF automatically registers the scoring function with the workspace
  RBTDLL_EXPORT RbtBaseSF *GetSF() const;
  RBTDLL_EXPORT void SetSF(RbtBaseSF *);

  // Get/Set the transform
  // SetTransform automatically registers the transform with the workspace
  RbtBaseTransform *GetTransform() const;
  RBTDLL_EXPORT void SetTransform(RbtBaseTransform *);

  // Run the simulation!
  virtual void Run();

  // Population handling
  // DM 30 Oct 2001 - provide a home for persistent populations
  void SetPopulation(RbtPopulationPtr population);
  RbtPopulationPtr GetPopulation() const;
  void ClearPopulation();

  // Docking site handling
  // DM 09 Apr 2002 - workspace now manages the docking site
  RBTDLL_EXPORT RbtDockingSitePtr GetDockingSite() const;
  RBTDLL_EXPORT void SetDockingSite(RbtDockingSitePtr spDockSite);

  // Filter handling
  // BGD 27 Jan 2003 - SetFilter automotically registers the Filter
  // with the workspace
  RbtFilterPtr GetFilter() const;
  RBTDLL_EXPORT void SetFilter(RbtFilterPtr spFilter);

protected:
  ////////////////////////////////////////
  // Protected methods
  ///////////////////

private:
  ////////////////////////////////////////
  // Private methods
  /////////////////

  RbtWorkSpace(const RbtWorkSpace &); // Copy constructor disabled by default

  RbtWorkSpace &
  operator=(const RbtWorkSpace &); // Copy assignment disabled by default

protected:
  ////////////////////////////////////////
  // Protected data
  ////////////////

private:
  ////////////////////////////////////////
  // Private data
  //////////////
  RbtModelList m_models;
  RbtMolecularFileSinkPtr m_spSink;
  RbtMolecularFileSinkPtr m_spHisSink;
  RbtBaseSF *m_SF;
  RbtBaseTransform *m_transform;
  RbtPopulationPtr m_population;
  RbtDockingSitePtr m_spDockSite;
  RbtFilterPtr m_spFilter;
};

// Useful typedefs
typedef SmartPtr<RbtWorkSpace> RbtWorkSpacePtr; // Smart pointer

#endif //_RBTWORKSPACE_H_
