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

// Base class for all scoring functions
#ifndef _RBTBASESF_H_
#define _RBTBASESF_H_

#include "rxdock/BaseObject.h"
#include "rxdock/Config.h"

namespace rxdock {

class SFAgg; // forward declaration

class BaseSF : public BaseObject {
public:
  // Class type string
  static std::string _CT;
  // Parameter names
  static std::string _WEIGHT;
  // DM 09 Apr 2002 - max distance range over which scoring function operates
  // Subclasses are free not to use if not required
  static std::string _RANGE;
  // DM 17 Jan 2006 - hardcoded name of scoring function branch under which to
  // save "system" scoring function terms (e.g. intra-receptor, intra-solvent).
  // This is a dirty (and hopefully temporary) solution!
  static std::string _SYSTEM_SF;
  static std::string _INTRA_SF; // Ditto, for ligand intramolecular terms

  ////////////////////////////////////////
  // Constructors/destructors
  virtual ~BaseSF();

  // Give aggregates access to BaseSF private methods/data
  friend class SFAgg;

  ////////////////////////////////////////
  // Public methods
  ////////////////
  // Fully qualified name, prefixed by all ancestors (e.g.
  // rxdock.score.inter.hbond)
  std::string GetFullName() const;
  double GetWeight() const;
  void SetWeight(double);

  double GetRange() const;
  RBTDLL_EXPORT void SetRange(double);

  // Main public method - returns current weighted score
  RBTDLL_EXPORT double Score() const;
  // Returns all child component scores as a string-variant map
  // Key = fully qualified component name, value = weighted score
  //(for saving in a Model's data fields)
  virtual void ScoreMap(StringVariantMap &scoreMap) const;

  // Aggregate handling methods
  virtual void Add(BaseSF *);
  virtual void Remove(BaseSF *);
  virtual bool isAgg() const;
  virtual unsigned int GetNumSF() const;
  virtual BaseSF *GetSF(unsigned int iSF) const;
  void Orphan(); // Force removal from the parent aggregate
  BaseSF *GetParentSF() const;

protected:
  ////////////////////////////////////////
  // Protected methods
  ///////////////////
  BaseSF(const std::string &strClass, const std::string &strName);
  BaseSF();
  // PURE VIRTUAL - DERIVED CLASSES MUST OVERRIDE
  virtual double RawScore() const = 0;
  // DM 25 Oct 2000 - track changes to parameter values in local data members
  // ParameterUpdated is invoked by ParamHandler::SetParameter
  void ParameterUpdated(const std::string &strName);
  // Helper method for ScoreMap
  void AddToParentMapEntry(StringVariantMap &scoreMap, double rs) const;

private:
  ////////////////////////////////////////
  // Private methods
  /////////////////
  BaseSF(const BaseSF &);            // Copy constructor disabled by default
  BaseSF &operator=(const BaseSF &); // Copy assignment disabled by
                                     // default

protected:
  ////////////////////////////////////////
  // Protected data
  ////////////////

private:
  ////////////////////////////////////////
  // Private data
  //////////////
  BaseSF *m_parent;
  double m_weight;
  double m_range;
};

// Useful typedefs
typedef std::vector<BaseSF *> BaseSFList; // Vector of smart pointers
typedef BaseSFList::iterator BaseSFListIter;
typedef BaseSFList::const_iterator BaseSFListConstIter;

} // namespace rxdock

#endif //_RBTBASESF_H_
