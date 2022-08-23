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

// Abstract base class for molecular file-based data sources

#ifndef _RBTBASEMOLECULARFILESOURCE_H_
#define _RBTBASEMOLECULARFILESOURCE_H_

#include "rxdock/BaseFileSource.h"
#include "rxdock/BaseMolecularDataSource.h"
#include "rxdock/ParameterFileSource.h"

#include <nlohmann/json.hpp>

using json = nlohmann::json;

namespace rxdock {

class BaseMolecularFileSource : public BaseFileSource,
                                public BaseMolecularDataSource {
public:
  // Constructors
  // BaseMolecularFileSource(const char* fileName, const char* sourceName);
  // Single-record version
  RBTDLL_EXPORT BaseMolecularFileSource(const std::string &fileName,
                                        const std::string &sourceName);
  // Multi-record version, with record delimiter passed as an argument
  BaseMolecularFileSource(const std::string &fileName,
                          const std::string &strRecDelim,
                          const std::string &sourceName);

  // Destructor
  virtual ~BaseMolecularFileSource();

  friend void to_json(json &j, const BaseMolecularFileSource &baseMolFilSrc);
  friend void from_json(const json &j, BaseMolecularFileSource &baseMolFilSrc);

  // Pure virtual methods from BaseMolecularDataSource
  void Reset();
  int GetNumTitles();
  int GetNumAtoms();
  int GetNumBonds();
  int GetNumSegments();

  std::vector<std::string> GetTitleList();
  AtomList GetAtomList();
  BondList GetBondList();
  SegmentMap GetSegmentMap();

  // DM 12 May 1999 - support for data records (e.g. SD file)
  // Get number of data fields
  virtual int GetNumData();
  // Get list of field names as string list
  virtual std::vector<std::string> GetDataFieldList();
  // Get all data as map of key=field name, value=variant (double,string or
  // string list)
  virtual StringVariantMap GetDataMap();
  // Query as to whether a particular data field name is present
  virtual bool isDataFieldPresent(const std::string &strDataField);
  // Get a particular data value
  virtual Variant GetDataValue(const std::string &strDataField);

  // These methods allow filtering of the data source by segment name
  // So for example, we could just read the segment named TAR from the source
  // The filter modifies the behaviour of GetAtomList, GetBondList etc
  SegmentMap GetSegmentFilterMap();
  RBTDLL_EXPORT void SetSegmentFilterMap(const SegmentMap &segmentFilterMap);
  void ClearSegmentFilterMap();
  bool isSegmentFilterMapDefined();

protected:
  // Protected methods
  void ClearMolCache();            // Clear the atom, bond and title lists
  void RemoveAtom(AtomPtr spAtom); // DM 16 Feb 2000 - moved from MdlFileSource
  void
  RenumberAtomsAndBonds(); // DM 30 Oct 2000 - now independent from RemoveAtom
  // Moved from MOL2FileSource - sets partial charges on a per-residue basis
  // Automatically determines protonation states of ASP/GLU/HIS using Mandatory
  // and Forbidden atoms in spParamSource
  void SetupPartialIonicGroups(AtomList &atoms,
                               ParameterFileSourcePtr spParamSource);

private:
  // Private methods
  BaseMolecularFileSource(); // Disable default constructor
  BaseMolecularFileSource(
      const BaseMolecularFileSource &); // Copy constructor disabled by default
  BaseMolecularFileSource &operator=(
      const BaseMolecularFileSource &); // Copy assignment disabled by default

  // Get methods with segment filter applied
  int GetNumAtomsWithFilter();
  int GetNumBondsWithFilter();
  AtomList GetAtomListWithFilter();
  BondList GetBondListWithFilter();

protected:
  // Protected data
  std::vector<std::string> m_titleList;
  AtomList m_atomList;
  BondList m_bondList;
  SegmentMap m_segmentMap;
  StringVariantMap m_dataMap; // DM 12 May 1999 - for storing data records

private:
  // Private data
  SegmentMap m_segmentFilterMap;
};

void to_json(json &j, const BaseMolecularFileSource &baseMolFilSrc);
void from_json(const json &j, BaseMolecularFileSource &baseMolFilSrc);

// useful typedefs
typedef SmartPtr<BaseMolecularFileSource>
    MolecularFileSourcePtr; // Smart pointer

} // namespace rxdock

#endif //_RBTBASEMOLECULARFILESOURCE_H_
