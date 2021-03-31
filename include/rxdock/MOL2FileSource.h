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

#ifndef _RBTMOL2FILESOURCE_H_
#define _RBTMOL2FILESOURCE_H_

#include "rxdock/BaseMolecularFileSource.h"
#include "rxdock/ElementFileSource.h"
#include "rxdock/ParameterFileSource.h"

namespace rxdock {

// Simple struct for holding MOL2 SUBSTRUCTURE records
class MOL2Substructure {
public:
  MOL2Substructure() : root_atom(0), chain("UNK"), sub_type("UNK") {}
  MOL2Substructure(const std::string &name, int rr, const std::string &cc,
                   const std::string &ss)
      : subst_name(name), root_atom(rr), chain(cc), sub_type(ss) {}
  std::string GetName() const { return subst_name; };
  int GetRoot() const { return root_atom; };
  std::string GetChain() const { return chain; };
  std::string GetType() const { return sub_type; };

private:
  std::string subst_name;
  int root_atom;
  std::string chain;
  std::string sub_type;
};

// A MOL2SubstructureMap stores MOL2Substructures keyed by SUBSTRUCTURE ID
// (subst_id)
typedef std::map<int, MOL2Substructure> MOL2SubstructureMap;
typedef MOL2SubstructureMap::iterator MOL2SubstructureMapIter;
typedef MOL2SubstructureMap::const_iterator MOL2SubstructureMapConstIter;

class MOL2FileSource : public BaseMolecularFileSource {
public:
  static const std::string _CT;
  // record delimiter strings
  static const std::string _TRIPOS_DELIM;
  static const unsigned int _TRIPOS_DELIM_SIZE;
  static const std::string _IDS_MOL2_RECDELIM;

  MOL2FileSource(const std::string &fileName, bool bImplHydrogens = true);
  ~MOL2FileSource();
  bool isTitleListSupported() { return true; }
  bool isAtomListSupported() { return true; }
  bool isCoordinatesSupported() { return true; }
  bool isBondListSupported() { return true; }
  bool isDataSupported() { return false; }

protected:
  virtual void Parse();

private:
  MOL2FileSource();
  MOL2FileSource(const MOL2FileSource &);
  MOL2FileSource &operator=(const MOL2FileSource &);

  // Helper functions to clean up atom attributes after parsing the file
  void SetupAtomParams();
  void FixImplicitHydrogenCount();
  void FixHybridState();
  void FixTriposTypes();
  void SetupVdWRadii();
  void RemoveNonPolarHydrogens();

  // Splits subst_name into ID and name components (e.g. ALA27 => ALA and 27)
  void GetSSIDandName(const std::string &subst_name, int subst_id,
                      std::string &sID, std::string &sName);

  // parser functions (switched by pointers)
  void ParseRecordMOLECULE(const std::string &aLine);
  void ParseRecordATOM(const std::string &aLine);
  void ParseRecordBOND(const std::string &aLine);
  void ParseRecordSUBSTRUCTURE(const std::string &aLine);
  void ParseRecordUNSUPPORTED(const std::string &aLine);
  // get delimiter tag
  std::string GetMOL2Tag(const std::string &aLine);
  // parse/tokenize second line (number of  atoms, bonds etc) in MOLECULE
  void ParseCountFields(const std::string &aLine);
  void Tokenize(const std::string &aString,
                std::vector<std::string> &aTokensBuf);

  unsigned int m_NL; // No. of lines
  // data fields in MOLECULE
  unsigned int nAtoms;         // number of atoms
  unsigned int nBonds;         //           bonds
  unsigned int nSubstructures; // substructures
  unsigned int nFeatures;      // features
  unsigned int nSets;          // sets

  // helper data structures used to lookup chain names for each substructure
  // and to sort atom list into a more sensible order (ordered by substructure
  // ID)
  MOL2SubstructureMap m_ssInfo;
  std::map<int, AtomList> m_ssAtoms;

  ElementFileSourcePtr
      m_spElementData; // Elemental data source for vdW radii etc
  ParameterFileSourcePtr
      m_spParamSource;    // Parameter file source for partial charges etc
  TriposAtomType m_typer; // Tripos atom typing object
  bool m_bImplHydrogens;
};

// usual smart pointer
typedef SmartPtr<MOL2FileSource> MOL2FileSourcePtr;
// function pointer (callback) for Parse()
typedef void (MOL2FileSource::*fcPtr)(const std::string &);

} // namespace rxdock

#endif // _RBTMOL2FILESOURCE_H_
