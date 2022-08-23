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

#ifndef _RBTMODEL_H_
#define _RBTMODEL_H_

#include "rxdock/Atom.h"
#include "rxdock/Bond.h"
#include "rxdock/Config.h"
#include "rxdock/ModelMutator.h"
#include "rxdock/PrincipalAxes.h"
#include "rxdock/PseudoAtom.h"
#include "rxdock/Variant.h"

#include <nlohmann/json.hpp>

using json = nlohmann::json;

namespace rxdock {

class BaseMolecularFileSource;
class FlexData;
class ChromElement;

class Model {
public:
  //////////////////////
  // Constructors
  //////////////////////
  // Read topology from MolecularFileSource, and coordinates if supported by
  // source
  RBTDLL_EXPORT Model(BaseMolecularFileSource *pMolSource);

  //(Fairly) temporary constructor taking arbitrary atom and bond lists
  // Use with caution
  RBTDLL_EXPORT Model(AtomList &atomList, BondList &bondList);

  Model(json j);

  // Default destructor
  virtual ~Model();

  friend void to_json(json &j, const Model &model);
  friend void from_json(const json &j, Model &model);

  //////////////////////////
  // Friend functions/classes
  //////////////////////////
  // DM 30 Jul 2001 - allow access to atom list by helper class for flexible
  // models
  friend class ModelMutator;

  //////////////////////
  // Public accessor functions
  //////////////////////

  // Only Get functions are currently provided
  // Models are designed to be CreateOnce (from file), ReadMany

  // Model name
  std::string GetName() const { return m_strName; }

  // Titles
  int GetNumTitles() const { return m_titleList.size(); }
  std::vector<std::string> GetTitleList() const { return m_titleList; }
  void SetTitle(std::size_t index, const std::string &strTitle) {
    m_titleList[index] = strTitle;
  }

  // Atoms
  int GetNumAtoms() const { return m_atomList.size(); }
  AtomList GetAtomList() const { return m_atomList; }

  // Bonds
  int GetNumBonds() const { return m_bondList.size(); }
  BondList GetBondList() const { return m_bondList; }

  // Segments
  int GetNumSegments() const { return m_segmentMap.size(); }
  SegmentMap GetSegmentMap() const { return m_segmentMap; }

  // Rings
  int GetNumRings() const { return m_ringList.size(); }
  AtomListList GetRingAtomLists() const { return m_ringList; }

  // DM 12 May 1999 - support for associated model data (e.g. SD file, or
  // generated by rbdock etc) Get number of data fields
  int GetNumData() const { return m_dataMap.size(); }
  // Get list of field names as string list
  std::vector<std::string> GetDataFieldList() const;
  // Get all data as map of key=field name, value=variant (double,string or
  // string list)
  StringVariantMap GetDataMap() const { return m_dataMap; }
  // Query as to whether a particular data field name is present
  bool isDataFieldPresent(const std::string &strDataField) const;
  // Get a particular data value
  RBTDLL_EXPORT Variant GetDataValue(const std::string &strDataField) const;
  // Set a data value (replaces existing value if field name already exists)
  RBTDLL_EXPORT void SetDataValue(const std::string &strDataField,
                                  const Variant &dataValue);
  // Removes a data field completely from the data map
  void ClearDataField(const std::string &strDataField);
  // Removes all data fields starting with a given prefix from the data map
  RBTDLL_EXPORT void ClearAllDataFields(const std::string &strDataFieldPrefix);
  // Removes all data fields from the data map
  void ClearAllDataFields();

  // DM 11 Jul 2000 - pseudoatom handling
  PseudoAtomPtr AddPseudoAtom(const AtomList &atomList);
  void ClearPseudoAtoms();
  void UpdatePseudoAtoms(); // Updates pseudoatom coords
  unsigned int GetNumPseudoAtoms() const;
  PseudoAtomList GetPseudoAtomList() const;

  // DM 1 Jul 2002 - tethered atom handling
  unsigned int GetNumTetheredAtoms() const;
  AtomList GetTetheredAtomList() const;

  // DM 7 June 2006 - occupancy handling
  // The model occupancy is represented in two ways:
  // 1) Occupancy variable (0.0 to 1.0)
  // 2) Enabled state (true/false)
  // If occupancy >= threshold, enabled=true
  // If occupancy < threshold, enabled=false
  // The occupancy variable is managed by the occupancy chromosome element
  // Scoring functions read the enabled state to determined if a model is
  // considered as present or absent
  double GetOccupancy() const { return m_occupancy; }
  bool GetEnabled() const { return m_enabled; }
  // Sets the occupancy and enabled state simultaneously
  RBTDLL_EXPORT void SetOccupancy(double occupancy, double threshold = 0.5);

  //////////////////////
  // Public methods
  //////////////////////

  // Update coords from a data source
  void UpdateCoords(BaseMolecularFileSource *pMolSource);

  // DM 07 Jan 1999
  // Translate molecule by the given vector
  void Translate(const Vector &vector);
  // Rotate molecule around the given axis (through the center of mass) by theta
  // degrees
  void Rotate(const Vector &axis, double thetaDeg);
  // DM 09 Feb 1999 - Rotate molecule around the given axis (through the given
  // coordinate) by theta degrees
  void Rotate(const Vector &axis, double thetaDeg, const Coord &center);
  // Rotate around a given bond by theta degrees (spins both ends of the bond in
  // opposite directions)
  RBTDLL_EXPORT void RotateBond(BondPtr spBond, double thetaDeg);
  // DM 09 Feb 1999 - Rotate around a given bond by theta degrees, keeping the
  // given atom fixed (only spins the other end of the bond)
  void RotateBond(BondPtr spBond, double thetaDeg, AtomPtr spFixedAtom);
  // DM 25 Feb 1999 - Rotate around a given bond by theta degrees, only spinning
  // one end of the bond If bSwap is false, spins the end bonded to atom2 in the
  // bond If bSwap is true, spins the end bonded to atom1 in the bond
  RBTDLL_EXPORT void RotateBond(BondPtr spBond, double thetaDeg, bool bSwap);

  // DM 8 Feb 1999
  void SaveCoords(const std::string &coordName = "");
  RBTDLL_EXPORT void RevertCoords(const std::string &coordName = "");
  std::map<std::string, int> GetSavedCoordNames() const { return m_coordNames; }
  int GetNumSavedCoords() const { return m_coordNames.size(); }
  int GetCurrentCoords() const { return m_currentCoord; }
  RBTDLL_EXPORT void RevertCoords(int);

  // Returns center of mass of model
  Coord GetCenterOfMass() const;
  // Translates the model so that the center of mass is at the given coord
  void SetCenterOfMass(const Coord &c) { Translate(c - GetCenterOfMass()); }

  // DM 9 Nov 1999
  // Returns total atomic mass (molecular weight) for the model
  RBTDLL_EXPORT double GetTotalMass() const;

  // DM 14 Apr 1999 - principal axes methods
  // Return principal axes and center of mass for the model
  PrincipalAxes GetPrincipalAxes() const;
  // Aligns the principal axes of the model to lie along alignAxes
  // If required (bAlignCOM=true), also aligns the center of mass with
  // alignAxes.com Default is to align with X,Y,Z Cartesian axes centered at
  // origin
  void AlignPrincipalAxes(const PrincipalAxes &alignAxes = PrincipalAxes(),
                          bool bAlignCOM = true);

  // DM 19 Oct 2005 - new chromosome handling
  // Sets the flexibility type for this model
  // To revert to a rigid model, call SetFlexData(nullptr)
  RBTDLL_EXPORT void SetFlexData(FlexData *pFlexData);
  // Gets the flexibility type of this model
  // Returns nullptr for a rigid model
  RBTDLL_EXPORT FlexData *GetFlexData() const;
  // Returns a clone of the current chromosome for this model
  // The caller has the responsibility for mem management of the clone
  // Returns nullptr for a rigid model
  RBTDLL_EXPORT ChromElement *GetChrom() const;

  bool isFlexible() const;
  const AtomRList &GetFlexIntns(Atom *pAtom) const;
  BondList GetFlexBonds() const;
  // Selects all atoms that are rotated by at least one rotable bond
  void SelectFlexAtoms();
  void SelectFlexAtoms(Atom *pAtom);

  ////////////////////////////////////////////
  // Atom list functions (provided for convenience, as user could just as well
  // call the rxdock namespace functions with Model::GetAtomList)
  // e.g. AtomList atomList =
  // GetSelectedAtomList(spModel->GetAtomList); is equivalent to
  // AtomList atomList = spModel->GetSelectedAtomList();
  ////////////////////////////////////////////

  // Unary
  // Generic template version of GetNumAtoms, passing in your own predicate
  // template<class Predicate> UInt GetNumAtoms(const Predicate& pred);
  // Generic template version of GetAtomList, passing in your own predicate
  // template<class Predicate> AtomList GetAtomList(const Predicate& pred);

  // Selected atoms
  void SetAtomSelectionFlags(bool bSelected = true);
  unsigned int GetNumSelectedAtoms();
  AtomList GetSelectedAtomList();

  // Cyclic atoms
  void SetAtomCyclicFlags(bool bCyclic = true);
  unsigned int GetNumCyclicAtoms();
  AtomList GetCyclicAtomList();

  // DM 21 Jul 1999 User1 flag
  void SetAtomUser1Flags(bool bUser1 = true);
  // DM 29 Jan 2000 User1 double
  void SetAtomUser1Values(double dUser1 = 0.0);
  // DM 27 Jul 2000 User2 double
  void SetAtomUser2Values(double dUser2 = 0.0);

  // Hydrogen bond acceptor atoms
  unsigned int GetNumHBondAcceptorAtoms();
  AtomList GetHBondAcceptorAtomList();

  // Hydrogen bond donor atoms
  unsigned int GetNumHBondDonorAtoms();
  AtomList GetHBondDonorAtomList();

  //(Formally) charged atoms
  unsigned int GetNumChargedAtoms();
  AtomList GetChargedAtomList();

  // Planar atoms
  unsigned int GetNumPlanarAtoms();
  AtomList GetPlanarAtomList();

  // Binary

  // Atoms with atomic no = nAtomicNo
  unsigned int GetNumAtomsWithAtomicNo(int nAtomicNo);
  AtomList GetAtomListWithAtomicNo(int nAtomicNo);

  // Atoms with FFType = strFFType
  unsigned int GetNumAtomsWithFFType(std::string strFFType);
  AtomList GetAtomListWithFFType(std::string strFFType);

  ////////////////////////////////////////////
  // Bond list functions (provided for convenience, as user could just as well
  // call the rxdock namespace functions with Model::GetBondList)
  // e.g. BondList bondList =
  // GetSelectedBondList(spModel->GetBondList); is equivalent to
  // BondList bondList = spModel->GetSelectedBondList();
  ////////////////////////////////////////////

  // Unary

  // Generic template version of GetNumBonds, passing in your own predicate
  template <class Predicate> unsigned int GetNumBonds(const Predicate &pred);
  // Generic template version of GetBondList, passing in your own predicate
  template <class Predicate> BondList GetBondList(const Predicate &pred);

  // Selected bonds
  void SetBondSelectionFlags(bool bSelected = true);
  unsigned int GetNumSelectedBonds();
  BondList GetSelectedBondList();

  // Cyclic bonds
  void SetBondCyclicFlags(bool bCyclic = true);
  unsigned int GetNumCyclicBonds();
  BondList GetCyclicBondList();

  // DM 10 Dec 1998 - at some point these functions should be implemented as
  // generic rxdock namespace functions operating on arbitrary atom and bond
  // lists

  // D Morley, 2 Dec 1998 - go back to the old way, it's more convenient to get
  // all donors in the same list and separate them later void
  // GetHBondDonorLists(AtomList& donorList, AtomList& donorHList);

  // Get min and max x,y,z coords
  // bInit = true => initialize min, max
  // bInit = false => don't initialize min, max (useful for accumulating min,max
  // over several Models)
  void GetMinMaxCoords(Coord &minCoord, Coord &maxCoord, bool bInit = true);

  // Get map of (key=force field atom type string, value=no. of occurrences)
  std::map<std::string, int> GetAtomTypeMap();

  // Get map of (key=force field bond type (atom type pair) string, value=no. of
  // occurrences)
  std::map<std::string, int> GetBondTypeMap();

protected:
  //////////////////////
  // Private methods
  //////////////////////

private:
  Model();                         // Disable default constructor
  Model(const Model &);            // Copy constructor disabled by default
  Model &operator=(const Model &); // Copy assignment disabled by default

  // Create a new model from a data source
  void Create(BaseMolecularFileSource *pMolSource);
  void Clear();                      // Clear the current model
  void AddAtoms(AtomList &atomList); // Register an atom list with the model

  //////////////////////
  // Private data
  //////////////////////
private:
  std::string m_strName;                // Model name
  std::vector<std::string> m_titleList; // Title list (read from file)
  AtomList m_atomList;                  // atom list
  BondList m_bondList;                  // bond list
  SegmentMap m_segmentMap; // map of (key=segment name, value=atom count)
  AtomListList m_ringList; //(DM 7 Dec 1998) list of atom lists for each ring
  std::map<std::string, int>
      m_coordNames;           //(DM 8 Feb 1999) map of named coord sets
                              //(key=name, value=index into Atom coord map)
  int m_currentCoord;         // DM 11 Jul 2003 - which coord set is current
  StringVariantMap m_dataMap; // DM 12 May 1999 - associated data (e.g. from
                              // SD file or generated by rbdock)
  PseudoAtomList
      m_pseudoAtomList; // DM 11 Jul 2000 - store associated pseudoatoms
  ModelMutatorPtr
      m_spMutator; // DM 30 Jul 2001 - helper object for fast dihedral mutations
  FlexData *m_pFlexData;  // Stores the flexibility data for this model
  ChromElement *m_pChrom; // Reference chromosome. GetChrom() returns a clone
                          // of this object.
  double m_occupancy; // Occupancy value (0->1), in support of solvent occupancy
  bool m_enabled;     // Enabled state, depends on occupancy value and threshold
};

void to_json(json &j, const Model &model);
void from_json(const json &j, Model &model);

// Useful typedefs
typedef SmartPtr<Model> ModelPtr;        // Smart pointer
typedef std::vector<ModelPtr> ModelList; // Vector of smart pointers
typedef ModelList::iterator ModelListIter;
typedef ModelList::const_iterator ModelListConstIter;

} // namespace rxdock

#endif //_RBTMODEL_H_
