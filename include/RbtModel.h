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

#include "RbtAtom.h"
#include "RbtBond.h"
#include "RbtConfig.h"
#include "RbtModelMutator.h"
#include "RbtPrincipalAxes.h"
#include "RbtPseudoAtom.h"
#include "RbtVariant.h"

class RbtBaseMolecularFileSource;
class RbtFlexData;
class RbtChromElement;

class RbtModel {
public:
  //////////////////////
  // Constructors
  //////////////////////
  // Read topology from RbtMolecularFileSource, and coordinates if supported by
  // source
  RbtModel(RbtBaseMolecularFileSource *pMolSource);

  //(Fairly) temporary constructor taking arbitrary atom and bond lists
  // Use with caution
  RbtModel(RbtAtomList &atomList, RbtBondList &bondList);

  // Default destructor
  virtual ~RbtModel();

  //////////////////////////
  // Friend functions/classes
  //////////////////////////
  // DM 30 Jul 2001 - allow access to atom list by helper class for flexible
  // models
  friend class RbtModelMutator;

  //////////////////////
  // Public accessor functions
  //////////////////////

  // Only Get functions are currently provided
  // RbtModels are designed to be CreateOnce (from file), ReadMany

  // Model name
  std::string GetName() const { return m_strName; }

  // Titles
  RbtInt GetNumTitles() const { return m_titleList.size(); }
  RbtStringList GetTitleList() const { return m_titleList; }

  // Atoms
  RbtInt GetNumAtoms() const { return m_atomList.size(); }
  RbtAtomList GetAtomList() const { return m_atomList; }

  // Bonds
  RbtInt GetNumBonds() const { return m_bondList.size(); }
  RbtBondList GetBondList() const { return m_bondList; }

  // Segments
  RbtInt GetNumSegments() const { return m_segmentMap.size(); }
  RbtSegmentMap GetSegmentMap() const { return m_segmentMap; }

  // Rings
  RbtInt GetNumRings() const { return m_ringList.size(); }
  RbtAtomListList GetRingAtomLists() const { return m_ringList; }

  // DM 12 May 1999 - support for associated model data (e.g. SD file, or
  // generated by rbdock etc) Get number of data fields
  RbtInt GetNumData() const { return m_dataMap.size(); }
  // Get list of field names as string list
  RbtStringList GetDataFieldList() const;
  // Get all data as map of key=field name, value=variant (double,string or
  // string list)
  RbtStringVariantMap GetDataMap() const { return m_dataMap; }
  // Query as to whether a particular data field name is present
  RbtBool isDataFieldPresent(const std::string &strDataField) const;
  // Get a particular data value
  RbtVariant GetDataValue(const std::string &strDataField) const;
  // Set a data value (replaces existing value if field name already exists)
  void SetDataValue(const std::string &strDataField,
                    const RbtVariant &dataValue);
  // Removes a data field completely from the data map
  void ClearDataField(const std::string &strDataField);
  // Removes all data fields starting with a given prefix from the data map
  void ClearAllDataFields(const std::string &strDataFieldPrefix);
  // Removes all data fields from the data map
  void ClearAllDataFields();

  // DM 11 Jul 2000 - pseudoatom handling
  RbtPseudoAtomPtr AddPseudoAtom(const RbtAtomList &atomList);
  void ClearPseudoAtoms();
  void UpdatePseudoAtoms(); // Updates pseudoatom coords
  RbtUInt GetNumPseudoAtoms() const;
  RbtPseudoAtomList GetPseudoAtomList() const;

  // DM 1 Jul 2002 - tethered atom handling
  RbtUInt GetNumTetheredAtoms() const throw(RbtError);
  RbtAtomList GetTetheredAtomList() const throw(RbtError);

  // DM 7 June 2006 - occupancy handling
  // The model occupancy is represented in two ways:
  // 1) Occupancy variable (0.0 to 1.0)
  // 2) Enabled state (true/false)
  // If occupancy >= threshold, enabled=true
  // If occupancy < threshold, enabled=false
  // The occupancy variable is managed by the occupancy chromosome element
  // Scoring functions read the enabled state to determined if a model is
  // considered as present or absent
  RbtDouble GetOccupancy() const { return m_occupancy; }
  RbtBool GetEnabled() const { return m_enabled; }
  // Sets the occupancy and enabled state simultaneously
  void SetOccupancy(RbtDouble occupancy, RbtDouble threshold = 0.5);

  //////////////////////
  // Public methods
  //////////////////////

  // Update coords from a data source
  void UpdateCoords(RbtBaseMolecularFileSource *pMolSource) throw(RbtError);

  // DM 07 Jan 1999
  // Translate molecule by the given vector
  void Translate(const RbtVector &vector);
  // Rotate molecule around the given axis (through the center of mass) by theta
  // degrees
  void Rotate(const RbtVector &axis, RbtDouble thetaDeg);
  // DM 09 Feb 1999 - Rotate molecule around the given axis (through the given
  // coordinate) by theta degrees
  void Rotate(const RbtVector &axis, RbtDouble thetaDeg,
              const RbtCoord &center);
  // Rotate around a given bond by theta degrees (spins both ends of the bond in
  // opposite directions)
  void RotateBond(RbtBondPtr spBond, RbtDouble thetaDeg);
  // DM 09 Feb 1999 - Rotate around a given bond by theta degrees, keeping the
  // given atom fixed (only spins the other end of the bond)
  void RotateBond(RbtBondPtr spBond, RbtDouble thetaDeg,
                  RbtAtomPtr spFixedAtom);
  // DM 25 Feb 1999 - Rotate around a given bond by theta degrees, only spinning
  // one end of the bond If bSwap is false, spins the end bonded to atom2 in the
  // bond If bSwap is true, spins the end bonded to atom1 in the bond
  void RotateBond(RbtBondPtr spBond, RbtDouble thetaDeg, RbtBool bSwap);

  // DM 8 Feb 1999
  void SaveCoords(const std::string &coordName = "");
  void RevertCoords(const std::string &coordName = "") throw(RbtError);
  RbtStringIntMap GetSavedCoordNames() const { return m_coordNames; }
  RbtInt GetNumSavedCoords() const { return m_coordNames.size(); }
  RbtInt GetCurrentCoords() const { return m_currentCoord; }
  void RevertCoords(RbtInt);

  // Returns center of mass of model
  RbtCoord GetCenterOfMass() const;
  // Translates the model so that the center of mass is at the given coord
  void SetCenterOfMass(const RbtCoord &c) { Translate(c - GetCenterOfMass()); }

  // DM 9 Nov 1999
  // Returns total atomic mass (molecular weight) for the model
  RbtDouble GetTotalAtomicMass() const;

  // DM 14 Apr 1999 - principal axes methods
  // Return principal axes and center of mass for the model
  RbtPrincipalAxes GetPrincipalAxes() const;
  // Aligns the principal axes of the model to lie along alignAxes
  // If required (bAlignCOM=true), also aligns the center of mass with
  // alignAxes.com Default is to align with X,Y,Z Cartesian axes centered at
  // origin
  void
  AlignPrincipalAxes(const RbtPrincipalAxes &alignAxes = RbtPrincipalAxes(),
                     RbtBool bAlignCOM = true);

  // DM 19 Oct 2005 - new chromosome handling
  // Sets the flexibility type for this model
  // To revert to a rigid model, call SetFlexData(NULL)
  void SetFlexData(RbtFlexData *pFlexData);
  // Gets the flexibility type of this model
  // Returns NULL for a rigid model
  RbtFlexData *GetFlexData() const;
  // Returns a clone of the current chromosome for this model
  // The caller has the responsibility for mem management of the clone
  // Returns NULL for a rigid model
  RbtChromElement *GetChrom() const;

  RbtBool isFlexible() const;
  const RbtAtomRList &GetFlexIntns(RbtAtom *pAtom) const throw(RbtError);
  RbtBondList GetFlexBonds() const throw(RbtError);
  // Selects all atoms that are rotated by at least one rotable bond
  void SelectFlexAtoms();
  void SelectFlexAtoms(RbtAtom *pAtom);

  ////////////////////////////////////////////
  // Atom list functions (provided for convenience, as user could just as well
  // call the Rbt:: functions with RbtModel::GetAtomList)
  // e.g. RbtAtomList atomList = Rbt::GetSelectedAtomList(spModel->GetAtomList);
  // is equivalent to
  // RbtAtomList atomList = spModel->GetSelectedAtomList();
  ////////////////////////////////////////////

  // Unary
  // Generic template version of GetNumAtoms, passing in your own predicate
  // template<class Predicate> RbtUInt GetNumAtoms(const Predicate& pred);
  // Generic template version of GetAtomList, passing in your own predicate
  // template<class Predicate> RbtAtomList GetAtomList(const Predicate& pred);

  // Selected atoms
  void SetAtomSelectionFlags(RbtBool bSelected = true);
  RbtUInt GetNumSelectedAtoms();
  RbtAtomList GetSelectedAtomList();

  // Cyclic atoms
  void SetAtomCyclicFlags(RbtBool bCyclic = true);
  RbtUInt GetNumCyclicAtoms();
  RbtAtomList GetCyclicAtomList();

  // DM 21 Jul 1999 User1 flag
  void SetAtomUser1Flags(RbtBool bUser1 = true);
  // DM 29 Jan 2000 User1 double
  void SetAtomUser1Values(RbtDouble dUser1 = 0.0);
  // DM 27 Jul 2000 User2 double
  void SetAtomUser2Values(RbtDouble dUser2 = 0.0);

  // Hydrogen bond acceptor atoms
  RbtUInt GetNumHBondAcceptorAtoms();
  RbtAtomList GetHBondAcceptorAtomList();

  // Hydrogen bond donor atoms
  RbtUInt GetNumHBondDonorAtoms();
  RbtAtomList GetHBondDonorAtomList();

  //(Formally) charged atoms
  RbtUInt GetNumChargedAtoms();
  RbtAtomList GetChargedAtomList();

  // Planar atoms
  RbtUInt GetNumPlanarAtoms();
  RbtAtomList GetPlanarAtomList();

  // Binary

  // Atoms with atomic no = nAtomicNo
  RbtUInt GetNumAtomsWithAtomicNo_eq(RbtInt nAtomicNo);
  RbtAtomList GetAtomListWithAtomicNo_eq(RbtInt nAtomicNo);

  // Atoms with FFType = strFFType
  RbtUInt GetNumAtomsWithFFType_eq(std::string strFFType);
  RbtAtomList GetAtomListWithFFType_eq(std::string strFFType);

  ////////////////////////////////////////////
  // Bond list functions (provided for convenience, as user could just as well
  // call the Rbt:: functions with RbtModel::GetBondList)
  // e.g. RbtBondList bondList = Rbt::GetSelectedBondList(spModel->GetBondList);
  // is equivalent to
  // RbtBondList bondList = spModel->GetSelectedBondList();
  ////////////////////////////////////////////

  // Unary

  // Generic template version of GetNumBonds, passing in your own predicate
  template <class Predicate> RbtUInt GetNumBonds(const Predicate &pred);
  // Generic template version of GetBondList, passing in your own predicate
  template <class Predicate> RbtBondList GetBondList(const Predicate &pred);

  // Selected bonds
  void SetBondSelectionFlags(RbtBool bSelected = true);
  RbtUInt GetNumSelectedBonds();
  RbtBondList GetSelectedBondList();

  // Cyclic bonds
  void SetBondCyclicFlags(RbtBool bCyclic = true);
  RbtUInt GetNumCyclicBonds();
  RbtBondList GetCyclicBondList();

  // DM 10 Dec 1998 - at some point these functions should be implemented as
  // generic Rbt:: functions operating on arbitrary atom and bond lists

  // D Morley, 2 Dec 1998 - go back to the old way, it's more convenient to get
  // all donors in the same list and separate them later void
  // GetHBondDonorLists(RbtAtomList& donorList, RbtAtomList& donorHList);

  // Get min and max x,y,z coords
  // bInit = true => initialize min, max
  // bInit = false => don't initialize min, max (useful for accumulating min,max
  // over several RbtModels)
  void GetMinMaxCoords(RbtCoord &minCoord, RbtCoord &maxCoord,
                       RbtBool bInit = true);

  // Get map of (key=force field atom type string, value=no. of occurrences)
  RbtStringIntMap GetAtomTypeMap();

  // Get map of (key=force field bond type (atom type pair) string, value=no. of
  // occurrences)
  RbtStringIntMap GetBondTypeMap();

protected:
  //////////////////////
  // Private methods
  //////////////////////

private:
  RbtModel();                            // Disable default constructor
  RbtModel(const RbtModel &);            // Copy constructor disabled by default
  RbtModel &operator=(const RbtModel &); // Copy assignment disabled by default

  // Create a new model from a data source
  void Create(RbtBaseMolecularFileSource *pMolSource) throw(RbtError);
  void Clear();                         // Clear the current model
  void AddAtoms(RbtAtomList &atomList); // Register an atom list with the model

  //////////////////////
  // Private data
  //////////////////////
private:
  std::string m_strName;      // Model name
  RbtStringList m_titleList;  // Title list (read from file)
  RbtAtomList m_atomList;     // atom list
  RbtBondList m_bondList;     // bond list
  RbtSegmentMap m_segmentMap; // map of (key=segment name, value=atom count)
  RbtAtomListList m_ringList; //(DM 7 Dec 1998) list of atom lists for each ring
  RbtStringIntMap m_coordNames; //(DM 8 Feb 1999) map of named coord sets
                                //(key=name, value=index into RbtAtom coord map)
  RbtInt m_currentCoord;        // DM 11 Jul 2003 - which coord set is current
  RbtStringVariantMap m_dataMap; // DM 12 May 1999 - associated data (e.g. from
                                 // SD file or generated by rbdock)
  RbtPseudoAtomList
      m_pseudoAtomList; // DM 11 Jul 2000 - store associated pseudoatoms
  RbtModelMutatorPtr
      m_spMutator; // DM 30 Jul 2001 - helper object for fast dihedral mutations
  RbtFlexData *m_pFlexData;  // Stores the flexibility data for this model
  RbtChromElement *m_pChrom; // Reference chromosome. GetChrom() returns a clone
                             // of this object.
  RbtDouble
      m_occupancy;   // Occupancy value (0->1), in support of solvent occupancy
  RbtBool m_enabled; // Enabled state, depends on occupancy value and threshold
};

// Useful typedefs
typedef SmartPtr<RbtModel> RbtModelPtr;   // Smart pointer
typedef vector<RbtModelPtr> RbtModelList; // Vector of smart pointers
typedef RbtModelList::iterator RbtModelListIter;
typedef RbtModelList::const_iterator RbtModelListConstIter;

#endif //_RBTMODEL_H_
