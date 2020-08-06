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

// Atom typing class for Tripos 5.2 forcefield

#ifndef _RBTTRIPOSATOMTYPE_H_
#define _RBTTRIPOSATOMTYPE_H_

#include "rxdock/Config.h"

namespace rxdock {

class Atom;

class TriposAtomType {
public:
  ///////////////////////////////////////////////
  // Enums
  ///////////////////////////////////////////////

  // Atom types
  enum eType {
    UNDEFINED = 0,
    Al,
    Br,
    C_cat,
    C_1,
    C_1_H1,
    C_2,
    C_2_H1,
    C_2_H2,
    C_3,
    C_3_H1,
    C_3_H2,
    C_3_H3,
    C_ar,
    C_ar_H1,
    Ca,
    Cl,
    Du,
    F,
    H,
    H_P,
    I,
    K,
    Li,
    LP,
    N_1,
    N_2,
    N_3,
    N_4,
    N_am,
    N_ar,
    N_pl3,
    Na,
    O_2,
    O_3,
    O_co2,
    P_3,
    S_2,
    S_3,
    S_o,
    S_o2,
    Si,
    MAXTYPES // KEEP AS LAST TYPE: used to size the atom name string list
  };

  struct info {
    std::string name;
    int atomicNo;
    int hybrid;
    info() : atomicNo(0), hybrid(0) {}
    info(const std::string &n, int a, int h)
        : name(n), atomicNo(a), hybrid(h) {}
  };

  ////////////////////////////////////////
  // Constructors/destructors

  RBTDLL_EXPORT TriposAtomType();
  RBTDLL_EXPORT virtual ~TriposAtomType();

  ////////////////////////////////////////
  // Public methods
  ////////////////

  // Returns Tripos atom type for given atom
  // If useExtendedTypes is true, pseudo types for extended carbons and polar Hs
  // are used
  eType operator()(Atom *pAtom, bool useExtendedTypes = false) const;
  // Converts Tripos type to string
  RBTDLL_EXPORT std::string Type2Str(eType) const;
  // Get hybridisation from Tripos type
  int Type2Hybrid(eType) const;
  // Get atomic number from Tripos type
  int Type2AtomicNo(eType) const;
  // Converts string to Tripos type
  RBTDLL_EXPORT eType Str2Type(const std::string &) const;

protected:
  ////////////////////////////////////////
  // Protected methods
  ///////////////////

private:
  ////////////////////////////////////////
  // Private methods
  /////////////////
  void SetupTypeInfo();

protected:
  ////////////////////////////////////////
  // Protected data
  ////////////////

private:
  ////////////////////////////////////////
  // Private data
  //////////////
  std::vector<info> m_typeInfo;
};

// Useful typedefs
typedef std::vector<TriposAtomType::info> TriposInfoList;
typedef TriposInfoList::iterator TriposInfoListIter;
typedef TriposInfoList::const_iterator TriposInfoListConstIter;

typedef std::vector<TriposAtomType::eType> TriposAtomTypeList;
typedef TriposAtomTypeList::iterator TriposAtomTypeListIter;
typedef TriposAtomTypeList::const_iterator TriposAtomTypeListConstIter;

} // namespace rxdock

#endif //_RBTTRIPOSATOMTYPE_H_
