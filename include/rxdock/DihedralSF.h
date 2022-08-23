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

// Base implementation class for all dihedral scoring functions. Provides
// Dihedral struct for dihedral definitions and score, plus a method for
// generating lists of dihedral definitions from a bond list.
// This version hardwired to use Tripos 5.2 atom types and params.

#ifndef _RBTDIHEDRALSF_H_
#define _RBTDIHEDRALSF_H_

#include "rxdock/Atom.h"
#include "rxdock/BaseSF.h"
#include "rxdock/Bond.h"
#include "rxdock/ParameterFileSource.h"
#include "rxdock/TriposAtomType.h"

#include <nlohmann/json.hpp>

using json = nlohmann::json;

namespace rxdock {

// Class for holding dihedral atom specifiers and force field params
// Main method is operator() - calculates single dihedral score for current
// conformation DM (9 Dec 2003): Extended to allow additional terms to be added
// to the dihedral potential (AddTerm) Two main uses for AddTerm:
//  i) to define a Fourier series expansion, rather than a simple cosine
//  potential
// ii) implicit hydrogen correction (ghost terms are added for the missing
// hydrogen dihedrals
//                                   defined as offsets from the real heavy-atom
//                                   dihedral specifiers)

class DihedralElement {
public:
  // Simple struct for storing the force field params for a single term of the
  // dihedral potential
  class prms {
  public:
    prms(double ss = 3.0, double kk = 0.2, double o = 0.0)
        : s(std::fabs(ss)), k(kk), sign((ss > 0.0) ? +1.0 : -1.0), offset(o) {}
    double s;      // rotational degeneracy
    double k;      // barrier height
    double sign;   //+1 or -1, to invert potential
    double offset; // dihedral offset in degrees (useful for ghost implicit
                   // hydrogens)
  };

  // Constructor takes the real atom specifiers, plus the first term of the
  // potential
  DihedralElement(Atom *pAtom1, Atom *pAtom2, Atom *pAtom3, Atom *pAtom4,
                  const prms &dihprms);

  friend void to_json(json &j, const DihedralElement &dihedralElem);
  friend void from_json(const json &j, DihedralElement &dihedralElem);

  double operator()() const; // Calculate dihedral score for this interaction
  Atom *GetAtom1Ptr() const { return m_pAtom1; }
  Atom *GetAtom2Ptr() const { return m_pAtom2; }
  Atom *GetAtom3Ptr() const { return m_pAtom3; }
  Atom *GetAtom4Ptr() const { return m_pAtom4; }
  // Add additional terms to the potential
  void AddTerm(const prms &dihprms);

private:
  DihedralElement(); // Forbid default constructor
  Atom *m_pAtom1;
  Atom *m_pAtom2;
  Atom *m_pAtom3;
  Atom *m_pAtom4;
  std::vector<prms> m_prms;
};

void to_json(json &j, const DihedralElement &dihedralElem);
void from_json(const json &j, DihedralElement &dihedralElem);

// Useful typedefs
typedef std::vector<DihedralElement *> DihedralList; // Vector of smart pointers
typedef DihedralList::iterator DihedralListIter;
typedef DihedralList::const_iterator DihedralListConstIter;

class DihedralSF : public virtual BaseSF {
public:
  // Class type string
  static const std::string _CT;
  // Parameter names
  static const std::string _IMPL_H_CORR; // DM 6 Aug 2002 - option to correct
                                         // barrier heights for implicit H's

  virtual ~DihedralSF();

  friend void to_json(json &j, const DihedralSF &dihedralSF);
  friend void from_json(const json &j, DihedralSF &dihedralSF);

protected:
  DihedralSF();
  // Creates a vector of pointers to dihedral objects from a supplied bond list
  // Dihedral params are set using Tripos 5.2 atom types + params
  // NOTE: It is the responsibility of the subclass to delete the dihedral
  // objects which are created
  DihedralList CreateDihedralList(const BondList &bondList);

private:
  // Lookup the dihedral params for a given set of types
  DihedralElement::prms FindDihedralParams(TriposAtomType::eType t1,
                                           TriposAtomType::eType t2,
                                           TriposAtomType::eType t3,
                                           TriposAtomType::eType t4);

  // Determines the list of bonded atoms for pAtom1, and the list of dihedral
  // offsets for any implicit hydrogens pAtom2 is the other central atom in the
  // bond, is excluded from the returned list of bonded atoms
  void CalcBondedAtoms(Atom *pAtom1, Atom *pAtom2, AtomList &bondedAtoms,
                       std::vector<double> &offsets);

  ParameterFileSourcePtr m_spDihedralSource;
  std::vector<std::string> m_centralPairs;
  TriposAtomType m_triposType;
};

void to_json(json &j, const DihedralSF &dihedralSF);
void from_json(const json &j, DihedralSF &dihedralSF);

} // namespace rxdock

#endif //_RBTDIHEDRALSF_H_
