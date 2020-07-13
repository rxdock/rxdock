// Smarts class. Defines a function that takes a model and a smarts string as
// input, and returns a list of lists of atoms that match the pattern
// !!!!!!!!!!!!!!!!!!!       WARNING        !!!!!!!!!!!!!!!!!!!!!!
// This subroutines assume that atomId starts from 1 and goes incrementally
// It will fail if this is not true.

#include "RbtSmarts.h"
#include "dt_smarts.h"
#include "dt_smiles.h"

void print_atoms_in_path(dt_Handle path);

RbtAtomListList DT::QueryModel(RbtModelPtr spModel, const std::string &strSmart,
                               std::string &strSmiles) {
  Rbt::isAtomicNo_eq isO(8);
  Rbt::isAtomicNo_eq isH(1);
  RbtAtomList atomList = spModel->GetAtomList();
  RbtBondList bondList = spModel->GetBondList();
  RbtIntIntMap r2d; // track rDock to Daylight numbering
  RbtIntIntMap d2r; // track Daylight to rDock numbering

  // Convert rDock RbtModel object to daylight molecule (with all hydrogens
  // implicit, even the polar ones)
  dt_Handle m = dt_alloc_mol();
  // keep a vector of all the Daylight atom handles we create
  // so we can create the bonds correctly later
  vector<dt_Handle> atoms;
  dt_mod_on(m);
  // First create the atoms
  for (RbtAtomListConstIter at = atomList.begin(); at != atomList.end(); at++) {
    // ignore all hydrogens
    if (isH(*at)) {
      // std::cout << "Atom: " << (*at)->GetName() << "; ignored" <<
      // std::endl;
      continue;
    }
    int rID = at - atomList.begin(); // The rDock atom numbering (sequence in
                                     // atomList vector, from zero)
    int dID = atoms.size(); // The Daylight numbering (sequence in atoms
                            // vector, from zero)
    r2d[rID] = dID;
    d2r[dID] = rID;
    atoms.push_back(dt_addatom(m, (*at)->GetAtomicNo(),
                               (*at)->GetNumImplicitHydrogens() +
                                   (*at)->GetCoordinationNumber(1)));
    dt_setarborder(atoms.back(), dID);
    // Workaround for terminal -ve charged acids
    // We want to use the SMILES string we generate from the dt_mol
    // as a SMARTS query that will deal with the symmetry of C(=O)[O-]
    // i.e. both oxygens are intended to be equivalent.
    // Do this by not defining the charges or bond orders for the oxygens.
    // Check for NOT (oxygen with negative partial "group" charge) before
    // defining formal charge "group" charge is the rDock distributed formal
    // charge
    if (!(isO(*at) && ((*at)->GetGroupCharge() < 0.0) &&
          ((*at)->GetGroupCharge() > -1.0))) {
      dt_setcharge(atoms.back(), (*at)->GetFormalCharge());
    }
    // std::cout << "Atom: " << (*at)->GetName() << "; rID = " << rID << ";
    // dID = " << dID << std::endl;
  }

  // Now the bonds
  for (RbtBondListConstIter bd = bondList.begin(); bd != bondList.end(); bd++) {
    RbtAtomPtr at1 = (*bd)->GetAtom1Ptr();
    RbtAtomPtr at2 = (*bd)->GetAtom2Ptr();
    // ignore all bonds to hydrogen
    if (isH(at1) || isH(at2)) {
      // std::cout << "Bond: " << at1->GetName() << "-" <<
      // at2->GetName() <<
      // "; ignored" << std::endl;
      continue;
    }
    // Try and find the bonded atoms in the rDock atom list (comparison is on
    // RbtAtom* pointers)
    RbtAtomListConstIter iter1 = std::find_if(
        atomList.begin(), atomList.end(),
        std::bind(Rbt::isAtomPtr_eq(), std::placeholders::_1, at1));
    RbtAtomListConstIter iter2 = std::find_if(
        atomList.begin(), atomList.end(),
        std::bind(Rbt::isAtomPtr_eq(), std::placeholders::_1, at2));
    Assert<RbtAssert>((iter1 != atomList.end()) && (iter2 != atomList.end()));
    // rDock numbering
    int rID1 = iter1 - atomList.begin();
    int rID2 = iter2 - atomList.begin();
    Assert<RbtAssert>((r2d.find(rID1) != r2d.end()) &&
                      (r2d.find(rID2) != r2d.end()));
    // Use lookup to get the Daylight numbering
    int dID1 = r2d[rID1];
    int dID2 = r2d[rID2];
    // std::cout << "Bond: " << at1->GetName() << "-" << at2->GetName()
    //	 << "; rID = " << rID1 << "-" << rID2 << "; dID = " << dID1 << "-" <<
    // dID2 << std::endl; Finally create the Daylight bond Define formal bond
    // order except for terminal charged acids, see above
    dt_Handle b =
        dt_addbond(atoms[dID1], atoms[dID2], (*bd)->GetFormalBondOrder());
    if ((isO(at1) && (at1->GetGroupCharge() < 0.0) &&
         (at1->GetGroupCharge() > -1.0)) ||
        (isO(at2) && (at2->GetGroupCharge() < 0.0) &&
         (at2->GetGroupCharge() > -1.0))) {
      dt_setbondorder(b, DX_BTY_UNKNOWN);
    }
  }

  dt_mod_off(m);

  // Determine SMILES string for molecule
  char *arbsmi;
  int len;
  arbsmi = dt_arbsmiles(&len, m, 0);
  strSmiles = arbsmi; // strSmiles is returned by reference

  // Convert input SMARTS query string to daylight smarts object
  // If query string is empty, then use the SMILES string for the molecule
  // i.e. query the molecule with itself, to determine all symmetry-related
  // pathsets
  char *smart;
  if (strSmart.empty()) {
    smart = new char[strSmiles.length() + 1];
    strSmiles.copy(smart, std::string::npos);
    smart[strSmiles.length()] = 0;
  } else {
    smart = new char[strSmart.length() + 1];
    strSmart.copy(smart, std::string::npos);
    smart[strSmart.length()] = 0;
  }

  // Run the query and convert all the pathsets to RbtAtom lists
  dt_Handle sm = dt_smartin(strlen(smart), smart);
  dt_Handle pathset = dt_match(sm, m, 0);
  dt_Handle path;
  RbtAtomListList atLL;
  if (NULL_OB != pathset) {
    dt_Handle paths = dt_stream(pathset, TYP_PATH);
    while (NULL_OB != (path = dt_next(paths))) {
      // std::printf("match .......");
      // print_atoms_in_path(path);
      // std::printf("\n");
      dt_Handle atoms_in_path = dt_stream(path, TYP_ATOM);
      RbtAtomList atL;
      dt_Handle atom;
      while (NULL_OB != (atom = dt_next(atoms_in_path))) {
        int dID = dt_arborder(atom);
        Assert<RbtAssert>((d2r.find(dID) != d2r.end()));
        int rID = d2r[dID];
        atL.push_back(atomList[rID]);
      }
      atLL.push_back(atL);
    }
  }

  return atLL;
  dt_dealloc(m);
  dt_dealloc(pathset);
  dt_dealloc(sm);
  // REMEMBER TO DELETE CHAR* strings
}

/*============================================================================
 *  print_atoms_in_path() -- print atoms in the given path
 *
 *  Actually, this is polymorphic; works for any atom-containing object.
 */

void print_atoms_in_path(dt_Handle path) {
  dt_Handle mol, atom, atoms;
  char *str, symbol[10];
  int lens;

  atoms = dt_stream(path, TYP_ATOM);
  while (NULL_OB != (atom = dt_next(atoms))) {
    str = dt_symbol(&lens, atom);
    strncpy(symbol, str, lens);
    if (dt_aromatic(atom))
      *symbol = *symbol + 'a' - 'A';
    std::printf(" %.*s(%d)", lens, str, dt_arborder(atom));
  }
  dt_dealloc(atoms);
}
