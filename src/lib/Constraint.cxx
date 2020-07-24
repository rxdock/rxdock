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

#include "Constraint.h"
#include "LigandError.h"
#include <functional>
#include <sstream>

using namespace rxdock;

// initialization of the static data of Constraint
std::string Constraint::_CT("Constraint");
unsigned int HeavyConstraint::counter =
    0; // 7 Feb 2005 (DM) new constraint type
unsigned int HBAConstraint::counter = 0;
unsigned int HBDConstraint::counter = 0;
unsigned int HydroConstraint::counter = 0;
unsigned int HydroAliphaticConstraint::counter = 0;
unsigned int HydroAromaticConstraint::counter = 0;
unsigned int NegChargeConstraint::counter = 0;
unsigned int PosChargeConstraint::counter = 0;
unsigned int RingAromaticConstraint::counter = 0;

Constraint::Constraint(Coord c, double t) {
  m_atomList.clear();
  coord = c;
  tolerance = t;
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

Constraint::Constraint(const Constraint &c) {
  coord = c.coord;
  tolerance = c.tolerance;
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

Constraint::Constraint() {
  coord = Coord(0, 0, 0);
  tolerance = 0.0;
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

Constraint::~Constraint() { _RBTOBJECTCOUNTER_DESTR_(_CT); }

void Constraint::copy(const Constraint &c) {
  coord = c.coord;
  tolerance = c.tolerance;
}

std::ostream &Constraint::Print(std::ostream &o) const {
  o << coord << "\t" << tolerance;
  return o;
}

std::ostream &operator<<(std::ostream &s, const Constraint &c) {
  return (c.Print(s));
}

// DM 17 Oct 2003 - rewrite to avoid step increase in score at the tolerance
// distance Allow for fact that atomlist may be empty
double Constraint::Score() const {
  if (m_atomList.empty()) {
    // Return artifically high score if the ligand has none of the desired
    // features This will only apply to the optional constraints (mandatory
    // constraints are prefiltered)
    return 99.0;
  }

  // Compile a list of all the distance**2 to the constraint center
  std::vector<double> d2list;
  for (AtomListConstIter iter = m_atomList.begin(); iter != m_atomList.end();
       iter++) {
    d2list.push_back(Length2((*iter)->GetCoords(), coord));
  }
  // Min distance**2
  double d2min = *(std::min_element(d2list.begin(), d2list.end()));
  // Score is zero up to the tolerance range, then quadratic beyond.
  double dr = std::sqrt(d2min) - tolerance;
  return (dr > 0.0) ? dr * dr : 0.0;
}

ConstraintPtr rxdock::CreateConstraint(Coord &c, double &t, std::string &n,
                                       bool bCount) {
  if (n == "Any") {
    if (bCount)
      HeavyConstraint::counter++;
    return new HeavyConstraint(c, t);
  } else if (n == "Don") {
    if (bCount)
      HBDConstraint::counter++;
    return new HBDConstraint(c, t);
  } else if (n == "Acc") {
    if (bCount)
      HBAConstraint::counter++;
    return new HBAConstraint(c, t);
  } else if (n == "Hyd") {
    if (bCount)
      HydroConstraint::counter++;
    return new HydroConstraint(c, t);
  } else if (n == "Ani") {
    if (bCount)
      NegChargeConstraint::counter++;
    return new NegChargeConstraint(c, t);
  } else if (n == "Cat") {
    if (bCount)
      PosChargeConstraint::counter++;
    return new PosChargeConstraint(c, t);
  } else if (n == "Aro") {
    if (bCount)
      RingAromaticConstraint::counter++;
    return new RingAromaticConstraint(c, t);
  } else if (n == "Hal") {
    if (bCount)
      HydroAliphaticConstraint::counter++;
    return new HydroAliphaticConstraint(c, t);
  } else if (n == "Har") {
    if (bCount)
      HydroAromaticConstraint::counter++;
    return new HydroAromaticConstraint(c, t);
  } else
    throw Error(_WHERE_, "Constraint " + n + " not recognized");
}

void rxdock::ZeroCounters() {
  HeavyConstraint::counter = 0;
  HBAConstraint::counter = 0;
  HBDConstraint::counter = 0;
  HydroConstraint::counter = 0;
  HydroAliphaticConstraint::counter = 0;
  HydroAromaticConstraint::counter = 0;
  NegChargeConstraint::counter = 0;
  PosChargeConstraint::counter = 0;
  RingAromaticConstraint::counter = 0;
}

void rxdock::ReadConstraint(std::istream &ifile, ConstraintPtr &cnt,
                            bool bCount) {
  std::string n;
  Coord c;
  double t;
  ifile >> c;
  ifile >> t;
  ifile >> n;
  std::cout << c << "\t" << t << "\t" << n << std::endl;
  if (ifile.fail())
    throw Error(_WHERE_, "Problems reading constraint " + n);
  cnt = CreateConstraint(c, t, n, bCount);
}

void rxdock::ReadConstraintFromMoe(std::istream &ifile, ConstraintPtr &cnt,
                                   bool bCount) {
  std::string n, extra;
  double x, y, z, t;
  ifile >> n;
  std::cout << n << std::endl;
  ifile >> extra;
  ifile >> x >> y >> z >> t;
  ifile >> extra >> extra;
  if (ifile.fail())
    throw Error(_WHERE_, "Problems reading constraint " + n);
  Coord c(x, y, z);
  std::cout << c << "\t" << t << "\t" << n << std::endl;
  cnt = CreateConstraint(c, t, n, bCount);
}

void rxdock::ReadStartMoe(std::istream &ifile) {
  const int MAXLINELENGTH = 100; // I'm quite sure the max length for any
                                 // line in moe is 72, but
                                 // prefer to be on the safe side
  char *line = new char[MAXLINELENGTH + 1];
  ConstraintPtr cons;
  char ch;
  ifile.getline(line, MAXLINELENGTH); // get rid of first four lines
  ifile.getline(line, MAXLINELENGTH);
  ifile.getline(line, MAXLINELENGTH);
  // This line has some data that could be neccesary later on, like
  // partial match information also, it contains the
  // values for the fields comment and title, that can be as
  // long as necessary. Need to make sure I read all the lines
  // before the next commentted line (starting with a #)
  while ((ch = ifile.peek()) != '#')
    ifile.getline(line, MAXLINELENGTH);
  ifile.getline(line, MAXLINELENGTH); // next comment line
}

void rxdock::ReadConstraints(std::istream &ifile, ConstraintList &cl,
                             bool bCount) {
  char ch;
  bool moe = false;
  cl.clear();
  ConstraintPtr cons;
  ch = ifile.peek();
  if (ch == '#') // moe file
  {
    moe = true;
    ReadStartMoe(ifile);
  }
  while ((!moe && !ifile.eof()) || (moe && ((ch = ifile.peek()) != '#'))) {
    if (moe)
      ReadConstraintFromMoe(ifile, cons, bCount);
    else
      ReadConstraint(ifile, cons, bCount);
    cl.push_back(cons);
    while (ifile.get(ch)) {
      if (!isspace(ch)) {
        ifile.putback(ch);
        break;
      }
    }
  }
  if (cl.size() == 0)
    std::cout << Constraint::_CT
              << "** WARNING: number of constraints read is 0\n";
}

// 07 Feb 2005 (DM) - new constraint type, any heavy atom
void HeavyConstraint::AddAtomList(ModelPtr lig, bool bCheck) {
  m_atomList =
      GetAtomListWithPredicate(lig->GetAtomList(), std::not1(isAtomicNo_eq(1)));
  if (bCheck && (m_atomList.size() < counter)) {
    std::ostringstream ostr;
    ostr << "The ligand has only " << m_atomList.size() << " heavy atom(s) ("
         << counter << " required)" << std::ends;
    throw LigandError(_WHERE_, ostr.str());
  }
}

// 16 May 2003 (DM) - limit to neutral H-bond donor hydrogens
void HBDConstraint::AddAtomList(ModelPtr lig, bool bCheck) {
  m_atomList = GetAtomListWithPredicate(lig->GetAtomList(), isAtomHBondDonor());
  m_atomList =
      GetAtomListWithPredicate(m_atomList, std::not1(isAtomCationic()));
  if (bCheck && (m_atomList.size() < counter)) {
    std::ostringstream ostr;
    ostr << "The ligand has only " << m_atomList.size()
         << " neutral H-bond donor hydrogens(s) (" << counter << " required)"
         << std::ends;
    throw LigandError(_WHERE_, ostr.str());
  }
}

// 16 May 2003 (DM) - limit to neutral H-bond acceptors
void HBAConstraint::AddAtomList(ModelPtr lig, bool bCheck) {
  m_atomList =
      GetAtomListWithPredicate(lig->GetAtomList(), isAtomHBondAcceptor());
  m_atomList = GetAtomListWithPredicate(m_atomList, std::not1(isAtomAnionic()));
  if (bCheck && (m_atomList.size() < counter)) {
    std::ostringstream ostr;
    ostr << "The ligand has only " << m_atomList.size()
         << " neutral H-bond acceptor(s) (" << counter << " required)"
         << std::ends;
    throw LigandError(_WHERE_, ostr.str());
  }
}

void HydroConstraint::AddAtomList(ModelPtr lig, bool bCheck) {
  m_atomList = GetAtomListWithPredicate(lig->GetAtomList(), isAtomLipophilic());
  if (bCheck && (m_atomList.size() < counter)) {
    std::ostringstream ostr;
    ostr << "The ligand has only " << m_atomList.size()
         << " hydrophobic atom(s) (" << counter << " required)" << std::ends;
    throw LigandError(_WHERE_, ostr.str());
  }
}

void HydroAliphaticConstraint::AddAtomList(ModelPtr lig, bool bCheck) {
  m_atomList = GetAtomListWithPredicate(lig->GetAtomList(), isAtomLipophilic());
  m_atomList =
      GetAtomListWithPredicate(m_atomList, isHybridState_eq(Atom::SP3));
  if (bCheck && (m_atomList.size() < counter)) {
    std::ostringstream ostr;
    ostr << "The ligand has only " << m_atomList.size()
         << " hydrophobic aliphatic atom(s) (" << counter << " required)"
         << std::ends;
    throw LigandError(_WHERE_, ostr.str());
  }
}

void HydroAromaticConstraint::AddAtomList(ModelPtr lig, bool bCheck) {
  m_atomList = GetAtomListWithPredicate(lig->GetAtomList(), isAtomLipophilic());
  m_atomList =
      GetAtomListWithPredicate(m_atomList, isHybridState_eq(Atom::AROM));
  if (bCheck && (m_atomList.size() < counter)) {
    std::ostringstream ostr;
    ostr << "The ligand has only " << m_atomList.size()
         << " hydrophobic aromatic atom(s) (" << counter << " required)"
         << std::ends;
    throw LigandError(_WHERE_, ostr.str());
  }
}

void NegChargeConstraint::AddAtomList(ModelPtr lig, bool bCheck) {
  m_atomList = GetAtomListWithPredicate(lig->GetAtomList(), isAtomAnionic());
  if (bCheck && (m_atomList.size() < counter)) {
    std::ostringstream ostr;
    ostr << "The ligand has only " << m_atomList.size()
         << " negatively charged atom(s) (" << counter << " required)"
         << std::ends;
    throw LigandError(_WHERE_, ostr.str());
  }
}

void PosChargeConstraint::AddAtomList(ModelPtr lig, bool bCheck) {
  m_atomList = GetAtomListWithPredicate(lig->GetAtomList(), isAtomCationic());
  if (bCheck && (m_atomList.size() < counter)) {
    std::ostringstream ostr;
    ostr << "The ligand has only " << m_atomList.size()
         << " positively charged atom(s) (" << counter << " required)"
         << std::ends;
    throw LigandError(_WHERE_, ostr.str());
  }
}

void RingAromaticConstraint::AddAtomList(ModelPtr lig, bool bCheck) {
  m_atomList.clear();
  AtomList at = lig->GetAtomList();
  AtomListList ligRingLists = lig->GetRingAtomLists();
  for (AtomListListConstIter rIter = ligRingLists.begin();
       rIter != ligRingLists.end(); rIter++) {
    if (GetNumAtomsWithPredicate(*rIter, isPiAtom()) == (*rIter).size()) {
      PseudoAtomPtr spPseudoAtom = lig->AddPseudoAtom(*rIter);
      m_atomList.push_back(spPseudoAtom);
    }
  }
  if (bCheck && (m_atomList.size() < counter)) {
    std::ostringstream ostr;
    ostr << "The ligand has only " << m_atomList.size() << " aromatic ring(s) ("
         << counter << " required)" << std::ends;
    throw LigandError(_WHERE_, ostr.str());
  }
}
