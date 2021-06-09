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

#include "rxdock/PharmaSF.h"
#include "rxdock/FileError.h"
#include "rxdock/LigandError.h"
#include "rxdock/MdlFileSink.h"
#include "rxdock/WorkSpace.h"

#include <loguru.hpp>

using namespace rxdock;

// Static data members
const std::string PharmaSF::_CT = "PharmaSF";
const std::string PharmaSF::_CONSTRAINTS_FILE = "constraints-file";
const std::string PharmaSF::_OPTIONAL_FILE = "optional-file";
const std::string PharmaSF::_NOPT = "number-of-optional-parameters";
const std::string PharmaSF::_WRITE_ERRORS = "write-errors";

// NB - Virtual base class constructor (BaseSF) gets called first,
// implicit constructor for BaseInterSF is called second
PharmaSF::PharmaSF(const std::string &strName)
    : BaseSF(_CT, strName), m_nopt(0), m_bWriteErrors(false) {
  LOG_F(2, "PharmaSF parameterised constructor");
  // Add parameters It gets the right name in SetupReceptor
  AddParameter(_CONSTRAINTS_FILE, ".const");
  AddParameter(_OPTIONAL_FILE, "_opt.const");
  AddParameter(_NOPT, m_nopt);
  AddParameter(_WRITE_ERRORS, m_bWriteErrors);
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

PharmaSF::~PharmaSF() {
  LOG_F(2, "PharmaSF destructor");
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

void PharmaSF::SetupReceptor() {
  m_constrList.clear();
  m_optList.clear();
  m_conScores.clear();
  m_optScores.clear();

  if (GetReceptor().Null())
    return;
  std::string strWSName = GetWorkSpace()->GetName();
  std::string strConstraintFile =
      GetDataFileName("", GetParameter(_CONSTRAINTS_FILE));
  std::string strOptFile = GetDataFileName("", GetParameter(_OPTIONAL_FILE));

  // Create an output sink for ligands that do not have sufficient
  // pharmacophore features. These ligands will be written to an error SD file
  // but will not be docked.
  if (m_bWriteErrors) {
    MolecularFileSinkPtr p = GetWorkSpace()->GetSink();
    if (!p) {
      // This error will be thrown if the receptor model is connected
      // to the workspace BEFORE the ligand output sink
      // Should never happen with the current rbdock.cxx executable
      throw Error(_WHERE_, "Ligand output file sink undefined");
    }
    std::string strOutputFile(p->GetFileName());
    strOutputFile.erase(strOutputFile.find(".sd"), 3);
    strOutputFile += "_errors.sd";
    LOG_F(INFO,
          "PharmaSF::SetupReceptor: Reading mandatory ph4 constraints from {}",
          strConstraintFile);
    LOG_F(INFO,
          "PharmaSF::SetupReceptor: Ligands missing the mandatory features "
          "will be written to {}",
          strOutputFile);
    m_spErrorFile = MolecularFileSinkPtr(
        new MdlFileSink(strOutputFile.c_str(), GetLigand()));
  }
  LOG_F(INFO,
        "PharmaSF::SetupReceptor: Reading mandatory ph4 constraints from {}",
        strConstraintFile);
  std::ifstream constrFile(strConstraintFile.c_str(), std::ios_base::in);
  if (!constrFile) {
    throw FileReadError(_WHERE_, "cannot open mandatory constraints file " +
                                     strConstraintFile);
  }
  ZeroCounters();
  ReadConstraints(constrFile, m_constrList, true);
  constrFile.close();

  // Optional constraints
  std::ifstream optFile(strOptFile.c_str(), std::ios_base::in);
  if (optFile.good()) {
    LOG_F(INFO,
          "PharmaSF::SetupReceptor: Reading optional ph4 constraints from ",
          strOptFile);
    ReadConstraints(optFile, m_optList, false);
    // Keep m_nopt within range
    SetParameter(_NOPT, std::min(m_nopt, int(m_optList.size())));
    SetParameter(_NOPT, std::max(m_nopt, 0));
    LOG_F(INFO,
          "PharmaSF::SetupReceptor: {} of the optional constraints are "
          "required to be satisfied",
          m_nopt);
    optFile.close();
  } else {
    LOG_F(INFO,
          "PharmaSF::SetupReceptor: No optional ph4 constraints file found");
  }

  // Initialise the component score vectors
  m_conScores = std::vector<double>(m_constrList.size(), 0.0);
  m_optScores = std::vector<double>(m_optList.size(), 0.0);
}

void PharmaSF::SetupLigand() {
  if (GetLigand().Null())
    return;
  try {
    LOG_F(INFO, "PharmaSF::SetupReceptor: Checking mandatory ph4 features...");
    for (ConstraintListIter iter = m_constrList.begin();
         iter != m_constrList.end(); iter++) {
      (*iter)->AddAtomList(GetLigand(), true);
    }
    LOG_F(INFO, "PharmaSF::SetupReceptor: All mandatory features found");
  } catch (LigandError &e) {
    if (m_bWriteErrors) {
      m_spErrorFile->SetModel(GetLigand());
      m_spErrorFile->Render();
    }
    throw;
  }

  // Update the optional features
  for (ConstraintListIter iter = m_optList.begin(); iter != m_optList.end();
       iter++) {
    (*iter)->AddAtomList(GetLigand(), false);
  }
}

void PharmaSF::SetupScore() {
  // No further setup required
}

double PharmaSF::RawScore() const {
  double total(0.0);
  // Store and sum all the mandatory terms
  int i = 0;
  for (ConstraintListConstIter iter = m_constrList.begin();
       iter != m_constrList.end(); iter++, i++) {
    m_conScores[i] = (*iter)->Score();
  }
  total = std::accumulate(m_conScores.begin(), m_conScores.end(), total);

  // Store the optional scores for sorting
  i = 0;
  for (ConstraintListConstIter iter = m_optList.begin();
       iter != m_optList.end(); iter++, i++) {
    m_optScores[i] = (*iter)->Score();
  }
  // partial_sort_copy copies the N lowest scores
  std::vector<double> lowest(m_nopt);
  std::partial_sort_copy(m_optScores.begin(), m_optScores.end(), lowest.begin(),
                         lowest.end());
  LOG_F(1, "PharmaSF::RawScore: {} lowest optional scores: {}", m_nopt,
        fmt::join(lowest, ", "));
  total = std::accumulate(lowest.begin(), lowest.end(), total);
  return total;
}

// DM 25 Oct 2000 - track changes to parameter values in local data members
// ParameterUpdated is invoked by ParamHandler::SetParameter
void PharmaSF::ParameterUpdated(const std::string &strName) {
  if (strName == _WRITE_ERRORS) {
    m_bWriteErrors = GetParameter(_WRITE_ERRORS);
  } else if (strName == _NOPT) {
    m_nopt = GetParameter(_NOPT);
  }
  BaseSF::ParameterUpdated(strName);
}

// Override BaseSF::ScoreMap to provide additional raw descriptors
void PharmaSF::ScoreMap(StringVariantMap &scoreMap) const {
  if (isEnabled()) {
    // Copied from BaseSF
    double rs = RawScore();
    std::string name = GetFullName();
    scoreMap[name] = rs;
    AddToParentMapEntry(scoreMap, rs);
    // Store the mandatory constraint scores
    for (unsigned int i = 0; i < m_conScores.size(); i++) {
      std::ostringstream field;
      field << name << ".con_" << i + 1;
      scoreMap[field.str()] = m_conScores[i];
    }
    // Store the optional constraint scores (unsorted)
    for (unsigned int i = 0; i < m_optScores.size(); i++) {
      std::ostringstream field;
      field << name << ".opt_" << i + 1;
      scoreMap[field.str()] = m_optScores[i];
    }
  }
}
