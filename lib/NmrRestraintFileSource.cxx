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

#include "rxdock/NmrRestraintFileSource.h"
#include "rxdock/FileError.h"

using namespace rxdock;

////////////////////////////////////////
// Constructors/destructors
NmrRestraintFileSource::NmrRestraintFileSource(const std::string &fileName)
    : BaseFileSource(fileName) {
  _RBTOBJECTCOUNTER_CONSTR_("NmrRestraintFileSource");
}

// Destructor
NmrRestraintFileSource::~NmrRestraintFileSource() {
  ClearRestraintCache();
  _RBTOBJECTCOUNTER_DESTR_("NmrRestraintFileSource");
}

////////////////////////////////////////
// Public methods
////////////////

// Number of NOE restraints in file
unsigned int NmrRestraintFileSource::GetNumNoeRestraints() {
  Parse();
  return m_noeRestraintList.size();
}

// List of NOE restraints
NoeRestraintNamesList NmrRestraintFileSource::GetNoeRestraintList() {
  Parse();
  return m_noeRestraintList;
}

// Number of STD restraints in file
unsigned int NmrRestraintFileSource::GetNumStdRestraints() {
  Parse();
  return m_stdRestraintList.size();
}

// List of STD restraints
StdRestraintNamesList NmrRestraintFileSource::GetStdRestraintList() {
  Parse();
  return m_stdRestraintList;
}

////////////////////////////////////////
// Private methods
/////////////////

// Pure virtual in BaseFileSource - needs to be defined here
void NmrRestraintFileSource::Parse() {
  const std::string strDelimiter(","); // DM 23 Sept 1999

  // Only parse if we haven't already done so
  if (!m_bParsedOK) {
    ClearRestraintCache(); // Clear current cache
    Read();                // Read the file

    try {
      for (FileRecListIter fileIter = m_lineRecs.begin();
           fileIter != m_lineRecs.end(); fileIter++) {
        // Ignore blank lines and comment lines
        if (((*fileIter).length() == 0) || ((*fileIter).at(0) == '#')) {
          continue;
        }
        // Read the restraint
        std::string strAtomNames1;
        std::string strAtomNames2;
        double maxDist(0.0);
        std::istringstream istr(*fileIter);
        istr >> strAtomNames1 >> strAtomNames2 >> maxDist;
        // Check if we read all the fields OK
        if (!istr)
          throw FileParseError(_WHERE_, "Missing field(s) at line " +
                                            (*fileIter) + " in " +
                                            GetFileName());

        // DM 10 Dec 2002 - support for STD restraints. The first atom name
        // string is replaced by the keyword STD (or std) in the file. The
        // restraint is applied to the second atom list and implies that the
        // specified atoms should be within maxDist of *any* part of the
        // receptor surface
        if (strAtomNames1 == "std" || strAtomNames1 == "STD") {
          // Add new restraint to the list
          StdRestraintNames restraint;
          restraint.from.type = NoeRestraintType(strAtomNames2);
          // Check for bad syntax (restraint type=UNDEFINED)
          if (restraint.from.type == NOE_UNDEFINED)
            throw FileParseError(_WHERE_,
                                 "Unmatched or misplaced brackets at line " +
                                     (*fileIter) + " in " + GetFileName());
          std::vector<std::string> fromList =
              ConvertDelimitedStringToList(strAtomNames2, strDelimiter);
          restraint.from.names = fromList;
          restraint.maxDist = maxDist;
          m_stdRestraintList.push_back(restraint);
        }

        // NOE restraints
        else {
          // Add new restraint to the list
          NoeRestraintNames restraint;
          restraint.from.type = NoeRestraintType(strAtomNames1);
          restraint.to.type = NoeRestraintType(strAtomNames2);
          // Check for bad syntax (restraint type=UNDEFINED)
          if ((restraint.from.type == NOE_UNDEFINED) ||
              (restraint.to.type == NOE_UNDEFINED))
            throw FileParseError(_WHERE_,
                                 "Unmatched or misplaced brackets at line " +
                                     (*fileIter) + " in " + GetFileName());
          std::vector<std::string> fromList =
              ConvertDelimitedStringToList(strAtomNames1, strDelimiter);
          restraint.from.names = fromList;
          std::vector<std::string> toList =
              ConvertDelimitedStringToList(strAtomNames2, strDelimiter);
          restraint.to.names = toList;
          restraint.maxDist = maxDist;
          m_noeRestraintList.push_back(restraint);
        }
      }
      //////////////////////////////////////////////////////////
      // If we get this far everything is OK
      m_bParsedOK = true;
    }

    catch (Error &error) {
      ClearRestraintCache();
      throw; // Rethrow the Error
    }
  }
}

void NmrRestraintFileSource::ClearRestraintCache() {
  m_noeRestraintList.clear();
  m_stdRestraintList.clear();
}

// Returns NOE restraint type and modifies the atom name string accordingly
// Returns UNDEFINED if the atom name string has bad syntax
eNoeType NmrRestraintFileSource::NoeRestraintType(std::string &strAtomNames) {
  eNoeType restraintType(NOE_UNDEFINED);

  // Restraint type: NOE_MEAN
  // Syntax: (at1,at2,at3...)
  // Check if atom names are bracketed correctly, if so change type to MEAN, and
  // remove the brackets
  std::string::size_type ob = strAtomNames.find("(");
  std::string::size_type cb = strAtomNames.find(")");
  std::string::size_type iLast = strAtomNames.length() - 1;
  if ((ob == 0) && (cb == iLast)) {
    strAtomNames.erase(iLast, 1);
    strAtomNames.erase(0, 1);
    restraintType = NOE_MEAN;
  } else if ((ob == std::string::npos) && (cb == std::string::npos)) {
    // Restraint type: NOE_AND
    // Syntax: [at1,at2,at3...]
    // Check if atom names are bracketed correctly, if so change type to AND,
    // and remove the brackets
    std::string::size_type osqb = strAtomNames.find("[");
    std::string::size_type csqb = strAtomNames.find("]");
    if ((osqb == 0) && (csqb == iLast)) {
      strAtomNames.erase(iLast, 1);
      strAtomNames.erase(0, 1);
      restraintType = NOE_AND;
    }
    // Restraint type: NOE_OR
    // Syntax: at1,at2,at3...
    // Check for no brackets present
    else if ((osqb == std::string::npos) && (csqb == std::string::npos)) {
      restraintType = NOE_OR;
    }
  }
  return restraintType;
}
