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

#include "rxdock/SFFactory.h"
// Precalculated-grid scoring functions
#include "rxdock/CavityFillSF.h"
#include "rxdock/CavityGridSF.h"
#include "rxdock/VdwGridSF.h"

// Indexed-grid scoring functions
#include "rxdock/AromIdxSF.h"
#include "rxdock/PolarIdxSF.h"
#include "rxdock/VdwIdxSF.h"

// Intramolecular scoring functions
#include "rxdock/DihedralIntraSF.h"
#include "rxdock/DihedralTargetSF.h"
#include "rxdock/PolarIntraSF.h"
#include "rxdock/VdwIntraSF.h"

#include "rxdock/NmrSF.h"

// PMF scoring
#include "rxdock/PMFGridSF.h"
#include "rxdock/PMFIdxSF.h"
#include "rxdock/SetupPMFSF.h"

// SA scoring
#include "rxdock/SAIdxSF.h"
#include "rxdock/SetupSASF.h"

// Tethered Atoms scoring
#include "rxdock/TetherSF.h"

// Pharmacophore scoring
#include "rxdock/PharmaSF.h"

// Misc scoring functions
#include "rxdock/ConstSF.h"
#include "rxdock/RotSF.h"

// Pseudo-scoring functions
#include "rxdock/SetupPolarSF.h"

#include "rxdock/FileError.h"

using namespace rxdock;

// Parameter name which identifies a scoring function definition
std::string SFFactory::_SF("SCORING_FUNCTION");

////////////////////////////////////////
// Constructors/destructors
SFFactory::SFFactory() {}

SFFactory::~SFFactory() {}

////////////////////////////////////////
// Public methods
////////////////
BaseSF *SFFactory::Create(const std::string &_strSFClass,
                          const std::string &strName) {
  std::string strSFClass = _strSFClass;
  // compatibility with old format Rbt*
  if (strSFClass.substr(0, 3) == "Rbt")
    strSFClass = _strSFClass.substr(3, std::string::npos);

  // Precalculated-grid scoring functions
  if (strSFClass == VdwGridSF::_CT)
    return new VdwGridSF(strName);
  if (strSFClass == CavityGridSF::_CT)
    return new CavityGridSF(strName);
  if (strSFClass == CavityFillSF::_CT)
    return new CavityFillSF(strName);

  // Indexed-grid scoring functions
  if (strSFClass == AromIdxSF::_CT)
    return new AromIdxSF(strName);
  if (strSFClass == PolarIdxSF::_CT)
    return new PolarIdxSF(strName);
  if (strSFClass == VdwIdxSF::_CT)
    return new VdwIdxSF(strName);

  // Intramolecular scoring functions
  if (strSFClass == VdwIntraSF::_CT)
    return new VdwIntraSF(strName);
  if (strSFClass == PolarIntraSF::_CT)
    return new PolarIntraSF(strName);
  if (strSFClass == DihedralIntraSF::_CT)
    return new DihedralIntraSF(strName);
  if (strSFClass == DihedralTargetSF::_CT)
    return new DihedralTargetSF(strName);

  // Restraint scoring functions
  if (strSFClass == NmrSF::_CT)
    return new NmrSF(strName);
  if (strSFClass == "NmrIntnSF")
    return new NmrSF(strName);

  // Tethered Atoms scoring
  if (strSFClass == TetherSF::_CT)
    return new TetherSF(strName);

  // Pharmacophore scoring
  if (strSFClass == PharmaSF::_CT)
    return new PharmaSF(strName);

  // Potential Mean Force (PMF) scoring functions
  if (strSFClass == SetupPMFSF::_CT)
    return new SetupPMFSF(strName);
  if (strSFClass == PMFIdxSF::_CT)
    return new PMFIdxSF(strName);
  if (strSFClass == PMFGridSF::_CT)
    return new PMFGridSF(strName);

  // Solvent Accessible SF
  if (strSFClass == SetupSASF::_CT)
    return new SetupSASF(strName);
  if (strSFClass == SAIdxSF::_CT)
    return new SAIdxSF(strName);

  // Misc scoring functions
  if (strSFClass == ConstSF::_CT)
    return new ConstSF(strName);
  if (strSFClass == RotSF::_CT)
    return new RotSF(strName);

  // Pseudo-scoring functions
  if (strSFClass == SetupPolarSF::_CT)
    return new SetupPolarSF(strName);

  // Aggregate scoring functions
  if (strSFClass == SFAgg::_CT)
    return new SFAgg(strName);

  throw BadArgument(_WHERE_, "Unknown scoring function " + strSFClass);
}

// Creates an aggregate scoring function from a parameter file source
// Each component SF is in a named section, which should minimally contain a
// SCORING_FUNCTION parameter whose value is the class name to instantiate
// strSFClasses contains a comma-delimited list of SF class names to instantiate
// If strSFClasses is empty, all named sections in spPrmSource are scanned for
// valid scoring function definitions SF parameters are set from the list of
// parameters in each named section
SFAgg *SFFactory::CreateAggFromFile(ParameterFileSourcePtr spPrmSource,
                                    const std::string &strName,
                                    const std::string &strSFClasses) {
  // Get list of scoring function objects to create
  std::vector<std::string> sfList = ConvertDelimitedStringToList(strSFClasses);
  // If strSFClasses is empty, then default to reading all sections of the
  // parameter file for valid scoring function definitions
  // In this case we do not throw an error if a particular section
  // is not a scoring function, we simply skip it
  bool bThrowError(true);
  if (sfList.empty()) {
    sfList = spPrmSource->GetSectionList();
    bThrowError = false;
  }

  // Create empty aggregate
  SFAgg *pSFAgg(new SFAgg(strName));

  for (std::vector<std::string>::const_iterator sfIter = sfList.begin();
       sfIter != sfList.end(); sfIter++) {
    spPrmSource->SetSection(*sfIter);
    // Check if this section is a valid scoring function definition
    if (spPrmSource->isParameterPresent(_SF)) {
      std::string strSFClass(spPrmSource->GetParameterValueAsString(_SF));
      // Create new scoring function according to the string value of _SF
      // parameter
      BaseSF *pSF = Create(strSFClass, *sfIter);
      // Set all the scoring function parameters from the rest of the parameters
      // listed
      std::vector<std::string> prmList = spPrmSource->GetParameterList();
      for (std::vector<std::string>::const_iterator prmIter = prmList.begin();
           prmIter != prmList.end(); prmIter++) {
        if ((*prmIter) != _SF) { // Skip _SF parameter
          pSF->SetParameter(*prmIter,
                            spPrmSource->GetParameterValueAsString(*prmIter));
        }
      }
      pSFAgg->Add(pSF);
    } else if (bThrowError) {
      throw FileMissingParameter(
          _WHERE_, "Missing " + _SF + " parameter in section " + (*sfIter));
    }
  }
  return pSFAgg;
}
