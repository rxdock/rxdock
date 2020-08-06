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

#include "rxdock/TransformFactory.h"
// Component transforms
#include "rxdock/AlignTransform.h"
#include "rxdock/GATransform.h"
#include "rxdock/NullTransform.h"
#include "rxdock/RandLigTransform.h"
#include "rxdock/RandPopTransform.h"
#include "rxdock/SimAnnTransform.h"
#include "rxdock/SimplexTransform.h"

#include "rxdock/FileError.h"
#include "rxdock/SFRequest.h"

using namespace rxdock;

// Parameter name which identifies a scoring function definition
std::string TransformFactory::_TRANSFORM("TRANSFORM");

////////////////////////////////////////
// Constructors/destructors
TransformFactory::TransformFactory() {}

TransformFactory::~TransformFactory() {}

////////////////////////////////////////
// Public methods
////////////////
// Creates a single transform object of type strTransformClass, and name strName
// e.g. strTransformClass = SimAnnTransform
BaseTransform *TransformFactory::Create(const std::string &strTransformClass,
                                        const std::string &strName) {
  // Component transforms
  if (strTransformClass == SimAnnTransform::_CT)
    return new SimAnnTransform(strName);
  if (strTransformClass == GATransform::_CT)
    return new GATransform(strName);
  if (strTransformClass == AlignTransform::_CT)
    return new AlignTransform(strName);
  if (strTransformClass == NullTransform::_CT)
    return new NullTransform(strName);
  if (strTransformClass == RandLigTransform::_CT)
    return new RandLigTransform(strName);
  if (strTransformClass == RandPopTransform::_CT)
    return new RandPopTransform(strName);
  if (strTransformClass == SimplexTransform::_CT)
    return new SimplexTransform(strName);
  // Aggregate transforms
  if (strTransformClass == TransformAgg::_CT)
    return new TransformAgg(strName);

  throw BadArgument(_WHERE_, "Unknown transform " + strTransformClass);
}

// Creates an aggregate transform from a parameter file source
// Each component transform is in a named section, which should minimally
// contain a TRANSFORM parameter whose value is the class name to instantiate
// strTransformClasses contains a comma-delimited list of transform class names
// to instantiate If strTransformClasses is empty, all named sections in
// spPrmSource are scanned for valid transform definitions Transform parameters
// and scoring function requests are set from the list of parameters in each
// named section
TransformAgg *
TransformFactory::CreateAggFromFile(ParameterFileSourcePtr spPrmSource,
                                    const std::string &strName,
                                    const std::string &strTransformClasses) {
  // Get list of transform objects to create
  std::vector<std::string> transformList =
      ConvertDelimitedStringToList(strTransformClasses);
  // If strTransformClasses is empty, then default to reading all sections of
  // the parameter file for valid transform definitions In this case we do not
  // throw an error if a particular section is not a transform, we simply skip
  // it
  bool bThrowError(true);
  if (transformList.empty()) {
    transformList = spPrmSource->GetSectionList();
    bThrowError = false;
  }

  // Create empty aggregate
  TransformAgg *pTransformAgg(new TransformAgg(strName));

  for (std::vector<std::string>::const_iterator tIter = transformList.begin();
       tIter != transformList.end(); tIter++) {
    spPrmSource->SetSection(*tIter);
    // Check if this section is a valid scoring function definition
    if (spPrmSource->isParameterPresent(_TRANSFORM)) {
      std::string strTransformClass(
          spPrmSource->GetParameterValueAsString(_TRANSFORM));
      // Create new transform according to the string value of _TRANSFORM
      // parameter
      BaseTransform *pTransform = Create(strTransformClass, *tIter);
      // Set all the transform parameters from the rest of the parameters listed
      std::vector<std::string> prmList = spPrmSource->GetParameterList();
      for (std::vector<std::string>::const_iterator prmIter = prmList.begin();
           prmIter != prmList.end(); prmIter++) {
        // Look for scoring function request (PARAM@SF)
        // Only SetParamRequest currently supported
        std::vector<std::string> compList =
            ConvertDelimitedStringToList(*prmIter, "@");
        if (compList.size() == 2) {
          RequestPtr spReq(new SFSetParamRequest(
              compList[1], compList[0],
              spPrmSource->GetParameterValueAsString(*prmIter)));
          pTransform->AddSFRequest(spReq);
        } else if ((*prmIter) != _TRANSFORM) { // Skip _TRANSFORM parameter
          pTransform->SetParameter(
              *prmIter, spPrmSource->GetParameterValueAsString(*prmIter));
        }
      }
      pTransformAgg->Add(pTransform);
    } else if (bThrowError) {
      throw FileMissingParameter(_WHERE_, "Missing " + _TRANSFORM +
                                              " parameter in section " +
                                              (*tIter));
    }
  }
  return pTransformAgg;
}
