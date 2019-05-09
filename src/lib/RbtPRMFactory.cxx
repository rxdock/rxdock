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

#include "RbtPRMFactory.h"
#include "RbtChromPositionRefData.h"
#include "RbtCrdFileSource.h"
#include "RbtLigandFlexData.h"
#include "RbtMOL2FileSource.h"
#include "RbtMdlFileSource.h"
#include "RbtParameterFileSource.h"
#include "RbtPdbFileSource.h"
#include "RbtPsfFileSource.h"
#include "RbtReceptorFlexData.h"
#include "RbtSolventFlexData.h"

const std::string &RbtPRMFactory::_CT = "RbtPRMFactory";
const std::string &RbtPRMFactory::_REC_SECTION = "";
const std::string &RbtPRMFactory::_REC_FILE = "RECEPTOR_FILE";
const std::string &RbtPRMFactory::_REC_TOPOL_FILE = "RECEPTOR_TOPOL_FILE";
const std::string &RbtPRMFactory::_REC_COORD_FILE = "RECEPTOR_COORD_FILE";
const std::string &RbtPRMFactory::_REC_NUM_COORD_FILES =
    "RECEPTOR_NUM_COORD_FILES";
const std::string &RbtPRMFactory::_REC_FLEX_DISTANCE = "RECEPTOR_FLEX";
const std::string &RbtPRMFactory::_REC_DIHEDRAL_STEP = "RECEPTOR_DIHEDRAL_STEP";
const std::string &RbtPRMFactory::_LIG_SECTION = "LIGAND";
const std::string &RbtPRMFactory::_SOLV_SECTION = "SOLVENT";
const std::string &RbtPRMFactory::_SOLV_FILE = "FILE";

RbtPRMFactory::RbtPRMFactory(RbtParameterFileSource *pParamSource)
    : m_pParamSource(pParamSource), m_pDS(NULL), m_iTrace(0) {}

RbtPRMFactory::RbtPRMFactory(RbtParameterFileSource *pParamSource,
                             RbtDockingSite *pDS)
    : m_pParamSource(pParamSource), m_pDS(pDS), m_iTrace(0) {}
RbtModelPtr RbtPRMFactory::CreateReceptor() {
  RbtModelPtr retVal;
  m_pParamSource->SetSection(_REC_SECTION);
  // Detect if we have an ensemble of receptor coordinate files defined
  bool bEnsemble =
      m_pParamSource->isParameterPresent(_REC_NUM_COORD_FILES) &&
      (m_pParamSource->GetParameterValue(_REC_NUM_COORD_FILES) > 0);

  // Read topology and coordinates from a single molecular file source
  if (m_pParamSource->isParameterPresent(_REC_FILE)) {
    if (m_iTrace > 0) {
      std::cout << std::endl
                << "Using " << _REC_FILE
                << " as combined source of topology and 3D coordinates"
                << std::endl
                << std::endl;
    }
    std::string strFile = m_pParamSource->GetParameterValueAsString(_REC_FILE);
    RbtMolecularFileSourcePtr theSource = CreateMolFileSource(strFile);
    bool isOK = theSource->isAtomListSupported() &&
                theSource->isBondListSupported() &&
                theSource->isCoordinatesSupported();
    if (!isOK) {
      if (m_iTrace > 0) {
        std::cout << std::endl
                  << "Incompatible file type for " << _REC_FILE << " option"
                  << std::endl;
        std::cout
            << "File type must provide atom list, bond list and coordinates "
               "in a single file"
            << std::endl;
        std::cout << "Consider using " << _REC_TOPOL_FILE << " and "
                  << _REC_COORD_FILE
                  << " to specify topology (atoms/bonds) and coordinate files "
                     "separately"
                  << std::endl
                  << std::endl;
      }
      throw RbtBadReceptorFile(_WHERE_,
                               "Inappropriate molecular file type for " +
                                   _REC_FILE + " option");
    }
    retVal = new RbtModel(theSource);
  }
  // Read topology and coordinates separately from two molecular file sources
  else {
    if (m_iTrace > 0) {
      std::cout << std::endl
                << "Using " << _REC_TOPOL_FILE << " as topology source"
                << std::endl
                << std::endl;
    }
    std::string strTopolFile =
        m_pParamSource->GetParameterValueAsString(_REC_TOPOL_FILE);
    RbtMolecularFileSourcePtr theTopolSource =
        CreateMolFileSource(strTopolFile);
    bool isOK = theTopolSource->isAtomListSupported() &&
                theTopolSource->isBondListSupported();
    if (!isOK) {
      if (m_iTrace > 0) {
        std::cout << std::endl
                  << "Incompatible file type for " << _REC_TOPOL_FILE
                  << " option" << std::endl;
        std::cout
            << "File type must provide topology information (atom list and "
               "bond list)"
            << std::endl
            << std::endl;
      }
      throw RbtBadReceptorFile(_WHERE_,
                               "Inappropriate molecular file type for " +
                                   _REC_TOPOL_FILE + " option");
    }
    retVal = new RbtModel(theTopolSource);
    if (!bEnsemble) {
      // Now read coord file
      if (m_iTrace > 0) {
        std::cout << std::endl
                  << "Using " << _REC_COORD_FILE
                  << " as source of 3D coordinates" << std::endl
                  << std::endl;
      }
      std::string strCoordFile =
          m_pParamSource->GetParameterValueAsString(_REC_COORD_FILE);
      RbtMolecularFileSourcePtr theCoordSource =
          CreateMolFileSource(strCoordFile);
      isOK = theCoordSource->isAtomListSupported() &&
             theCoordSource->isCoordinatesSupported();
      if (!isOK) {
        if (m_iTrace > 0) {
          std::cout << std::endl
                    << "Incompatible file type for " << _REC_COORD_FILE
                    << " option" << std::endl;
          std::cout
              << "File type must provide coordinate information (atom list "
                 "with 3D coords)"
              << std::endl
              << std::endl;
        }
        throw RbtBadReceptorFile(_WHERE_,
                                 "Inappropriate molecular file type for " +
                                     _REC_COORD_FILE + " option");
      }
      retVal->UpdateCoords(theCoordSource);
    }
  }

  // Optional read of multiple receptor conformations (numbered from 1 thru N)
  if (bEnsemble) {
    int n = m_pParamSource->GetParameterValue(_REC_NUM_COORD_FILES);
    if (m_iTrace > 0) {
      std::cout << std::endl
                << "Using ensemble of " << _REC_COORD_FILE
                << "'s as source of 3D coordinates (N=" << n << ")" << std::endl
                << std::endl;
    }
    for (int i = 1; i <= n; i++) {
      std::ostringstream ostr;
      ostr << _REC_COORD_FILE << "_" << i;
      std::string paramName(ostr.str());
      std::string strCoordFile =
          m_pParamSource->GetParameterValueAsString(paramName);
      if (m_iTrace > 0) {
        std::cout << std::endl << "I=" << i << std::endl;
      }
      RbtMolecularFileSourcePtr theCoordSource =
          CreateMolFileSource(strCoordFile);
      bool isOK = theCoordSource->isAtomListSupported() &&
                  theCoordSource->isCoordinatesSupported();
      if (!isOK) {
        if (m_iTrace > 0) {
          std::cout << std::endl
                    << "Incompatible file type for " << paramName << " option"
                    << std::endl;
          std::cout
              << "File type must provide coordinate information (atom list "
                 "with 3D coords)"
              << std::endl
              << std::endl;
        }
        throw RbtBadReceptorFile(_WHERE_,
                                 "Inappropriate molecular file type for " +
                                     paramName + " option");
      }
      retVal->UpdateCoords(theCoordSource);
      retVal->SaveCoords(strCoordFile);
    }
    int nCoords = retVal->GetNumSavedCoords() - 1;
    std::cout << "Total number of receptor conformations read = " << nCoords
              << std::endl;
  }

  // If the docking site is defined, then we can define the
  // receptor flexibility
  if (m_pDS) {
    AttachReceptorFlexData(retVal);
  }

  return retVal;
}

RbtModelPtr RbtPRMFactory::CreateLigand(RbtBaseMolecularFileSource *pSource) {
  RbtModelPtr retVal;
  // TODO: Move some of the status checks from rbdock to here.
  retVal = new RbtModel(pSource);
  // If the docking site is defined, then we can define the
  // ligand flexibility
  if (m_pDS) {
    AttachLigandFlexData(retVal);
  }
  return retVal;
}

RbtModelList RbtPRMFactory::CreateSolvent() {
  RbtModelList retVal;
  m_pParamSource->SetSection(_SOLV_SECTION);
  if (m_pParamSource->isParameterPresent(_SOLV_FILE)) {
    std::string strFile = m_pParamSource->GetParameterValueAsString(_SOLV_FILE);
    if (m_iTrace > 0) {
      std::cout << std::endl
                << "Reading solvent from " << strFile << std::endl
                << std::endl;
    }
    RbtMolecularFileSourcePtr theSource = CreateMolFileSource(strFile);
    bool isOK =
        theSource->isAtomListSupported() && theSource->isCoordinatesSupported();
    if (!isOK) {
      if (m_iTrace > 0) {
        std::cout << std::endl
                  << "Incompatible file type for " << _SOLV_SECTION
                  << "::" << _SOLV_FILE << " option" << std::endl;
        std::cout
            << "File type must provide atom list and coordinates in a single "
               "file"
            << std::endl;
      }
      throw RbtBadReceptorFile(
          _WHERE_, "Inappropriate molecular file type for " + _SOLV_SECTION +
                       "::" + _SOLV_FILE + " option");
    }
    // Find all the water molecules within the file source and convert each to a
    // separate RbtModel
    Rbt::isAtomName_eq bIsO("OW");
    Rbt::isAtomName_eq bIsH1("H1");
    Rbt::isAtomName_eq bIsH2("H2");
    Rbt::isCoordinationNumber_eq bIsIsolated(0);
    RbtTriposAtomType triposType; // Tripos atom typer
    int nAtomId(0);               // incremental counter of all atoms created
    int nBondId(0);               // incremental counter of all bonds created
    RbtAtomList oAtomList = Rbt::GetAtomList(theSource->GetAtomList(), bIsO);
    RbtAtomList h1AtomList = Rbt::GetAtomList(theSource->GetAtomList(), bIsH1);
    RbtAtomList h2AtomList = Rbt::GetAtomList(theSource->GetAtomList(), bIsH2);
    for (RbtAtomListConstIter iter = oAtomList.begin(); iter != oAtomList.end();
         ++iter) {
      RbtAtomPtr oAtom = (*iter);
      // std::cout << "Processing " << (*oAtom) << "..." << std::endl;
      // Isolated oxygen atom:
      // Need to find its matching hydrogens in the file and create the bonds
      if (bIsIsolated(oAtom)) {
        // std::cout << "...is isolated" << std::endl;
        Rbt::isSubunitId_eq bInSameSubunit(oAtom->GetSubunitId());
        RbtAtomList matchingH1AtomList =
            Rbt::GetAtomList(h1AtomList, bInSameSubunit);
        RbtAtomList matchingH2AtomList =
            Rbt::GetAtomList(h2AtomList, bInSameSubunit);
        // std::cout << "..." << matchingH1AtomList.size() << " matching H1" <<
        // std::endl; std::cout << "..." << matchingH2AtomList.size() << "
        // matching H2" << std::endl; Found a water!
        if ((matchingH1AtomList.size() == 1) &&
            (matchingH2AtomList.size() == 1)) {
          RbtAtomPtr h1Atom = matchingH1AtomList.front();
          RbtAtomPtr h2Atom = matchingH2AtomList.front();
          // define atomic attributes
          oAtom->SetAtomId(++nAtomId);
          oAtom->SetAtomicNo(8);
          oAtom->SetAtomicMass(16.000);
          oAtom->SetHybridState(RbtAtom::SP3);
          oAtom->SetVdwRadius(1.35);
          h1Atom->SetAtomId(++nAtomId);
          h1Atom->SetAtomicNo(1);
          h1Atom->SetAtomicMass(1.008);
          h1Atom->SetHybridState(RbtAtom::SP3);
          h1Atom->SetVdwRadius(0.5);
          h2Atom->SetAtomId(++nAtomId);
          h2Atom->SetAtomicNo(1);
          h2Atom->SetAtomicMass(1.008);
          h2Atom->SetHybridState(RbtAtom::SP3);
          h2Atom->SetVdwRadius(0.5);
          // Create the atom list and bond list for this water
          RbtAtomList waterAtomList;
          RbtBondList waterBondList;
          waterAtomList.push_back(oAtom);
          waterAtomList.push_back(h1Atom);
          waterAtomList.push_back(h2Atom);
          waterBondList.push_back(new RbtBond(++nBondId, oAtom, h1Atom, 1));
          waterBondList.push_back(new RbtBond(++nBondId, oAtom, h2Atom, 1));
          // Now we've created the bonds we can assign the Tripos types
          // automatically
          oAtom->SetTriposType(triposType(oAtom, true));
          h1Atom->SetTriposType(triposType(h1Atom, true));
          h2Atom->SetTriposType(triposType(h2Atom, true));
          if (m_iTrace > 0) {
            std::cout << "Creating water model # " << retVal.size()
                      << std::endl;
            std::cout << (*oAtom) << std::endl
                      << (*h1Atom) << std::endl
                      << (*h2Atom) << std::endl;
          }
          RbtModelPtr solvent(new RbtModel(waterAtomList, waterBondList));
          // If the docking site is defined, then we can define the
          // solvent flexibility
          if (m_pDS) {
            AttachSolventFlexData(solvent);
          }
          retVal.push_back(solvent);
        }
      }
    }
  }

  return retVal;
}

void RbtPRMFactory::AttachReceptorFlexData(RbtModel *pReceptor) {
  m_pParamSource->SetSection(_REC_SECTION);
  // Check whether flexible receptor is requested (terminal OH/NH3)
  // Parameter value is the range from the docking volume to include
  if (m_pParamSource->isParameterPresent(_REC_FLEX_DISTANCE)) {
    double flexDist = m_pParamSource->GetParameterValue(_REC_FLEX_DISTANCE);
    if (m_iTrace > 0) {
      std::cout << std::endl
                << "Target OH/NH3 groups within " << flexDist
                << " A of docking site requested as flexible" << std::endl;
    }
    RbtFlexData *pFlexData = new RbtReceptorFlexData(m_pDS);
    pFlexData->SetParameter(RbtReceptorFlexData::_FLEX_DISTANCE, flexDist);
    if (m_pParamSource->isParameterPresent(_REC_DIHEDRAL_STEP)) {
      double dihedralStepSize =
          m_pParamSource->GetParameterValue(_REC_DIHEDRAL_STEP);
      pFlexData->SetParameter(RbtReceptorFlexData::_DIHEDRAL_STEP,
                              dihedralStepSize);
      if (m_iTrace > 0) {
        std::cout << std::endl
                  << _REC_DIHEDRAL_STEP << " = " << dihedralStepSize
                  << std::endl;
      }
    }
    pReceptor->SetFlexData(pFlexData);
    if (m_iTrace > 0) {
      std::cout << std::endl
                << "RECEPTOR FLEXIBILITY PARAMETERS:" << std::endl
                << *pFlexData << std::endl;
    }
  }
}

void RbtPRMFactory::AttachLigandFlexData(RbtModel *pLigand) {
  m_pParamSource->SetSection(_LIG_SECTION);
  RbtStringList paramList = m_pParamSource->GetParameterList();
  RbtFlexData *pFlexData = new RbtLigandFlexData(m_pDS);
  for (RbtStringListConstIter iter = paramList.begin(); iter != paramList.end();
       ++iter) {
    std::string strValue = m_pParamSource->GetParameterValueAsString(*iter);
    pFlexData->SetParameter(*iter, strValue);
  }
  pLigand->SetFlexData(pFlexData);
  if (m_iTrace > 0) {
    std::cout << std::endl
              << "LIGAND FLEXIBILITY PARAMETERS:" << std::endl
              << *pFlexData << std::endl;
  }
}

void RbtPRMFactory::AttachSolventFlexData(RbtModel *pSolvent) {
  m_pParamSource->SetSection(_SOLV_SECTION);
  RbtStringList paramList = m_pParamSource->GetParameterList();
  RbtFlexData *pFlexData = new RbtSolventFlexData(m_pDS);
  // Determine the occupancy and flexibility mode on a per-solvent basis
  // from the User1 and User2 values of the first atom (oxygen)
  // User1 = occupancy [0,1]
  // User2 = mode [0,8]
  // mode 0 => TRANS=FIXED    ROT=FIXED
  // mode 1 => TRANS=FIXED    ROT=TETHERED
  // mode 2 => TRANS=FIXED    ROT=FREE
  // mode 3 => TRANS=TETHERED ROT=FIXED
  // mode 4 => TRANS=TETHERED ROT=TETHERED
  // mode 5 => TRANS=TETHERED ROT=FREE
  // mode 6 => TRANS=FREE     ROT=FIXED
  // mode 7 => TRANS=FREE     ROT=TETHERED
  // mode 8 => TRANS=FREE     ROT=FREE
  RbtAtomList atomList = pSolvent->GetAtomList();
  double occupancy(1.0);
  RbtChromElement::eMode transMode(RbtChromElement::FIXED);
  RbtChromElement::eMode rotMode(RbtChromElement::FIXED);
  if (!atomList.empty()) {
    RbtAtom *pFirstAtom = atomList[0];
    occupancy = pFirstAtom->GetUser1Value();
    occupancy = std::min(1.0, occupancy);
    occupancy = std::max(0.0, occupancy);
    int mode = (int)pFirstAtom->GetUser2Value();
    mode = std::min(8, mode);
    mode = std::max(0, mode);
    transMode = (RbtChromElement::eMode)(mode / 3);
    rotMode = (RbtChromElement::eMode)(mode % 3);
  }
  pFlexData->SetParameter(RbtLigandFlexData::_TRANS_MODE,
                          RbtChromElement::ModeToStr(transMode));
  pFlexData->SetParameter(RbtLigandFlexData::_ROT_MODE,
                          RbtChromElement::ModeToStr(rotMode));
  pFlexData->SetParameter(RbtSolventFlexData::_OCCUPANCY, occupancy);
  // Set the remaining flexibility parameters from the <receptor>.prm file
  // source
  for (RbtStringListConstIter iter = paramList.begin(); iter != paramList.end();
       ++iter) {
    if ((*iter) != _SOLV_FILE) {
      std::string strValue = m_pParamSource->GetParameterValueAsString(*iter);
      pFlexData->SetParameter(*iter, strValue);
    }
  }
  pSolvent->SetFlexData(pFlexData);
  if (m_iTrace > 0) {
    std::cout << std::endl
              << "SOLVENT FLEXIBILITY PARAMETERS:" << std::endl
              << *pFlexData << std::endl;
  }
}

RbtMolecularFileSourcePtr
RbtPRMFactory::CreateMolFileSource(const std::string &fileName) {
  RbtMolecularFileSourcePtr retVal;
  std::string fileType = Rbt::GetFileType(fileName);
  std::string fileTypeUpper;
  std::transform(fileType.begin(), fileType.end(),
                 std::back_inserter(fileTypeUpper), toupper);
  std::string fullFileName = Rbt::GetRbtFileName("", fileName);

  if (m_iTrace > 0) {
    std::cout << _CT << ": File name requested = " << fileName
              << "; type = " << fileTypeUpper << std::endl;
    std::cout << _CT << ": Reading file from " << fullFileName << std::endl;
  }

  bool bImplH = true;
  if (m_pParamSource->isParameterPresent("RECEPTOR_ALL_H")) {
    RbtVariant vAllH(
        m_pParamSource->GetParameterValueAsString("RECEPTOR_ALL_H"));
    bImplH = !vAllH.Bool();
  }

  if (fileTypeUpper == "MOL2") {
    retVal = new RbtMOL2FileSource(fullFileName, bImplH);
  } else if (fileTypeUpper == "PSF") {
    std::string strMassesFile =
        m_pParamSource->GetParameterValueAsString("RECEPTOR_MASSES_FILE");
    strMassesFile = Rbt::GetRbtFileName("data", strMassesFile);
    if (m_iTrace > 0) {
      std::cout << _CT << ": Using file " << strMassesFile
                << " to lookup Charmm atom type info" << std::endl;
    }
    retVal = new RbtPsfFileSource(fullFileName, strMassesFile, bImplH);
  } else if (fileTypeUpper == "CRD") {
    retVal = new RbtCrdFileSource(fullFileName);
    bImplH = false;
  } else if (fileTypeUpper == "PDB") {
    retVal = new RbtPdbFileSource(fullFileName);
    bImplH = false;
  } else if ((fileTypeUpper == "SD") || (fileTypeUpper == "SDF")) {
    retVal = new RbtMdlFileSource(fullFileName, false, false, bImplH);
  } else {
    throw RbtBadReceptorFile(
        _WHERE_, fileType + " is not a supported molecular file type");
  }

  if (m_pParamSource->isParameterPresent("RECEPTOR_SEGMENT_NAME")) {
    std::string strSegmentName =
        m_pParamSource->GetParameterValueAsString("RECEPTOR_SEGMENT_NAME");
    RbtSegmentMap segmentMap = Rbt::ConvertStringToSegmentMap(strSegmentName);
    retVal->SetSegmentFilterMap(segmentMap);
    if (m_iTrace > 0) {
      std::cout << _CT << ": Setting segment/chain filter to " << strSegmentName
                << std::endl;
    }
  } else {
    if (m_iTrace > 0) {
      std::cout << _CT << ": No segment/chain filter defined" << std::endl;
    }
  }

  if (m_iTrace > 0) {
    if (bImplH)
      std::cout << _CT << ": Removing non-polar hydrogens" << std::endl;
    else
      std::cout << _CT << ": Reading all hydrogens present" << std::endl;
  }

  RbtError status = retVal->Status();
  if (!status.isOK())
    throw status;
  int nAtoms = retVal->GetNumAtoms();
  if (nAtoms == 0) {
    if (m_iTrace > 0) {
      std::cout << _CT << ": File source contains zero atoms!" << std::endl;
      std::cout << _CT << ": Have you defined RECEPTOR_SEGMENT_NAME correctly?"
                << std::endl;
    }
    throw RbtBadReceptorFile(_WHERE_,
                             "Zero atoms provided by " + retVal->GetFileName());
  } else {
    if (m_iTrace > 0) {
      std::cout << _CT << ": File source contains " << nAtoms << " atoms"
                << std::endl;
    }
  }
  return retVal;
}
