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

#include "PRMFactory.h"
#include "ChromPositionRefData.h"
#include "CrdFileSource.h"
#include "LigandFlexData.h"
#include "MOL2FileSource.h"
#include "MdlFileSource.h"
#include "ParameterFileSource.h"
#include "PdbFileSource.h"
#include "PsfFileSource.h"
#include "ReceptorFlexData.h"
#include "SolventFlexData.h"

using namespace rxdock;

const std::string &PRMFactory::_CT = "PRMFactory";
const std::string &PRMFactory::_REC_SECTION = "";
const std::string &PRMFactory::_REC_FILE = "RECEPTOR_FILE";
const std::string &PRMFactory::_REC_TOPOL_FILE = "RECEPTOR_TOPOL_FILE";
const std::string &PRMFactory::_REC_COORD_FILE = "RECEPTOR_COORD_FILE";
const std::string &PRMFactory::_REC_NUM_COORD_FILES =
    "RECEPTOR_NUM_COORD_FILES";
const std::string &PRMFactory::_REC_FLEX_DISTANCE = "RECEPTOR_FLEX";
const std::string &PRMFactory::_REC_DIHEDRAL_STEP = "RECEPTOR_DIHEDRAL_STEP";
const std::string &PRMFactory::_LIG_SECTION = "LIGAND";
const std::string &PRMFactory::_SOLV_SECTION = "SOLVENT";
const std::string &PRMFactory::_SOLV_FILE = "FILE";

PRMFactory::PRMFactory(ParameterFileSource *pParamSource)
    : m_pParamSource(pParamSource), m_pDS(nullptr), m_iTrace(0) {}

PRMFactory::PRMFactory(ParameterFileSource *pParamSource, DockingSite *pDS)
    : m_pParamSource(pParamSource), m_pDS(pDS), m_iTrace(0) {}
ModelPtr PRMFactory::CreateReceptor() {
  ModelPtr retVal;
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
    MolecularFileSourcePtr theSource = CreateMolFileSource(strFile);
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
      throw BadReceptorFile(_WHERE_, "Inappropriate molecular file type for " +
                                         _REC_FILE + " option");
    }
    retVal = new Model(theSource);
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
    MolecularFileSourcePtr theTopolSource = CreateMolFileSource(strTopolFile);
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
      throw BadReceptorFile(_WHERE_, "Inappropriate molecular file type for " +
                                         _REC_TOPOL_FILE + " option");
    }
    retVal = new Model(theTopolSource);
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
      MolecularFileSourcePtr theCoordSource = CreateMolFileSource(strCoordFile);
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
        throw BadReceptorFile(_WHERE_,
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
      MolecularFileSourcePtr theCoordSource = CreateMolFileSource(strCoordFile);
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
        throw BadReceptorFile(_WHERE_,
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

ModelPtr PRMFactory::CreateLigand(BaseMolecularFileSource *pSource) {
  ModelPtr retVal;
  // TODO: Move some of the status checks from rbdock to here.
  retVal = new Model(pSource);
  // If the docking site is defined, then we can define the
  // ligand flexibility
  if (m_pDS) {
    AttachLigandFlexData(retVal);
  }
  return retVal;
}

ModelList PRMFactory::CreateSolvent() {
  ModelList retVal;
  m_pParamSource->SetSection(_SOLV_SECTION);
  if (m_pParamSource->isParameterPresent(_SOLV_FILE)) {
    std::string strFile = m_pParamSource->GetParameterValueAsString(_SOLV_FILE);
    if (m_iTrace > 0) {
      std::cout << std::endl
                << "Reading solvent from " << strFile << std::endl
                << std::endl;
    }
    MolecularFileSourcePtr theSource = CreateMolFileSource(strFile);
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
      throw BadReceptorFile(_WHERE_, "Inappropriate molecular file type for " +
                                         _SOLV_SECTION + "::" + _SOLV_FILE +
                                         " option");
    }
    // Find all the water molecules within the file source and convert each to a
    // separate Model
    isAtomName_eq bIsO("OW");
    isAtomName_eq bIsH1("H1");
    isAtomName_eq bIsH2("H2");
    isCoordinationNumber_eq bIsIsolated(0);
    TriposAtomType triposType; // Tripos atom typer
    int nAtomId(0);            // incremental counter of all atoms created
    int nBondId(0);            // incremental counter of all bonds created
    AtomList oAtomList =
        GetAtomListWithPredicate(theSource->GetAtomList(), bIsO);
    AtomList h1AtomList =
        GetAtomListWithPredicate(theSource->GetAtomList(), bIsH1);
    AtomList h2AtomList =
        GetAtomListWithPredicate(theSource->GetAtomList(), bIsH2);
    for (AtomListConstIter iter = oAtomList.begin(); iter != oAtomList.end();
         ++iter) {
      AtomPtr oAtom = (*iter);
      // std::cout << "Processing " << (*oAtom) << "..." << std::endl;
      // Isolated oxygen atom:
      // Need to find its matching hydrogens in the file and create the bonds
      if (bIsIsolated(oAtom)) {
        // std::cout << "...is isolated" << std::endl;
        isSubunitId_eq bInSameSubunit(oAtom->GetSubunitId());
        AtomList matchingH1AtomList =
            GetAtomListWithPredicate(h1AtomList, bInSameSubunit);
        AtomList matchingH2AtomList =
            GetAtomListWithPredicate(h2AtomList, bInSameSubunit);
        // std::cout << "..." << matchingH1AtomList.size() << " matching H1" <<
        // std::endl; std::cout << "..." << matchingH2AtomList.size() << "
        // matching H2" << std::endl; Found a water!
        if ((matchingH1AtomList.size() == 1) &&
            (matchingH2AtomList.size() == 1)) {
          AtomPtr h1Atom = matchingH1AtomList.front();
          AtomPtr h2Atom = matchingH2AtomList.front();
          // define atomic attributes
          oAtom->SetAtomId(++nAtomId);
          oAtom->SetAtomicNo(8);
          oAtom->SetAtomicMass(16.000);
          oAtom->SetHybridState(Atom::SP3);
          oAtom->SetVdwRadius(1.35);
          h1Atom->SetAtomId(++nAtomId);
          h1Atom->SetAtomicNo(1);
          h1Atom->SetAtomicMass(1.008);
          h1Atom->SetHybridState(Atom::SP3);
          h1Atom->SetVdwRadius(0.5);
          h2Atom->SetAtomId(++nAtomId);
          h2Atom->SetAtomicNo(1);
          h2Atom->SetAtomicMass(1.008);
          h2Atom->SetHybridState(Atom::SP3);
          h2Atom->SetVdwRadius(0.5);
          // Create the atom list and bond list for this water
          AtomList waterAtomList;
          BondList waterBondList;
          waterAtomList.push_back(oAtom);
          waterAtomList.push_back(h1Atom);
          waterAtomList.push_back(h2Atom);
          waterBondList.push_back(new Bond(++nBondId, oAtom, h1Atom, 1));
          waterBondList.push_back(new Bond(++nBondId, oAtom, h2Atom, 1));
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
          ModelPtr solvent(new Model(waterAtomList, waterBondList));
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

void PRMFactory::AttachReceptorFlexData(Model *pReceptor) {
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
    FlexData *pFlexData = new ReceptorFlexData(m_pDS);
    pFlexData->SetParameter(ReceptorFlexData::_FLEX_DISTANCE, flexDist);
    if (m_pParamSource->isParameterPresent(_REC_DIHEDRAL_STEP)) {
      double dihedralStepSize =
          m_pParamSource->GetParameterValue(_REC_DIHEDRAL_STEP);
      pFlexData->SetParameter(ReceptorFlexData::_DIHEDRAL_STEP,
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

void PRMFactory::AttachLigandFlexData(Model *pLigand) {
  m_pParamSource->SetSection(_LIG_SECTION);
  std::vector<std::string> paramList = m_pParamSource->GetParameterList();
  FlexData *pFlexData = new LigandFlexData(m_pDS);
  for (std::vector<std::string>::const_iterator iter = paramList.begin();
       iter != paramList.end(); ++iter) {
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

void PRMFactory::AttachSolventFlexData(Model *pSolvent) {
  m_pParamSource->SetSection(_SOLV_SECTION);
  std::vector<std::string> paramList = m_pParamSource->GetParameterList();
  FlexData *pFlexData = new SolventFlexData(m_pDS);
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
  AtomList atomList = pSolvent->GetAtomList();
  double occupancy(1.0);
  ChromElement::eMode transMode(ChromElement::FIXED);
  ChromElement::eMode rotMode(ChromElement::FIXED);
  if (!atomList.empty()) {
    Atom *pFirstAtom = atomList[0];
    occupancy = pFirstAtom->GetUser1Value();
    occupancy = std::min(1.0, occupancy);
    occupancy = std::max(0.0, occupancy);
    int mode = (int)pFirstAtom->GetUser2Value();
    mode = std::min(8, mode);
    mode = std::max(0, mode);
    transMode = (ChromElement::eMode)(mode / 3);
    rotMode = (ChromElement::eMode)(mode % 3);
  }
  pFlexData->SetParameter(LigandFlexData::_TRANS_MODE,
                          ChromElement::ModeToStr(transMode));
  pFlexData->SetParameter(LigandFlexData::_ROT_MODE,
                          ChromElement::ModeToStr(rotMode));
  pFlexData->SetParameter(SolventFlexData::_OCCUPANCY, occupancy);
  // Set the remaining flexibility parameters from the <receptor>.prm file
  // source
  for (std::vector<std::string>::const_iterator iter = paramList.begin();
       iter != paramList.end(); ++iter) {
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

MolecularFileSourcePtr
PRMFactory::CreateMolFileSource(const std::string &fileName) {
  MolecularFileSourcePtr retVal;
  std::string fileType = GetFileType(fileName);
  std::string fileTypeUpper;
  // std::toupper will not work https://stackoverflow.com/a/7131881
  std::transform(fileType.begin(), fileType.end(),
                 std::back_inserter(fileTypeUpper), ::toupper);
  std::string fullFileName = GetDataFileName("", fileName);

  if (m_iTrace > 0) {
    std::cout << _CT << ": File name requested = " << fileName
              << "; type = " << fileTypeUpper << std::endl;
    std::cout << _CT << ": Reading file from " << fullFileName << std::endl;
  }

  bool bImplH = true;
  if (m_pParamSource->isParameterPresent("RECEPTOR_ALL_H")) {
    Variant vAllH(m_pParamSource->GetParameterValueAsString("RECEPTOR_ALL_H"));
    bImplH = !vAllH.GetBool();
  }

  if (fileTypeUpper == "MOL2") {
    retVal = new MOL2FileSource(fullFileName, bImplH);
  } else if (fileTypeUpper == "PSF") {
    std::string strMassesFile =
        m_pParamSource->GetParameterValueAsString("RECEPTOR_MASSES_FILE");
    strMassesFile = GetDataFileName("data", strMassesFile);
    if (m_iTrace > 0) {
      std::cout << _CT << ": Using file " << strMassesFile
                << " to lookup Charmm atom type info" << std::endl;
    }
    retVal = new PsfFileSource(fullFileName, strMassesFile, bImplH);
  } else if (fileTypeUpper == "CRD") {
    retVal = new CrdFileSource(fullFileName);
    bImplH = false;
  } else if (fileTypeUpper == "PDB") {
    retVal = new PdbFileSource(fullFileName);
    bImplH = false;
  } else if ((fileTypeUpper == "SD") || (fileTypeUpper == "SDF")) {
    retVal = new MdlFileSource(fullFileName, false, false, bImplH);
  } else {
    throw BadReceptorFile(_WHERE_,
                          fileType + " is not a supported molecular file type");
  }

  if (m_pParamSource->isParameterPresent("RECEPTOR_SEGMENT_NAME")) {
    std::string strSegmentName =
        m_pParamSource->GetParameterValueAsString("RECEPTOR_SEGMENT_NAME");
    SegmentMap segmentMap = ConvertStringToSegmentMap(strSegmentName);
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

  Error status = retVal->Status();
  if (!status.isOK())
    throw status;
  int nAtoms = retVal->GetNumAtoms();
  if (nAtoms == 0) {
    if (m_iTrace > 0) {
      std::cout << _CT << ": File source contains zero atoms!" << std::endl;
      std::cout << _CT << ": Have you defined RECEPTOR_SEGMENT_NAME correctly?"
                << std::endl;
    }
    throw BadReceptorFile(_WHERE_,
                          "Zero atoms provided by " + retVal->GetFileName());
  } else {
    if (m_iTrace > 0) {
      std::cout << _CT << ": File source contains " << nAtoms << " atoms"
                << std::endl;
    }
  }
  return retVal;
}
