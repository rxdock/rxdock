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

#include <iomanip>

#include "RbtMdlFileSink.h"

using namespace rxdock;

////////////////////////////////////////
// Constructors/destructors
RbtMdlFileSink::RbtMdlFileSink(const std::string &fileName, RbtModelPtr spModel)
    : RbtBaseMolecularFileSink(fileName, spModel), m_bFirstRender(true) {
  // DM 27 Apr 1999 - default behaviour is for the first Render to overwrite any
  // existing file then subsequent Renders to append.
  SetAppend(false);
  // Open an Element data source
  m_spElementData = RbtElementFileSourcePtr(
      new RbtElementFileSource(GetRbtFileName("data", "RbtElements.dat")));
  _RBTOBJECTCOUNTER_CONSTR_("RbtMdlFileSink");
}

RbtMdlFileSink::~RbtMdlFileSink() {
  _RBTOBJECTCOUNTER_DESTR_("RbtMdlFileSink");
}

////////////////////////////////////////
// Public methods
////////////////
//
////////////////////////////////////////
// Override public methods from RbtBaseFileSink
void RbtMdlFileSink::Render() {
  try {
    RbtModelPtr spModel(GetModel());
    RbtAtomList modelAtomList(spModel->GetAtomList());
    RbtBondList modelBondList(spModel->GetBondList());
    RbtModelList solventList = GetSolvent();
    // Concatenate all solvent atoms and bonds into a single list
    // DM 7 June 2006 - only render the enabled solvent models
    RbtAtomList solventAtomList;
    RbtBondList solventBondList;
    for (RbtModelListConstIter iter = solventList.begin();
         iter != solventList.end(); ++iter) {
      if ((*iter)->GetEnabled()) {
        RbtAtomList atomList = (*iter)->GetAtomList();
        RbtBondList bondList = (*iter)->GetBondList();
        std::copy(atomList.begin(), atomList.end(),
                  std::back_inserter(solventAtomList));
        std::copy(bondList.begin(), bondList.end(),
                  std::back_inserter(solventBondList));
      }
    }

    // Title lines
    std::vector<std::string> titleList(spModel->GetTitleList());
    if (!titleList.empty())
      AddLine(titleList[0]);
    else
      AddLine(spModel->GetName());

    // Molfile allows up to 8 characters for program name
    std::string programName = GetProgramName();
    if (programName.length() > 8)
      programName = programName.substr(0, 8);
    programName.insert(programName.length(), 8 - programName.length(), ' ');

    // File timestamp is 10 characters
    std::time_t t = std::time(nullptr);
    char ftime[11];
    std::strftime(ftime, sizeof(ftime), "%m%d%y%H%M", std::localtime(&t));

    // First two characters are user initials, leave them empty for now
    AddLine("  " + programName + std::string(ftime) + "3D");
    AddLine(GetProduct() + "/" + GetProgramVersion());

    // Write number of atoms and bonds
    std::ostringstream ostr;
    ostr << std::setw(3) << modelAtomList.size() + solventAtomList.size()
         << std::setw(3) << modelBondList.size() + solventBondList.size()
         << std::setw(3) << 0 << std::setw(3) << 0 << std::setw(3) << 0
         << std::setw(3) << 0 << std::setw(3) << 0 << std::setw(3) << 0
         << std::setw(3) << 0 << std::setw(3) << 0 << std::setw(3) << 999
         << " V2000";
    AddLine(ostr.str());

    // DM 19 June 2006 - clear the map of logical atom IDs each time
    // we render a model
    m_atomIdMap.clear();
    RenderAtomList(modelAtomList);
    RenderAtomList(solventAtomList);
    RenderBondList(modelBondList);
    RenderBondList(solventBondList);
    RenderData(spModel->GetDataMap());

    AddLine("$$$$");

    Write();
    if (m_bFirstRender) {
      SetAppend(true);
      m_bFirstRender = false;
    }
  } catch (RbtError &error) {
    // Close();
    throw; // Rethrow the RbtError
  }
}

void RbtMdlFileSink::RenderAtomList(const RbtAtomList &atomList) {
  for (RbtAtomListConstIter aIter = atomList.begin(); aIter != atomList.end();
       aIter++) {
    RbtAtomPtr spAtom(*aIter);
    // DM 19 June 2006. Check if this atom has been rendered previously.
    if (m_atomIdMap.find(spAtom) == m_atomIdMap.end()) {
      unsigned int nextAtomId = m_atomIdMap.size() + 1;
      // std::cout << "RenderAtom " << spAtom->GetFullAtomName() << " (actual
      // ID=" << spAtom->GetAtomId()
      //	 << "); file ID=" << nextAtomId << std::endl;
      m_atomIdMap.insert(std::make_pair(spAtom, nextAtomId));
    } else {
      // Should never happen. Probably best to throw an error at this point.
      throw RbtBadArgument(_WHERE_, "Duplicate atom rendered: " +
                                        spAtom->GetFullAtomName());
    }
    RbtElementData elData =
        m_spElementData->GetElementData(spAtom->GetAtomicNo());
    int nFormalCharge = spAtom->GetFormalCharge();
    if (nFormalCharge != 0)
      nFormalCharge = 4 - nFormalCharge;
    std::ostringstream ostr;
    ostr.precision(4);
    ostr.setf(std::ios_base::fixed, std::ios_base::floatfield);
    ostr.setf(std::ios_base::right, std::ios_base::adjustfield);
    ostr << std::setw(10) << spAtom->GetX() << std::setw(10) << spAtom->GetY()
         << std::setw(10) << spAtom->GetZ(); // X,Y,Z coord
    ostr.setf(std::ios_base::left, std::ios_base::adjustfield);
    ostr << std::setw(0) << " " << std::setw(3)
         << elData.element.c_str(); // Element name
    ostr.setf(std::ios_base::right, std::ios_base::adjustfield);
    ostr << std::setw(2) << 0             // mass difference
         << std::setw(3) << nFormalCharge // charge
         << std::setw(3) << 0             // atom stereo parity
         << std::setw(3) << 0             // hydrogen count+1 (query CTABs only)
         << std::setw(3) << 0             // stereo care box (query CTABs only)
         << std::setw(3) << 0;            // valence (0 = no marking)
    // Mass diff, formal charge, stereo parity, num hydrogens,
    // center
    AddLine(ostr.str());
  }
}

void RbtMdlFileSink::RenderBondList(const RbtBondList &bondList) {
  for (RbtBondListConstIter bIter = bondList.begin(); bIter != bondList.end();
       bIter++) {
    RbtBondPtr spBond(*bIter);
    // DM 19 June 2006. Lookup the logical atom IDs for each atom in the bond
    RbtAtomIdMap::const_iterator aIter1 =
        m_atomIdMap.find(spBond->GetAtom1Ptr());
    RbtAtomIdMap::const_iterator aIter2 =
        m_atomIdMap.find(spBond->GetAtom2Ptr());
    if ((aIter1 != m_atomIdMap.end()) && (aIter2 != m_atomIdMap.end())) {
      unsigned int id1 = (*aIter1).second;
      unsigned int id2 = (*aIter2).second;
      // std::cout << "RenderBond " << spBond->GetAtom1Ptr()->GetFullAtomName()
      //	 << spBond->GetAtom2Ptr()->GetFullAtomName()
      //	 << "; file ID1=" << id1
      //	 << "; file ID2=" << id2 << std::endl;
      std::ostringstream ostr;
      ostr.setf(std::ios_base::right, std::ios_base::adjustfield);
      ostr << std::setw(3) << id1 << std::setw(3) << id2 << std::setw(3)
           << spBond->GetFormalBondOrder() << std::setw(3) << 0 << std::setw(3)
           << 0 << std::setw(3) << 0; // Atom1, Atom2, bond order, stereo
                                      // designator, unused, topology code
      AddLine(ostr.str());
    } else {
      // Should never happen. Probably best to throw an error at this point.
      throw RbtBadArgument(_WHERE_,
                           "Error rendering bond, logical atom IDs not found");
    }
  }
}

void RbtMdlFileSink::RenderData(const RbtStringVariantMap &dataMap) {
  if (!dataMap.empty()) {
    AddLine("M  END"); // End of molecule marker
    for (RbtStringVariantMapConstIter iter = dataMap.begin();
         iter != dataMap.end(); iter++) {
      AddLine(">  <" + (*iter).first + ">"); // Field name
      std::vector<std::string> sl = (*iter).second.StringList();
      for (std::vector<std::string>::const_iterator slIter = sl.begin();
           slIter != sl.end(); ++slIter) {
        AddLine(*slIter); // Field values
      }
      AddLine(""); // Blank line denotes end of field values
    }
  }
}
