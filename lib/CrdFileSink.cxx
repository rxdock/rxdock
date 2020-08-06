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

#include "rxdock/CrdFileSink.h"

using namespace rxdock;

////////////////////////////////////////
// Constructors/destructors
// CrdFileSink::CrdFileSink(const char* fileName, ModelPtr spModel) :
//  BaseMolecularFileSink(fileName,spModel) //Call base class constructor
//{
//  SetMultiConf(true);//Default to allowing multiple conformers/ models in crd
//  files _RBTOBJECTCOUNTER_CONSTR_("CrdFileSink");
//}

CrdFileSink::CrdFileSink(const std::string &fileName, ModelPtr spModel)
    : BaseMolecularFileSink(fileName, spModel) // Call base class constructor
{
  SetMultiConf(
      true); // Default to allowing multiple conformers/ models in crd files
  _RBTOBJECTCOUNTER_CONSTR_("CrdFileSink");
}

CrdFileSink::~CrdFileSink() { _RBTOBJECTCOUNTER_DESTR_("CrdFileSink"); }

////////////////////////////////////////
// Public methods
////////////////
//
////////////////////////////////////////
// Override public methods from BaseFileSink
void CrdFileSink::Render() {
  // Reset() has the intelligence to only reset the counters when required
  // so we can simply call it each time here
  Reset();

  try {
    ModelPtr spModel(GetModel());

    // Only write titles and number of atoms if this is the first conformer
    if (m_nAtomId == 0) {
      // 1. Write the titles
      std::vector<std::string> titleList(spModel->GetTitleList());
      // DM 8 June 1999 - limit number of title lines output to a maximum of 1
      std::vector<std::string>::const_iterator iter;
      int nLines;
      for (iter = titleList.begin(), nLines = 0;
           (iter != titleList.end()) && (nLines < 1); iter++, nLines++) {
        // Add an initial star if the title does not already commence with a
        // star
        if ((*iter).find("*") != 0)
          AddLine("*" + (*iter));
        else
          AddLine(*iter);
      }
      // Title delimiter record (* on its own)
      AddLine("*");

      // 2. Write number of atoms
      std::ostringstream ostr;
      ostr << std::setw(5) << spModel->GetNumAtoms();
      AddLine(ostr.str());
      // Remember the line number containing the number of atoms
      // as we'll need to update the total num atoms after each Render when we
      // are in multiconf mode m_numAtomsLineRec = titleList.size()+1;
      m_numAtomsLineRec = nLines + 1;
    }

    // In multiconf mode, we need to keep track of the current subunit and
    // segment IDs so we can detect when they change
    std::string strLastSubunitId("UNLIKELY_TO_MATCH_THIS_BY_CHANCE");
    std::string strLastSegmentName("UNLIKELY_TO_MATCH_THIS_BY_CHANCE");

    // 3. Write the atom list
    AtomList atomList(spModel->GetAtomList());
    for (AtomListConstIter aIter = atomList.begin(); aIter != atomList.end();
         aIter++) {
      AtomPtr spAtom(*aIter);

      m_nAtomId++;
      // New subunit detected
      if (spAtom->GetSubunitId() != strLastSubunitId) {
        m_nSubunitId++;
        strLastSubunitId = spAtom->GetSubunitId();
      }
      // New segment detected
      if (spAtom->GetSegmentName() != strLastSegmentName) {
        m_nSegmentId++;
        strLastSegmentName = spAtom->GetSegmentName();
      }

      // Render the atom to a string stream
      std::ostringstream ostr;
      ostr.precision(5);
      ostr.setf(std::ios_base::fixed, std::ios_base::floatfield);
      ostr.setf(std::ios_base::right, std::ios_base::adjustfield);
      // In multiconf mode, use the incremental atom ID from the sink
      if (GetMultiConf())
        ostr << std::setw(5) << m_nAtomId;
      // Else use the atom ID from the atom itself
      else
        ostr << std::setw(5) << spAtom->GetAtomId();
      ostr << std::setw(5)
           << m_nSubunitId // This subunit ID is unique in the file
           << std::setw(0) << " ";
      ostr.setf(std::ios_base::left, std::ios_base::adjustfield);
      ostr << std::setw(4) << spAtom->GetSubunitName().c_str() << std::setw(0)
           << " " << std::setw(4) << spAtom->GetName().c_str();
      ostr.setf(std::ios_base::right, std::ios_base::adjustfield);
      ostr << std::setw(10) << spAtom->GetX() << std::setw(10) << spAtom->GetY()
           << std::setw(10) << spAtom->GetZ() << std::setw(0) << " ";

      // DM 21 Apr 1999 - if we are not using the model segment names, use the
      // incremental segment ID from the sink DM 11 Feb 1999 - output leading
      // zeros for segID, so they sort correctly when read as text
      if (!m_bUseModelSegmentNames) {
        ostr << std::setw(3) << std::setfill('0') << m_nSegmentId
             << std::setfill(' ') << " ";
        ostr.setf(std::ios_base::left, std::ios_base::adjustfield);
      }
      // Else use the segment name from the atom itself
      else {
        ostr.setf(std::ios_base::left, std::ios_base::adjustfield);
        ostr << std::setw(4) << spAtom->GetSegmentName().c_str();
      }

      ostr << std::setw(0) << " " << std::setw(4)
           << spAtom->GetSubunitId()
                  .c_str(); // This subunit ID is unique in the segment
      ostr.setf(std::ios_base::right, std::ios_base::adjustfield);
      ostr << std::setw(10) << 0.0;
      AddLine(ostr.str());
    }

    // In multiconf mode, we need to update the total number of atoms in the
    // cache then write the file without clearing the cache
    if (GetMultiConf()) {
      std::ostringstream ostr;
      ostr << std::setw(5) << m_nAtomId;
      ReplaceLine(ostr.str(), m_numAtomsLineRec);
      Write(false);
    } else
      Write(); // Commit the cache to the file
  } catch (Error &error) {
    // Close();
    throw; // Rethrow the Error
  }
}
