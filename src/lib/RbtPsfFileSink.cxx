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

#include "RbtPsfFileSink.h"
#include <cstdio>
#include <iomanip>
#include <sstream>

////////////////////////////////////////
// Constructors/destructors

RbtPsfFileSink::RbtPsfFileSink(const std::string &fileName, RbtModelPtr spModel)
    : RbtBaseMolecularFileSink(fileName, spModel) // Call base class constructor
{
  SetMultiConf(
      false); // Default to allowing multiple conformers/ models in crd files
  _RBTOBJECTCOUNTER_CONSTR_("RbtPsfFileSink");
}

RbtPsfFileSink::~RbtPsfFileSink() {
  _RBTOBJECTCOUNTER_DESTR_("RbtPsfFileSink");
}

////////////////////////////////////////
// Public methods
////////////////
//
////////////////////////////////////////
// Override public methods from RbtBaseFileSink
void RbtPsfFileSink::Render() {
  // String constants in PSF files
  const std::string strPsfKey("PSF");
  const std::string strTitleKey("!NTITLE");
  const std::string strAtomKey("!NATOM");
  const std::string strBondKey("!NBOND:");

  // Reset() has the intelligence to only reset the counters when required
  // so we can simply call it each time here
  Reset();

  try {
    RbtModelPtr spModel(GetModel());

    // 1. Write the titles
    AddLine(strPsfKey);
    AddLine("");
    RbtStringList titleList(spModel->GetTitleList());
    std::ostringstream ostr;
    ostr << std::setw(8) << titleList.size() << std::setw(0) << " "
         << strTitleKey << std::ends;
    AddLine(ostr.str());
    // delete ostr.str();
    for (RbtStringListConstIter iter = titleList.begin();
         iter != titleList.end(); iter++) {
      // Add an initial star if the title does not already commence with a star
      if ((*iter).find("*") != 0)
        AddLine("*" + (*iter));
      else
        AddLine(*iter);
    }
    AddLine("");

    // 2. Write number of atoms
    std::ostringstream ostr2;
    ostr2 << std::setw(8) << spModel->GetNumAtoms() << std::setw(0) << " "
          << strAtomKey << std::ends;
    AddLine(ostr2.str());
    // delete ostr2.str();

    // 3. Write the atom list
    RbtAtomList atomList(spModel->GetAtomList());
    char *line = new char[128];
    for (RbtAtomListConstIter aIter = atomList.begin(); aIter != atomList.end();
         aIter++) {
      RbtAtomPtr spAtom(*aIter);

      m_nAtomId++;
      std::sprintf(
          line, "%8d %-4.4s %-4.4s %-4.4s %-4.4s %-4.4s%15.3f%10.3f%12.1f",
          spAtom->GetAtomId(), spAtom->GetSegmentName().c_str(),
          spAtom->GetSubunitId().c_str(), spAtom->GetSubunitName().c_str(),
          spAtom->GetName().c_str(), spAtom->GetFFType().c_str(),
          spAtom->GetPartialCharge(), spAtom->GetAtomicMass(), 0.0);
      AddLine(line);
      // Render the atom to a string stream
      // std::ostringstream ostr;
      // ostr.precision(5);
      // ostr.setf(std::ios_base::fixed,ios_base::floatfield);
      // ostr.setf(std::ios_base::right,ios_base::adjustfield);
      // ostr.fill('#');
      // ostr << std::setw(8) << spAtom->GetAtomId() << " ";
      // ostr.setf(std::ios_base::left,ios_base::adjustfield);
      // ostr << std::setw(4) << spAtom->GetSegmentName() << " ";
      // ostr << std::setw(4) << spAtom->GetSubunitId() << " ";
      // ostr << std::setw(4) << spAtom->GetSubunitName() << " ";
      // ostr << std::setw(4) << spAtom->GetName() << " ";
      // ostr << std::setw(4) << spAtom->GetFFType() << " ";
      // ostr.setf(std::ios_base::right,ios_base::adjustfield);
      // ostr << std::setw(14) << spAtom->GetPartialCharge();
      // ostr << std::setw(10) << spAtom->GetAtomicMass();
      // ostr << std::setw(12) << 0.0 << std::ends;
      // AddLine(ostr.str());
      // delete ostr.str();
    }
    delete[] line;
    // 4. Write number of bonds
    std::ostringstream ostr3;
    ostr3 << std::setw(8) << spModel->GetNumBonds() << std::setw(0) << " "
          << strBondKey << std::ends;
    AddLine("");
    AddLine(ostr3.str());
    // delete ostr3.str();

    // 5. Write the bond list
    RbtBondList bondList(spModel->GetBondList());
    int nDiv = bondList.size() / 4;
    int nRem = bondList.size() % 4;
    RbtBondListConstIter bIter = bondList.begin();
    for (int i = 0; i < nDiv; i++) {
      std::ostringstream ostr;
      ostr << std::setw(8) << (*bIter)->GetAtom1Ptr()->GetAtomId()
           << std::setw(8) << (*bIter)->GetAtom2Ptr()->GetAtomId();
      bIter++;
      ostr << std::setw(8) << (*bIter)->GetAtom1Ptr()->GetAtomId()
           << std::setw(8) << (*bIter)->GetAtom2Ptr()->GetAtomId();
      bIter++;
      ostr << std::setw(8) << (*bIter)->GetAtom1Ptr()->GetAtomId()
           << std::setw(8) << (*bIter)->GetAtom2Ptr()->GetAtomId();
      bIter++;
      ostr << std::setw(8) << (*bIter)->GetAtom1Ptr()->GetAtomId()
           << std::setw(8) << (*bIter)->GetAtom2Ptr()->GetAtomId();
      bIter++;
      ostr << std::ends;
      AddLine(ostr.str());
      // delete ostr.str();
    }
    // Remaining bonds on final incomplete line
    if (nRem > 0) {
      std::ostringstream ostr;
      for (int i = 0; i < nRem; i++, bIter++) {
        ostr << std::setw(8) << (*bIter)->GetAtom1Ptr()->GetAtomId()
             << std::setw(8) << (*bIter)->GetAtom2Ptr()->GetAtomId();
      }
      ostr << std::ends;
      AddLine(ostr.str());
      // delete ostr.str();
    }

    Write(); // Commit the cache to the file
  } catch (RbtError &error) {
    throw; // Rethrow the RbtError
  }
}
