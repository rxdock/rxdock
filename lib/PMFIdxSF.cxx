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

#include "rxdock/PMFIdxSF.h"
#include "rxdock/Atom.h"
#include "rxdock/PMFDirSource.h"
#include "rxdock/WorkSpace.h"

#include <loguru.hpp>

#include <functional>

using namespace rxdock;

// the PMFs are defined only in this range in the Muegge database
const double cPMFStart = 0.2; // closest PMF distance
const double cPMFRes = 0.2;   // PMF resolution
const double cPMFEnd = 12.0;  // farest point in PMFs
const unsigned int cPlStart =
    1; // index value to get where the PMF plateau starts
const unsigned int cPlVal = 2; // index value to get the PMF plateau value

const std::string PMFIdxSF::_CT = "PMFIdxSF";
const std::string PMFIdxSF::_PMFDIR = "PMFDIR";
const std::string PMFIdxSF::_CC_CUTOFF = "CC_CUTOFF";
const std::string PMFIdxSF::_SLOPE = "SLOPE";

double delta; // used for linear interpolation

PMFIdxSF::PMFIdxSF(const std::string &aName) : BaseSF(_CT, aName) {
  // see PMF-related .prm files for explanation
  AddParameter(_PMFDIR, "data/pmf");
  AddParameter(_CC_CUTOFF, 6.0);
  AddParameter(_SLOPE, -3.0);
  delta = cPMFRes /
          2.0; // half of the PMF grid resoluton: delta for linear interpolation
  // create the PMF pseudogrid
  int nTypes = 37; // must be changed when we are defining new types :(
  Coord thePMFGridMin(cPMFStart, 0.0, 0.0);
  Coord thePMFGridStep(cPMFRes, 1.0, 1.0);
  thePMFGrid = RealGridPtr(new RealGrid(
      thePMFGridMin, thePMFGridStep, (int)(cPMFEnd / cPMFRes), nTypes, nTypes));
  // we are building a grid for slopes, so we can get the plateau
  // info similarly like the PMF values
  Coord theSlopeGridCoords(1.0, 1.0, 1.0);
  theSlopeGrid = RealGridPtr(new RealGrid(
      theSlopeGridCoords, theSlopeGridCoords, 2, // cPlStart & cPlVal
      nTypes, nTypes));

  // setup the grids for every types
  // first read .pmf files from directory
  PMFDirSource theSrcDir(GetRoot() + "/" + (std::string)GetParameter(_PMFDIR));
  std::vector<std::vector<PMFValue>> thePMF; // vector helping to read PMF files
  std::vector<std::string> theFileNames; // vector storing filenames (needed to
                                         // find out types in thePMF)
  std::vector<PMFValue>
      theSlopeIndex; // contains the first value where plateau starts
  theSrcDir.ReadFiles(&thePMF, &theFileNames, &theSlopeIndex);
  // convert vector to grid
  for (unsigned int i = 0; i < thePMF.size(); i++) {
    LOG_F(1, "PMFIdxSF: {} {} {} {} {}", theFileNames[i],
          theFileNames[i].substr(0, 2),
          PMFStr2Type(theFileNames[i].substr(0, 2)),
          theFileNames[i].substr(2, 2),
          PMFStr2Type(theFileNames[i].substr(2, 2)));

    // iterate through each individual PMF
    for (unsigned int j = 0; j < thePMF[i].size(); j++) {
      PMFType rType =
          PMFStr2Type(theFileNames[i].substr(0, 2)); // receptor type
      PMFType lType = PMFStr2Type(theFileNames[i].substr(2, 2)); // ligand type
      thePMFGrid->SetValue(
          thePMFGrid->GetIX((double)thePMF[i][j].distance), // distance
          rType, lType,
          thePMF[i][j].density // the density value
      );
    }
  }
  // fill slope grid
  for (unsigned int i = 0; i < theSlopeIndex.size(); i++) {
    PMFType rType = PMFStr2Type(theFileNames[i].substr(0, 2)); // receptor type
    PMFType lType = PMFStr2Type(theFileNames[i].substr(2, 2)); // ligand type
    theSlopeGrid->SetValue(cPlStart,     // index 1 for the distance where the
                           rType, lType, // plateau starts
                           theSlopeIndex[i].distance);
    theSlopeGrid->SetValue(cPlVal,       // that is the actual value
                           rType, lType, // can be more than 3.0
                           theSlopeIndex[i].density);
  }
  // checking values for type cP only. This test assumes the original
  // Muegge file list where the very first is the cPcP.pmf
  for (int j = 0; j < thePMF[0].size(); j++) {
    LOG_F(
        1, "PMFIdxSF: {}",
        thePMFGrid->GetValue(thePMFGrid->GetIX(thePMF[0][j].distance), cP, cP) -
            thePMF[0][j].density);
  }
  LOG_F(1, "PMFIdxSF: {}", theFileNames[0]);

  //	 this is a printout of values NOT the grid
  LOG_F(1, "PMFIdxSF: number of PMF functions: {}", thePMF.size());
  for (long k = 0; k < thePMF.size(); k++) {
    for (int l = 0, m = 0; l < thePMF[k].size(); l++, m++) {
      LOG_F(1, "l = {}", l);
      LOG_F(1, "| {}  {} |", thePMF[k][l].distance, thePMF[k][l].density);
      if (5 == m) {
        LOG_F(1, "m = 5, resetting to 0");
        m = 0;
      }
    }
  }

  // other test for grid. PMFType enum starts with 1.
  for (int i = 1; i < PMF_UNDEFINED; i++) {
    for (int j = 1; j < PMF_UNDEFINED; j++) {
      LOG_F(1, "PMFIdxSF: PMF Pair: {} {}", PMFType2Str((PMFType)i),
            PMFType2Str((PMFType)j));
      for (double grid = cPMFStart; grid < cPMFEnd; grid += cPMFRes) {
        unsigned int idx = thePMFGrid->GetIX(grid);
        LOG_F(1, "{} {}", PMFType2Str((PMFType)i), PMFType2Str((PMFType)j));
        LOG_F(1, "{} {} {}", grid, idx,
              thePMFGrid->GetValue(idx, (PMFType)i, (PMFType)j));
      }
    }
  }

  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

PMFIdxSF::~PMFIdxSF() {
  LOG_F(2, "PMFIdxSF PMF destructor");
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

void PMFIdxSF::Update(Subject *theChangedSubject) {
  LOG_F(2, "PMFIdxSF PMF Update");
  BaseInterSF::Update(theChangedSubject);
}

/**
 * SetupReceptor:
 * determine PMF atom types for all receptor atoms
 * create NonBondedGrid of nearest neighbours
 */
void PMFIdxSF::SetupReceptor() {
  LOG_F(2, "PMFIdxSF PMF SetupReceptor");
  // clear to be on the safe side
  theReceptorList.clear();
  theReceptorRList.clear();
  theSurround = NonBondedGridPtr();
  if (GetReceptor().Null()) {
    LOG_F(WARNING, "PMFIdxSF::SetupReceptor: no receptor defined");
    return;
  } else { // load PMFs from table files
    theSurround = CreateNonBondedGrid();
    DockingSitePtr spDS = GetWorkSpace()->GetDockingSite();
    double range = GetRange() + GetMaxError();
    AtomList atomList = spDS->GetAtomList(GetReceptor()->GetAtomList(), 0.0,
                                          GetCorrectedRange());
    theReceptorList =
        GetAtomListWithPredicate(atomList, std::not1(isAtomicNo_eq(1)));
    LOG_F(1, "Receptor list size: {}", theReceptorList.size());
    // create non-bonded grid to get the atoms around each gridpoint
    for (AtomListIter sIter = theReceptorList.begin();
         sIter != theReceptorList.end(); ++sIter) {
      theSurround->SetAtomLists(
          *sIter, range); // get surround at cutoff+gridstep*std::sqrt(3)/2
      // initialise cumulative PMF values for annotation
      (*sIter)->SetUser2Value(0.0);
    }
  }
  // transform smartpointers into regular ones
  std::copy(theReceptorList.begin(), theReceptorList.end(),
            std::back_inserter(theReceptorRList));
}

void PMFIdxSF::SetupLigand() {
  LOG_F(2, "PMFIdxSF PMF SetupLigand");
  theLigandList = GetAtomListWithPredicate(GetLigand()->GetAtomList(),
                                           std::not1(isAtomicNo_eq(1)));

  // clean old list
  theLigandRList.clear();
  // transform smartpointers into regular ones
  std::copy(theLigandList.begin(), theLigandList.end(),
            std::back_inserter(theLigandRList));
}

void PMFIdxSF::SetupScore() { LOG_F(2, "PMFIdxSF PMF SetupScore"); }

double PMFIdxSF::RawScore() const {
  LOG_F(2, "PMFIdxSF PMF RawScore");
  //	return 0.0;
  double theScore = 0.0;

  // check for existence of  atom list grid
  if (theSurround.Null()) {
    LOG_F(INFO, "PMFIdxSF::RawScore: No index grid");
    return theScore;
  }
  // enable/disable annotations
  bool bAnnotate = isAnnotationEnabled();

  // for all ligand atoms:
  for (AtomRListConstIter lIter = theLigandRList.begin();
       lIter != theLigandRList.end(); ++lIter) {
    const Coord &ligCoord = (*lIter)->GetCoords();
    // get receptor atoms that are within the PMF radius - if there are any
    const AtomRList &rAtomList = theSurround->GetAtomList(ligCoord);
    if (rAtomList.empty())
      continue;
    LOG_F(1, "PMFIdxSF::RawScore: loop over the receptor atoms around the "
             "ligand atom");
    for (AtomRListConstIter rIter = rAtomList.begin(); rIter != rAtomList.end();
         rIter++) {
      // get distance of the atom
      double theDist = Length((*rIter)->GetCoords(), ligCoord);
      double i_score; // interpolated score

      if (theDist > GetRange()) // skip distances out of a given distance
        continue;
      LOG_F(1, "PMFIdxSF::RawScore: theDist = {}", theDist);
      // to make code readable introduce some local variables
      // hoping cc can optimize them out
      const PMFType rType = (*rIter)->GetPMFType();
      LOG_F(1, "PMFIdxSF::RawScore: receptor type {}", PMFType2Str(rType));
      const PMFType lType = (*lIter)->GetPMFType();
      // optimal distance for C-C interactions is
      // under 6A. Note NC is the next item in PMFType
      // after the carbon types
      if (theDist > (double)GetParameter(_CC_CUTOFF) && rType < NC &&
          lType < NC)
        continue;
      // if we are in the plateau region
      else if (theDist < theSlopeGrid->GetValue(1, rType, lType)) {
        i_score = GetLinearCloseRangeValue(theDist, rType, lType);
        LOG_F(1, "PMFIdxSF::RawScore: Plateau reached: {} {} {} {}", theDist,
              theSlopeGrid->GetValue(1, rType, lType), PMFType2Str(rType),
              PMFType2Str(lType));
      } else {
        if (HH == rType || HL == lType || HL == rType || HH == lType) {
          LOG_F(1, "{} {} should not happen, ligand type {}",
                PMFType2Str(rType), PMFType2Str(lType), PMFType2Str(lType));
        }
        // unsigned int idx = thePMFGrid->GetIX(theDist);
        // double score = thePMFGrid->GetValue(idx, rType, lType);
        // theScore += score;
        // make a linear interpolation
        unsigned int inf_idx = thePMFGrid->GetIX(theDist - delta);
        double inf_score = thePMFGrid->GetValue(inf_idx, rType, lType);
        unsigned int sup_idx = thePMFGrid->GetIX(theDist + delta);
        double sup_score = thePMFGrid->GetValue(sup_idx, rType, lType);
        // now calculate the distances from the gridpoints
        double inf_d = (theDist - thePMFGrid->GetXCoord(inf_idx)) / cPMFRes;
        double sup_d = 1.0 - inf_d;
        // weight the score with the distances from gridpoints
        i_score = inf_score * sup_d + sup_score * inf_d;
      }
      // store (increment) contribution of receptor atom
      (*rIter)->SetUser2Value((*rIter)->GetUser2Value() + i_score);
      theScore += i_score;
    }
  }
  // check PMF scores: the two values should be the same
  double receptorPMFscore = 0.0;
  for (AtomRListConstIter sIter = theReceptorRList.begin();
       sIter != theReceptorRList.end(); ++sIter) {
    receptorPMFscore += (*sIter)->GetUser2Value();
  }
  LOG_F(1, "PMFIdxSF::RawScore: Receptor score is: {}", receptorPMFscore);
  LOG_F(1, "PMFIdxSF::RawScore: PMF score is     : {}", theScore);

  // save annotation (if needed)
  // since there is no ligand side here we are assigning the the very first
  // ligand atom also no distance
  if (bAnnotate) {
    for (AtomRListConstIter sIter = theReceptorRList.begin();
         sIter != theReceptorRList.end(); ++sIter) {
      Atom *lAtom = (*theLigandRList.begin());
      AnnotationPtr spAnnotation(new Annotation(
          lAtom,                     // ligand
          (*sIter),                  // receptor
          0.0,                       // distance
          (*sIter)->GetUser2Value()) // cumulative receptor PMF contribution
      );
      AddAnnotation(spAnnotation);
      (*sIter)->SetUser2Value(0.0); // clean old values
    }
  }

  return theScore;
}

double PMFIdxSF::GetLinearCloseRangeValue(double aDist, PMFType aRecType,
                                          PMFType aLigType) const {
  double thePlateauStart = theSlopeGrid->GetValue(cPlStart, aRecType, aLigType);
  double thePlateauVal = theSlopeGrid->GetValue(cPlVal, aRecType, aLigType);
  double theSlope = (double)GetParameter(_SLOPE);

  return theSlope * aDist - theSlope * thePlateauStart + thePlateauVal;
}
