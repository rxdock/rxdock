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

#include "PMFDirSource.h"
// #include <stdlib>
//#include "FileError.h"

using namespace rxdock;

PMFDirSource::PMFDirSource(const std::string &aDir) : DirectorySource(aDir) {
#ifdef _DEBUG
  std::cout << _CT << " PMFDirSource constructor" << std::endl;
#endif
}

void PMFDirSource::ReadFiles(std::vector<std::vector<PMFValue>> *aVect,
                             std::vector<std::string> *aNameVect,
                             std::vector<PMFValue> *aSlopeVect) {
  std::cout << "Reading " << fNum << " PMF files..." << std::endl;
  std::string theFileName;             // will be the filename with full path
  std::string theLine;                 // one line from the file
  std::vector<std::string> theStrData; // data represented as strings in vector
  PMFValue theDummy = {0.0, 0.0};      // dummy initial value for plateau

  while (fNum--) {
    theFileName += thePath + "/";
    //         std::cout << "Does not work with Solaris" <<endl;
    //         std::exit(1);
    std::string theFileStr(fNameList[fNum]->d_name); // copy from C-string
    // std::string theFileStr("junk");
    if (theFileName.size() + theFileStr.size() >
        PATH_SIZE) // check size (though should be ok)
    {
      std::cout << _CT << theFileName << std::endl;
      std::cout << _CT << theFileStr << std::endl;
      throw StringTooLong(_WHERE_, "");
    } else
      theFileName += theFileStr; // concatenate to get final size
    // get dot index in extension
    std::string::size_type theExtIdx = theFileName.find_last_of(".pmf");
    // this is for storing values from one file
    std::vector<PMFValue> theValues;
    // when file has no .pmf extension, skip. Also check for existence and perm
    if (std::string::npos != theExtIdx && !stat(theFileName.c_str(), &fStat) &&
        S_ISREG(fStat.st_mode)) {
#ifdef _DEBUG
      std::cout << _CT << " Processing: " << theFileStr << std::endl;
#endif // debug
       // get rid of .pmf and put it into the vector that will be used later to
       // figure out types
      aNameVect->push_back(theFileStr.erase(theFileStr.find_last_of("."), 4));
      inFile.open(theFileName.c_str(), std::ifstream::in);
      // read file contents into vector
      while (inFile >> theLine)
        theStrData.push_back(theLine);
      inFile.close();
      ParseLines(theStrData, &theValues);
      aVect->push_back(theValues);
      // search for the start of the 3.0 plateau backwards in pmf
      aSlopeVect->push_back(theDummy); // create a dummy value
      for (int i = theValues.size(); i >= 0; i--) {
        if (theValues[i].density >= 3.0) {
          aSlopeVect->back() = theValues[i];
          // aSlopeVect->back().distance = theValues[i].distance;
          // aSlopeVect->back().density  = theValues[i].density;
          break;
        }
      }
#ifdef GETJUNK
      //#ifdef _DEBUG
      std::cout << _CT << " Starting slope at: " << aSlopeVect->back()
                << " for ";
      std::cout << theValues[i].density << " ";
      std::cout << theFileStr.substr(0, 2) << " " << theFileStr.substr(2, 2)
                << std::endl;
//#endif // _DEBUG
#endif // GETJUNK
      theStrData.erase(theStrData.begin(),
                       theStrData.end()); // delete vector after parsing
    }
#ifdef _DEBUG
    else {
      std::cout << _CT << " Skipping file: " << theFileName << std::endl;
    }
#endif //_DEBUG

    theFileName.erase();
  }
  std::cout << " done." << std::endl;
}

void PMFDirSource::ParseLines(std::vector<std::string> anStrVect,
                              std::vector<PMFValue> *aValueVect) {
  PMFValue theActualValue;

  for (unsigned int i = 0; i < anStrVect.size(); i++) {
    theActualValue.distance = (float)std::atof(anStrVect[i].data());
    i++; // get next number
    theActualValue.density = (float)std::atof(anStrVect[i].data());
    aValueVect->push_back(theActualValue);
  }
}
