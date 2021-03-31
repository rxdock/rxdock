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

//#include <cstdlib>				// for getenv
//#include "rxdock/FileError.h"
#include "rxdock/DirectorySource.h"

#include <loguru.hpp>

using namespace rxdock;

const std::string DirectorySource::_CT = "DirectorySource";

/**
 * Constructor checks wether the directory exists and gives
 * $RBT_ROOT to the path if required
 */
DirectorySource::DirectorySource(const std::string &aDirectory) {
  // std::string theFileStr;
  //// filename without path after scandir() std::string theLine;
  //// line read

  thePath = (std::string)std::getenv("RBT_ROOT"); // path string to the dir
  if (thePath.empty())                            // $RBT_ROOT not defined ...
    throw EnvNotDefined(_WHERE_, "$RBT_ROOT ");
  else if (aDirectory.size() > PATH_SIZE) // or string is too long
    throw StringTooLong(_WHERE_, " ");
  else if (stat(aDirectory.c_str(), &fStat) <
           0)             // check in the cwd (or absolute path)
    thePath = aDirectory; // well, we want to read from cwd/abs path
  else {
    // thePath =+ "/"+aDirectory;					//
    // concatenate
    thePath = aDirectory; // providing $RBT_ROOT is OK, we do not have to give /
    if (thePath.size() > PATH_SIZE) // again too long
      throw StringTooLong(_WHERE_, "");
    else if (stat(thePath.c_str(), &fStat) <
             0) // if not even in the extended path
      throw DirIsNotAccessible(_WHERE_, "");
  }
  // is it _really_ a directory or just an entry?
  CheckDirectory(thePath); // throws exceptions if something messy
}

void DirectorySource::CheckDirectory(std::string &aDir) {
  //#ifdef COMENT
  if (!S_ISDIR(fStat.st_mode))
    throw DirIsNotAccessible(_WHERE_, aDir + " : Not a directory. ");
  char *theDir =
      (char *)aDir.c_str(); // unbelivable annoying buglet in Solaris gcc 2.95.3
  fNum = scandir(theDir, &fNameList, 0, alphasort); // get dir entries
  //    LOG_F(1, "Not running on Slowaris");
  //    std::exit(1);
  if (fNum < 1)
    throw NoFileInDir(_WHERE_, " No files in directory. ");
  else
    LOG_F(INFO, "DirectorySource::CheckDirectory: Reading files from {}", aDir);
  //#endif
}

DirectorySource::~DirectorySource() { LOG_F(2, "DirectorySource destructor"); }
