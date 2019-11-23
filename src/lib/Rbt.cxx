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

// Misc non-member functions in Rbt namespace

#include <algorithm> //For sort
#include <climits>   //For PATH_MAX
#include <cstdlib>   //For getenv
#include <ctime>     //For time functions
#include <dirent.h>  //For directory handling
#include <fstream>   //For ifstream
#ifdef _WIN32
#include <direct.h>
#define getcwd _getcwd // CRT library getcwd
#else
#include <unistd.h> //For POSIX getcwd
#endif
//#include <ios>
#include "Rbt.h"
#include "RbtFileError.h"
#include "RbtResources.h"

// GetRbtRoot - returns value of RBT_ROOT env variable
std::string Rbt::GetRbtRoot() {
  char *szRbtRoot = std::getenv("RBT_ROOT");
  if (szRbtRoot != (char *)nullptr) {
    return std::string(szRbtRoot);
  } else {
    return GetCurrentWorkingDirectory();
  }
}

// DM 02 Aug 2000
// GetRbtHome - returns value of RBT_HOME env variable
// or HOME if RBT_HOME is not defined
// If HOME is undefined, returns current working directory
std::string Rbt::GetRbtHome() {
  char *szRbtHome = std::getenv("RBT_HOME");
  if (szRbtHome != (char *)nullptr) {
    return std::string(szRbtHome);
  } else {
    szRbtHome = std::getenv("HOME");
    if (szRbtHome != (char *)nullptr) {
      return std::string(szRbtHome);
    } else {
      return GetCurrentWorkingDirectory();
    }
  }
}

// Rbt::GetCopyright - returns legalese statement
std::string Rbt::GetCopyright() { return IDS_COPYRIGHT; }
// Rbt::GetVersion - returns current library version
std::string Rbt::GetVersion() { return IDS_VERSION; }
// GetBuildNum - returns current library build number
std::string Rbt::GetBuild() { return IDS_BUILDNUM; }
// GetProduct - returns library product name
std::string Rbt::GetProduct() { return IDS_PRODUCT; }
// GetTime - returns current time as an RbtString
std::string Rbt::GetTime() {
  std::time_t t = std::time(nullptr);       // Get time in seconds since 1970
  std::tm *pLocalTime = std::localtime(&t); // Convert to local time struct
  return std::string(std::asctime(pLocalTime)); // Convert to ascii string
}

// GetCurrentWorkingDirectory - returns current working directory
std::string Rbt::GetCurrentWorkingDirectory() {
  std::string strCwd(".");
  char *szCwd = new char[PATH_MAX + 1]; // Allocate a temp char* array
  if (::getcwd(szCwd, PATH_MAX) != (char *)nullptr) { // Get the cwd
    strCwd = szCwd;
    // strCwd += "/";
  }
  delete[] szCwd; // Delete the temp array
  return strCwd;
}

// Rbt::GetRbtDirName
// Returns the full path to a subdirectory in the rDock directory structure
//
// For example, if RBT_ROOT environment variable is ~dave/ribodev/molmod/ribodev
// then GetRbtDirName("data") would return ~dave/ribodev/molmod/ribodev/data/
//
std::string Rbt::GetRbtDirName(const std::string &strSubDir) {
  std::string strRbtDir = GetRbtRoot();
  if (strSubDir.size() > 0) {
    strRbtDir += "/";
    strRbtDir += strSubDir;
  }
  return strRbtDir;
}

// Rbt::GetRbtFileName
// DM 17 Dec 1998 - slightly different behaviour
// First check if the file exists in the CWD, if so return this path
// Next check RBT_HOME directory, if so return this path
// Finally, return the path to the file in the rDock directory structure
//(without checking if the file is actually present)
std::string Rbt::GetRbtFileName(const std::string &strSubdir,
                                const std::string &strFile) {
  // First see if the file exists in the current directory
  std::string strFullPathToFile(strFile);
  // Just open it, don't try and parse it (after all, we don't know what format
  // the file is)
  std::ifstream fileIn(strFullPathToFile.c_str(), std::ios_base::in);
  if (fileIn) {
    fileIn.close();
    return strFullPathToFile;
  } else {
    // DM 02 Aug 200 - check RBT_HOME directory
    strFullPathToFile = GetRbtHome() + "/" + strFile;
    // DM 27 Apr 2005 - under gcc 3.4.3 there are problems reusing the same
    // ifstream object after the first "file open" fails
    // fileIn.open(strFullPathToFile.c_str(),std::ios_base::in);
    std::ifstream fileIn2(strFullPathToFile.c_str(), std::ios_base::in);
    if (fileIn2) {
      fileIn2.close();
      return strFullPathToFile;
    } else {
      return GetRbtDirName(strSubdir) + "/" + strFile;
    }
  }
}

// GetFileType
// Returns the string following the last "." in the file name.
// e.g. GetFileType("receptor.psf") would return "psf"
// If no "." is present, returns the whole file name
std::string Rbt::GetFileType(const std::string &strFile) {
  return strFile.substr(strFile.rfind(".") + 1);
}

// Rbt::GetDirList
// Returns a list of files in a directory (strDir) whose names begin with
// strFilePrefix (optional) and whose type is strFileType (optional, as returned
// by GetFileType)
std::vector<std::string> Rbt::GetDirList(const std::string &strDir,
                                         const std::string &strFilePrefix,
                                         const std::string &strFileType) {
  std::vector<std::string> dirList;
  bool bMatchPrefix =
      (strFilePrefix.size() > 0); // Check if we need to match on file prefix
  bool bMatchType = (strFileType.size() >
                     0); // Check if we need to match on file type (suffix)

  DIR *pDir = opendir(strDir.c_str()); // Open directory
  if (pDir != nullptr) {

    // DM 6 Dec 1999 - dirent_t appears to be SGI MIPSPro specific
    // At least on Linux g++, it is called dirent
#ifdef __sgi
    dirent_t *pEntry;
#else
    dirent *pEntry;
#endif

    while ((pEntry = readdir(pDir)) != nullptr) { // Read each entry
      std::string strFile = pEntry->d_name;
      if ((strFile == ".") || (strFile == ".."))
        continue; // Don't need to consider . and ..
      bool bMatch = true;
      if (bMatchPrefix && (strFile.find(strFilePrefix) !=
                           0)) // Eliminate those that don't match the prefix
        bMatch = false;
      if (bMatchType &&
          (GetFileType(strFile) !=
           strFileType)) // Eliminate those that don't match the type
        bMatch = false;
      if (bMatch) // Only store the hits
        dirList.push_back(strFile);
    }
    closedir(pDir);
  }

  // Sort the file list into alphabetical order
  std::sort(dirList.begin(), dirList.end());

  return dirList;
}

// Converts (comma)-delimited string of segment names to segment map
RbtSegmentMap Rbt::ConvertStringToSegmentMap(const std::string &strSegments,
                                             const std::string &strDelimiter) {
#ifdef _DEBUG
  // std::cout << "ConvertStringToSegmentMap: " << strSegments << " delimiter="
  // << strDelimiter << std::endl;
#endif //_DEBUG

  std::string::size_type nDelimiterSize = strDelimiter.size();
  RbtSegmentMap segmentMap;

  // Check for null string or null delimiter
  if ((strSegments.size() > 0) && (nDelimiterSize > 0)) {
    std::string::size_type iBegin =
        0; // indicies into string (sort-of iterators)
    std::string::size_type iEnd;
    // Not often a do..while loop is used, but in this case it's what we want
    // Even if no delimiter is present, we still need to extract the whole
    // string
    do {
      iEnd = strSegments.find(strDelimiter, iBegin);
#ifdef _DEBUG
      // std::cout << strSegments.substr(iBegin, iEnd-iBegin) << std::endl;
#endif //_DEBUG
      segmentMap[strSegments.substr(iBegin, iEnd - iBegin)] = 0;
      iBegin = iEnd + nDelimiterSize;
    } while (iEnd !=
             std::string::npos); // This is the check for delimiter not found
  }

  return segmentMap;
}

// Converts segment map to (comma)-delimited string of segment names
std::string Rbt::ConvertSegmentMapToString(const RbtSegmentMap &segmentMap,
                                           const std::string &strDelimiter) {
  std::string strSegments;

  // Check for empty segment map
  if (segmentMap.size() > 0) {
    RbtSegmentMapConstIter iter = segmentMap.begin();
    strSegments += (*iter++).first; // Add first string
    // Now loop over remaining entries, adding delimiter before each string
    while (iter != segmentMap.end()) {
      strSegments += strDelimiter;
      strSegments += (*iter++).first;
    }
  }
  return strSegments;
}

// Returns a segment map containing the members of map1 which are not in map2
// I know, should really be a template so as to be more universal...one day
// maybe. Or maybe there is already an STL algorithm for doing this.
RbtSegmentMap Rbt::SegmentDiffMap(const RbtSegmentMap &map1,
                                  const RbtSegmentMap &map2) {
  RbtSegmentMap map3 = map1; // Init return value to map1
  for (RbtSegmentMapConstIter iter = map2.begin(); iter != map2.end(); iter++)
    map3.erase((*iter).first); // Now delete everything in map2
  return map3;
}

// DM 30 Mar 1999
// Converts (comma)-delimited string to string list (similar to
// ConvertStringToSegmentMap, but returns list not map)
std::vector<std::string>
Rbt::ConvertDelimitedStringToList(const std::string &strValues,
                                  const std::string &strDelimiter) {
  std::string::size_type nDelimiterSize = strDelimiter.size();
  std::vector<std::string> listOfValues;

  // Check for null string or null delimiter
  if ((strValues.size() > 0) && (nDelimiterSize > 0)) {
    std::string::size_type iBegin =
        0; // indicies into string (sort-of iterators)
    std::string::size_type iEnd;
    // Not often a do..while loop is used, but in this case it's what we want
    // Even if no delimiter is present, we still need to extract the whole
    // string
    do {
      iEnd = strValues.find(strDelimiter, iBegin);
      std::string strValue = strValues.substr(iBegin, iEnd - iBegin);
      listOfValues.push_back(strValue);
      iBegin = iEnd + nDelimiterSize;
    } while (iEnd !=
             std::string::npos); // This is the check for delimiter not found
  }
  return listOfValues;
}

// Converts string list to (comma)-delimited string (inverse of above)
std::string
Rbt::ConvertListToDelimitedString(const std::vector<std::string> &listOfValues,
                                  const std::string &strDelimiter) {
  std::string strValues;

  // Check for empty string list
  if (!listOfValues.empty()) {
    std::vector<std::string>::const_iterator iter = listOfValues.begin();
    strValues += *iter++; // Add first string
    // Now loop over remaining entries, adding delimiter before each string
    while (iter != listOfValues.end()) {
      strValues += strDelimiter;
      strValues += *iter++;
    }
  }
  return strValues;
}

////////////////////////////////////////////////////////////////
// I/O ROUTINES
//
// PrintStdHeader
// Print a standard header to an output stream (for log files etc)
// Contains copyright info, library version, date, time etc
// DM 19 Feb 1999 - include executable information
std::ostream &Rbt::PrintStdHeader(std::ostream &s,
                                  const std::string &strExecutable) {
  s << "***********************************************" << std::endl;
  s << GetCopyright() << std::endl;
  if (!strExecutable.empty())
    s << "Executable:\t" << strExecutable << std::endl;
  s << "Library:\t" << GetProduct() << "/" << GetVersion() << "/" << GetBuild()
    << std::endl;
  s << "RBT_ROOT:\t" << GetRbtRoot() << std::endl;
  s << "RBT_HOME:\t" << GetRbtHome() << std::endl;
  s << "Current dir:\t" << GetCurrentWorkingDirectory() << std::endl;
  s << "Date:\t\t" << GetTime();
  s << "***********************************************" << std::endl;
  return s;
}

// Helper functions to read/write chars from iostreams
// Throws error if stream state is not Good() before and after the read/write
// It appears the STL std::ios_base exception throwing is not yet implemented
// at least on RedHat 6.1, so this is a temporary workaround (yeah right)
void Rbt::WriteWithThrow(std::ostream &ostr, const char *p, streamsize n) {
  if (!ostr)
    throw RbtFileWriteError(_WHERE_, "Error writing to output stream");
  ostr.write(p, n);
  if (!ostr)
    throw RbtFileWriteError(_WHERE_, "Error writing to output stream");
}

void Rbt::ReadWithThrow(std::istream &istr, char *p, streamsize n) {
  if (!istr)
    throw RbtFileReadError(_WHERE_, "Error reading from input stream");
  istr.read(p, n);
  if (!istr)
    throw RbtFileReadError(_WHERE_, "Error reading from input stream");
}

// Used to read RbtCoord. The separator between x y z should be a
// ',', but this takes care of small mistakes, reading any white
// space or commas there is between each variable.
// If necessary, it can be modified to accept the ',' as a
// parameter to be able to use it when other separators are
// needed. Also, it should be possible to implemente it as a
// manipulator of std::istream.
std::istream &eatSeps(std::istream &is) {
  char c;
  while (is.get(c)) {
    if (!isspace(c) && (c != ',')) {
      is.putback(c);
      break;
    }
  }
  return is;
}
