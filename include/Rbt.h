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

// Miscellaneous non-member functions in Rbt namespace
#ifndef _RBT_H_
#define _RBT_H_

#include "RbtContainers.h"
#include "RbtError.h"

namespace Rbt {
////////////////////////////////////////////////////////////////
// RESOURCE HANDLING FUNCTIONS
//
// GetRbtRoot - returns value of RBT_ROOT env variable
std::string GetRbtRoot();
// GetRbtHome - returns value of RBT_HOME env variable
//(or HOME if RBT_HOME is undefined)
std::string GetRbtHome();
// GetCopyright - returns legalese statement
std::string GetCopyright();
// GetVersion - returns current library version
std::string GetVersion();
// GetBuildNum - returns current library build number
std::string GetBuild();
// GetProduct - returns library product name
std::string GetProduct();
// GetTime - returns current time and date as an RbtString
std::string GetTime();
// GetCurrentDirectory - returns current working directory
std::string GetCurrentDirectory();
//
////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////
// FILE/DIRECTORY HANDLING FUNCTIONS
//
// GetRbtDirName
// Returns the full path to a subdirectory in the rDock directory structure
//
// For example, if RBT_ROOT environment variable is ~dave/ribodev/molmod/ribodev
// then GetRbtDirName("data") would return ~dave/ribodev/molmod/ribodev/data/
//
// If RBT_ROOT is not set, then GetRbtDirName returns ./ irrespective of the
// subdirectory asked for. Thus the fall-back position is that parameter files
// are read from the current directory if RBT_ROOT is not defined.
std::string GetRbtDirName(const std::string &strSubdir = "");

// GetRbtFileName
// As GetRbtDirName but returns the full path to a file in the rDock directory
// structure
std::string GetRbtFileName(const std::string &strSubdir,
                           const std::string &strFile);

// GetFileType
// Returns the string following the last "." in the file name.
// e.g. GetFileType("receptor.psf") would return "psf"
std::string GetFileType(const std::string &strFile);

// GetDirList
// Returns a list of files in a directory (strDir) whose names begin with
// strFilePrefix (optional) and whose type is strFileType (optional, as returned
// by GetFileType)
RbtStringList GetDirList(const std::string &strDir,
                         const std::string &strFilePrefix = "",
                         const std::string &strFileType = "");
//
////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////
// CONVERSION ROUTINES
//
// Converts (comma)-delimited string of segment names to segment map
RbtSegmentMap ConvertStringToSegmentMap(const std::string &strSegments,
                                        const std::string &strDelimiter = ",");
// Converts segment map to (comma)-delimited string of segment names
std::string ConvertSegmentMapToString(const RbtSegmentMap &segmentMap,
                                      const std::string &strDelimiter = ",");

// Returns a segment map containing the members of map1 which are not in map2
// I know, should really be a template so as to be more universal...one day
// maybe. Or maybe there is already an STL algorithm for doing this.
RbtSegmentMap SegmentDiffMap(const RbtSegmentMap &map1,
                             const RbtSegmentMap &map2);

// DM 30 Mar 1999
// Converts (comma)-delimited string to string list (similar to
// ConvertStringToSegmentMap, but returns list not map)
RbtStringList
ConvertDelimitedStringToList(const std::string &strValues,
                             const std::string &strDelimiter = ",");
// Converts string list to (comma)-delimited string (inverse of above)
std::string ConvertListToDelimitedString(const RbtStringList &listOfValues,
                                         const std::string &strDelimiter = ",");

//
////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////
// I/O ROUTINES
//
// PrintStdHeader
// Print a standard header to an output stream (for log files etc)
// Contains copyright info, library version, date, time etc
// DM 19 Feb 1999 - include executable information
std::ostream &PrintStdHeader(std::ostream &s,
                             const std::string &strExecutable = "");

// Helper functions to read/write chars from iostreams
// Throws error if stream state is not Good() before and after the read/write
// It appears the STL std::ios_base exception throwing is not yet implemented
// at least on RedHat 6.1, so this is a temporary workaround (yeah right)
//
// DM 25 Sep 2000 - with MIPSPro CC compiler, streamsize is not defined,
// at least with the iostreams library we are using, so typedef it here
// it is not in gcc 3.2.1 either SJ
//#ifdef __sgi
typedef int streamsize;
//#endif

void WriteWithThrow(std::ostream &ostr, const char *p,
                    streamsize n) throw(RbtError);
void ReadWithThrow(std::istream &istr, char *p, streamsize n) throw(RbtError);

} // namespace Rbt

#endif //_RBT_H_
