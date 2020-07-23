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

// Rbt file-based exceptions

#ifndef _RBTFILEERROR_H_
#define _RBTFILEERROR_H_

#include "RbtError.h"

namespace rxdock {

const std::string IDS_FILE_ERROR = "RBT_FILE_ERROR";
const std::string IDS_FILE_READ_ERROR = "RBT_FILE_READ_ERROR";
const std::string IDS_FILE_WRITE_ERROR = "RBT_FILE_WRITE_ERROR";
const std::string IDS_FILE_PARSE_ERROR = "RBT_FILE_PARSE_ERROR";
const std::string IDS_FILE_MISSING_PARAMETER = "RBT_FILE_MISSING_PARAMETER";
const std::string IDS_DIR_NOACCESS = "RBT_DIR_NOACCESS";
const std::string IDS_NO_FILEINDIR = "RBT_NO_FILEINDIR";
const std::string IDS_NOENV = "RBT_NO_ENVIRONMENT";
const std::string IDS_STRINGTOOLONG = "RBT_SRINGTOOLONG";

// Unspecified file error
class RbtFileError : public RbtError {
public:
  RbtFileError(const std::string &strFile, int nLine,
               const std::string &strMessage = "")
      : RbtError(IDS_FILE_ERROR, strFile, nLine, strMessage) {}
  // Protected constructor to allow derived file error classes to set error name
protected:
  RbtFileError(const std::string &strName, const std::string &strFile,
               int nLine, const std::string &strMessage = "")
      : RbtError(strName, strFile, nLine, strMessage) {}
};

// File read error
class RbtFileReadError : public RbtFileError {
public:
  RbtFileReadError(const std::string &strFile, int nLine,
                   const std::string &strMessage = "")
      : RbtFileError(IDS_FILE_READ_ERROR, strFile, nLine, strMessage) {}
};

// File write error
class RbtFileWriteError : public RbtFileError {
public:
  RbtFileWriteError(const std::string &strFile, int nLine,
                    const std::string &strMessage = "")
      : RbtFileError(IDS_FILE_WRITE_ERROR, strFile, nLine, strMessage) {}
};

// File parse error
class RbtFileParseError : public RbtFileError {
public:
  RbtFileParseError(const std::string &strFile, int nLine,
                    const std::string &strMessage = "")
      : RbtFileError(IDS_FILE_PARSE_ERROR, strFile, nLine, strMessage) {}
};

// Missing parameter error
class RbtFileMissingParameter : public RbtFileError {
public:
  RbtFileMissingParameter(const std::string &strFile, int nLine,
                          const std::string &strMessage = "")
      : RbtFileError(IDS_FILE_MISSING_PARAMETER, strFile, nLine, strMessage) {}
};

// Exceptions for RbtDirectorySource and derived classes:
//
// Directory access error
class RbtDirIsNotAccessible : public RbtFileError {
public:
  RbtDirIsNotAccessible(const std::string &strFile, int nLine,
                        const std::string &strMessage = "")
      : RbtFileError(IDS_DIR_NOACCESS, strFile, nLine, strMessage) {}
};

// Directory contains no file with that extension
class RbtNoFileInDir : public RbtFileError {
public:
  RbtNoFileInDir(const std::string &strFile, int nLine,
                 const std::string &strMessage = "")
      : RbtFileError(IDS_NO_FILEINDIR, strFile, nLine, strMessage) {}
};

// Misc exceptions
//
// If path is missing
class RbtEnvNotDefined : public RbtFileError {
public:
  RbtEnvNotDefined(const std::string &strFile, int nLine,
                   const std::string &strMessage = "")
      : RbtFileError(IDS_NOENV, strFile, nLine, strMessage) {}
};

// If string is too long
class RbtStringTooLong : public RbtFileError {
public:
  RbtStringTooLong(const std::string &strFile, int nLine,
                   const std::string &strMessage = "")
      : RbtFileError(IDS_STRINGTOOLONG, strFile, nLine, strMessage) {}
};

} // namespace rxdock

#endif //_RBTFILEERROR_H_
