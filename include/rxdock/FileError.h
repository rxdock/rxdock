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

// rxdock file-based exceptions

#ifndef _RBTFILEERROR_H_
#define _RBTFILEERROR_H_

#include "rxdock/Error.h"

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
class FileError : public Error {
public:
  FileError(const std::string &strFile, int nLine,
            const std::string &strMessage = "")
      : Error(IDS_FILE_ERROR, strFile, nLine, strMessage) {}
  // Protected constructor to allow derived file error classes to set error name
protected:
  FileError(const std::string &strName, const std::string &strFile, int nLine,
            const std::string &strMessage = "")
      : Error(strName, strFile, nLine, strMessage) {}
};

// File read error
class FileReadError : public FileError {
public:
  FileReadError(const std::string &strFile, int nLine,
                const std::string &strMessage = "")
      : FileError(IDS_FILE_READ_ERROR, strFile, nLine, strMessage) {}
};

// File write error
class FileWriteError : public FileError {
public:
  FileWriteError(const std::string &strFile, int nLine,
                 const std::string &strMessage = "")
      : FileError(IDS_FILE_WRITE_ERROR, strFile, nLine, strMessage) {}
};

// File parse error
class FileParseError : public FileError {
public:
  FileParseError(const std::string &strFile, int nLine,
                 const std::string &strMessage = "")
      : FileError(IDS_FILE_PARSE_ERROR, strFile, nLine, strMessage) {}
};

// Missing parameter error
class FileMissingParameter : public FileError {
public:
  FileMissingParameter(const std::string &strFile, int nLine,
                       const std::string &strMessage = "")
      : FileError(IDS_FILE_MISSING_PARAMETER, strFile, nLine, strMessage) {}
};

// Exceptions for DirectorySource and derived classes:
//
// Directory access error
class DirIsNotAccessible : public FileError {
public:
  DirIsNotAccessible(const std::string &strFile, int nLine,
                     const std::string &strMessage = "")
      : FileError(IDS_DIR_NOACCESS, strFile, nLine, strMessage) {}
};

// Directory contains no file with that extension
class NoFileInDir : public FileError {
public:
  NoFileInDir(const std::string &strFile, int nLine,
              const std::string &strMessage = "")
      : FileError(IDS_NO_FILEINDIR, strFile, nLine, strMessage) {}
};

// Misc exceptions
//
// If path is missing
class EnvNotDefined : public FileError {
public:
  EnvNotDefined(const std::string &strFile, int nLine,
                const std::string &strMessage = "")
      : FileError(IDS_NOENV, strFile, nLine, strMessage) {}
};

// If string is too long
class StringTooLong : public FileError {
public:
  StringTooLong(const std::string &strFile, int nLine,
                const std::string &strMessage = "")
      : FileError(IDS_STRINGTOOLONG, strFile, nLine, strMessage) {}
};

} // namespace rxdock

#endif //_RBTFILEERROR_H_
