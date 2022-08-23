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

#include "rxdock/BaseFileSource.h"
#include "rxdock/FileError.h"
#include <cstring>
#include <sys/stat.h>

#include <loguru.hpp>

using namespace rxdock;

// Constructors
// BaseFileSource::BaseFileSource(const char* fileName)
//{
//  m_strFileName = fileName;
//  ClearCache();
//  _RBTOBJECTCOUNTER_CONSTR_("BaseFileSource");
//}

BaseFileSource::BaseFileSource(const std::string &fileName)
    : m_numReads(0), m_bytesRead(0), m_bFileOpen(false), m_bMultiRec(false) {
  m_strFileName = fileName;
  struct stat fileStat;
  if (stat(fileName.c_str(), &fileStat) == 0) {
    m_fileSize = fileStat.st_size;
  }
  m_szBuf = new char[MAXLINELENGTH + 1]; // DM 24 Mar - allocate line buffer
  ClearCache();
  _RBTOBJECTCOUNTER_CONSTR_("BaseFileSource");
}

// Multi-record constructor
BaseFileSource::BaseFileSource(const std::string &fileName,
                               const std::string &strRecDelim)
    : m_numReads(0), m_bytesRead(0), m_bFileOpen(false), m_bMultiRec(true),
      m_strRecDelim(strRecDelim) {
  m_strFileName = fileName;
  struct stat fileStat;
  if (stat(fileName.c_str(), &fileStat) == 0) {
    m_fileSize = fileStat.st_size;
  }
  m_szBuf = new char[MAXLINELENGTH + 1]; // DM 24 Mar - allocate line buffer
  ClearCache();
  _RBTOBJECTCOUNTER_CONSTR_("BaseFileSource");
}

// Default destructor
BaseFileSource::~BaseFileSource() {
  Close();
  ClearCache();
  delete[] m_szBuf; // DM 24 Mar - delete line buffer
  _RBTOBJECTCOUNTER_DESTR_("BaseFileSource");
}

// Public methods
std::string BaseFileSource::GetFileName() { return m_strFileName; }

// void BaseFileSource::SetFileName(const char* fileName)
//{
//  Close();
//  ClearCache();
//  m_strFileName = fileName;
//}

void BaseFileSource::SetFileName(const std::string &fileName) {
  Close();
  ClearCache();
  m_strFileName = fileName;
}

// Status and StatusOK parse the file to check for errors
bool BaseFileSource::StatusOK() { return Status().isOK(); }

Error BaseFileSource::Status() {
  // Try parsing the file and see what we catch
  try {
    Parse();
    // If we get here then everything is fine
    return Error();
  }

  // Got an Error
  catch (Error &error) {
    return error;
  }
}

// FileStatus and FileStatusOK just try and read the file
bool BaseFileSource::FileStatusOK() { return FileStatus().isOK(); }

Error BaseFileSource::FileStatus() {
  // Try reading the file and see what we catch
  try {
    Read();
    // If we get here then everything is fine
    return Error();
  }

  // Got an Error
  catch (Error &error) {
    return error;
  }
}

// Multi-record methods

// Force the reading of the next record by clearing the cache
// Doesn't actually read the record
void BaseFileSource::NextRecord() {
  if (m_bMultiRec) {
    ClearCache();
  }
}

// Rewind the file back to the first record
void BaseFileSource::Rewind() {
  if (m_bMultiRec) {
    Close();
    ClearCache();
  }
}

// Estimate the number of records in the file from the file size
std::size_t BaseFileSource::GetEstimatedNumRecords() {
  if (m_bMultiRec) {
    if (m_numReads > 0 && m_bytesRead > 0) {
      double avgBytesRead =
          static_cast<double>(m_bytesRead) / static_cast<double>(m_numReads);
      double estNumRecords = static_cast<double>(m_fileSize) / avgBytesRead;
      return static_cast<std::size_t>(estNumRecords);
    } else {
      return 0;
    }
  }
  return 1;
}

// Protected functions

void BaseFileSource::Read(bool aDelimiterAtEnd) {
  // If we haven't already read the file, do it now
  std::size_t bytesLastRead = 0;
  if (!m_bReadOK) {
    if (aDelimiterAtEnd) {
      ClearCache();
      try {
        Open();
        // Multi-record read
        // Only read up to record delimiter (or end of file)
        // and leave file open for next record
        if (m_bMultiRec) {
          const char *cszRecDelim = m_strRecDelim.c_str();
          int n = strlen(cszRecDelim);
          while ((m_fileIn.getline(m_szBuf, MAXLINELENGTH)) &&
                 (strncmp(m_szBuf, cszRecDelim, n) != 0)) {
            LOG_F(1, "File line read is {}", m_szBuf);
            bytesLastRead += std::strlen(m_szBuf) + 1;
            // adding 1 byte per getline for newline (LF) character, no need
            // to check for CRLF as they are considered unsupported
            m_lineRecs.push_back(m_szBuf);
          }
        }
        // Single-record read
        // Read entire file and close immediately
        else {
          while (m_fileIn.getline(m_szBuf, MAXLINELENGTH)) {
            bytesLastRead += std::strlen(m_szBuf) + 1;
            m_lineRecs.push_back(m_szBuf);
          }
          Close();
        }
        // DM 25 Mar 1999 - check for end of file (i.e. no lines read)
        if (m_lineRecs.empty())
          throw FileReadError(_WHERE_,
                              "End of file/empty record in " + m_strFileName);
      }
      // Catch exceptions so we can close the file, then rethrow it
      catch (Error &error) {
        Close();
        throw;
      }
    } else { // ie MOL2 file records starts (not ends) with a certain pattern
      ClearCache();
      try {
        Open();
        // Multi-record read
        // Only read up to next record delimiter (or end of file)
        // and leave file open for next record
        if (m_bMultiRec) {
          const char *cszRecDelim = m_strRecDelim.c_str();
          int n = strlen(cszRecDelim);
          // skip to the header stuff until the first record
          // AND the first delimiter line
          while ((m_fileIn.getline(m_szBuf, MAXLINELENGTH)) &&
                 (strncmp(m_szBuf, cszRecDelim, n) != 0))
            ;
          while ((m_fileIn.getline(m_szBuf, MAXLINELENGTH)) &&
                 (strncmp(m_szBuf, cszRecDelim, n) != 0)) {
            LOG_F(1, "File line read is {}", m_szBuf);
            bytesLastRead += std::strlen(m_szBuf) + 1;
            m_lineRecs.push_back(m_szBuf);
          }
        }
        // Single-record read
        // Read entire file and close immediately
        else {
          while (m_fileIn.getline(m_szBuf, MAXLINELENGTH)) {
            bytesLastRead += std::strlen(m_szBuf) + 1;
            m_lineRecs.push_back(m_szBuf);
          }
          Close();
        }
        // DM 25 Mar 1999 - check for end of file (i.e. no lines read)
        if (m_lineRecs.empty())
          throw FileReadError(_WHERE_,
                              "End of file/empty record in " + m_strFileName);
      }
      // Catch exceptions so we can close the file, then rethrow it
      catch (Error &error) {
        Close();
        throw;
      }
    }
    // If we get to here, we read the file OK
    m_bReadOK = true;
    m_numReads++;
    m_bytesRead += bytesLastRead;
  }
}

// Private functions
void BaseFileSource::Open() {
  // DM 23 Mar 1999 - check if file is already open, to allow Open() to be
  // called redundantly
  if (!m_bFileOpen)
    m_fileIn.open(m_strFileName.c_str(), std::ios_base::in);

  // If file did not open, throw an error
  if (!m_fileIn)
    throw FileReadError(_WHERE_, "Error opening " + m_strFileName);
  else
    m_bFileOpen = true;
}

void BaseFileSource::Close() {
  m_fileIn.close();
  m_bFileOpen = false;
}

void BaseFileSource::ClearCache() {
  m_lineRecs.clear();  // Get rid of the previous file records
  m_bReadOK = false;   // Indicate the cache is invalid
  m_bParsedOK = false; // Tell the Parse function in derived classes that
  // it will have to reparse the file
}

void rxdock::to_json(json &j, const BaseFileSource &baseFileSource) {
  j = json{{"file-name", baseFileSource.m_strFileName},
           {"file-size", baseFileSource.m_fileSize},
           {"read", baseFileSource.m_bReadOK},
           {"num-reads", baseFileSource.m_numReads},
           {"bytes-read", baseFileSource.m_bytesRead},
           // fileIn skipped
           {"buffer", *baseFileSource.m_szBuf},
           {"file-open", baseFileSource.m_bFileOpen},
           {"multi-rec", baseFileSource.m_bMultiRec},
           {"rec-delim", baseFileSource.m_strRecDelim}};
}

void rxdock::from_json(const json &j, BaseFileSource &baseFileSource) {
  j.at("file-name").get_to(baseFileSource.m_strFileName);
  j.at("file-size").get_to(baseFileSource.m_fileSize);
  j.at("read").get_to(baseFileSource.m_bReadOK);
  j.at("num-reads").get_to(baseFileSource.m_numReads);
  j.at("bytes-read").get_to(baseFileSource.m_bytesRead);
  // fileIn skipped
  j.at("buffer").get_to(*baseFileSource.m_szBuf);
  j.at("file-open").get_to(baseFileSource.m_bFileOpen);
  j.at("multi-rec").get_to(baseFileSource.m_bMultiRec);
  j.at("rec-delim").get_to(baseFileSource.m_strRecDelim);
}