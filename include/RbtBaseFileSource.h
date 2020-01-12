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

// Abstract base class for file-based input

#ifndef _RBTBASEFILESOURCE_H_
#define _RBTBASEFILESOURCE_H_

#include <fstream>

#include "RbtConfig.h"

// useful typedefs
typedef std::string RbtFileRec;
typedef std::vector<RbtFileRec> RbtFileRecList;
typedef RbtFileRecList::iterator RbtFileRecListIter;

// Max line length expected in file
const int MAXLINELENGTH = 255;

class RbtBaseFileSource {
public:
  // Constructors
  // RbtBaseFileSource(const char* fileName);
  RbtBaseFileSource(const std::string &fileName);
  RbtBaseFileSource(const std::string &fileName,
                    const std::string &strRecDelim);

  // Default destructor
  virtual ~RbtBaseFileSource();

  // Public methods

  RBTDLL_EXPORT std::string GetFileName();
  // void SetFileName(const char* fileName);
  void SetFileName(const std::string &fileName);

  // Status and StatusOK parse the file to check for errors
  bool StatusOK();
  RBTDLL_EXPORT RbtError Status();

  // FileStatus and FileStatusOK just try and read the file
  RBTDLL_EXPORT bool FileStatusOK();
  RbtError FileStatus();

  // Multi-record methods
  bool isMultiRecordSupported() { return m_bMultiRec; }
  RBTDLL_EXPORT void NextRecord();
  void Rewind();
  RBTDLL_EXPORT std::size_t GetEstimatedNumRecords();

protected:
  //////////////////////////////////////////////////////
  // Protected member functions
  // PURE VIRTUAL - MUST BE OVERRIDDEN IN DERIVED CLASSES
  virtual void Parse() = 0;
  // chance to give false when record delimiter at the beginning
  void Read(bool aDelimiterAtEnd = true);
  //////////////////////////////////////////////////////

  // Protected data
protected:
  bool m_bParsedOK;          // For use by Parse
  RbtFileRecList m_lineRecs; // Allow direct access to cache from derived
                             // classes

  // Private member functions
private:
  RbtBaseFileSource(); // Disable default constructor
  RbtBaseFileSource(
      const RbtBaseFileSource &); // Copy constructor disabled by default
  RbtBaseFileSource &
  operator=(const RbtBaseFileSource &); // Copy assignment disabled by default
  void Open();
  void Close();
  void ClearCache();

  // Private data
private:
  std::string m_strFileName;
  std::size_t m_fileSize;
  bool m_bReadOK; // For use by Read
  std::size_t m_numReads;
  std::size_t m_bytesRead;
  std::ifstream m_fileIn;
  char *m_szBuf;    // Line buffer
  bool m_bFileOpen; // Keep track of whether we've opened the file or not
  bool m_bMultiRec; // Is file multi-record ?
  std::string m_strRecDelim; // Record delimiter
};

#endif //_RBTBASEFILESOURCE_H_
