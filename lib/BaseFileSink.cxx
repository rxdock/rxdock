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

#include "rxdock/BaseFileSink.h"
#include "rxdock/FileError.h"

using namespace rxdock;

////////////////////////////////////////
// Constructors/destructors
// BaseFileSink::BaseFileSink(const char* fileName) :
//  m_strFileName(fileName)
//{
//  _RBTOBJECTCOUNTER_CONSTR_("BaseFileSink");
//}

BaseFileSink::BaseFileSink(const std::string &fileName)
    : m_strFileName(fileName), m_bAppend(false) {
  _RBTOBJECTCOUNTER_CONSTR_("BaseFileSink");
}

BaseFileSink::~BaseFileSink() {
  Write(); // Just in case there is anything in the cache
  _RBTOBJECTCOUNTER_DESTR_("BaseFileSink");
}

////////////////////////////////////////
// Public methods
////////////////
void BaseFileSink::SetFileName(const std::string &fileName) {
  Write(); // Just in case there is anything in the cache
  m_strFileName = fileName;
}

Error BaseFileSink::Status() {
  // For file sinks, all we can is try and open the file for writing and see
  // what we catch
  try {
    Open(true); // Open for append, so as not to overwrite the file
    Close();
    // If we get here then everything is fine
    return Error();
  }

  // Got an Error
  catch (Error &error) {
    Close();
    return error;
  }
}

////////////////////////////////////////
// Protected methods
///////////////////
// DM 11 Feb 1999 - add flag to allow the cache to be written without clearing
// it
void BaseFileSink::Write(bool bClearCache) {
  // Only write the file if there is anything in the cache
  if (isCacheEmpty())
    return;

  try {
    Open(m_bAppend); // DM 06 Apr 1999 - open for append or overwrite, depending
                     // on m_bAppend attribute
    for (std::vector<std::string>::const_iterator iter = m_lineRecs.begin();
         iter != m_lineRecs.end(); iter++) {
      // for some reason the << overload is screwed up in some sstream
      // implementations so it is worth to pay this "pointless" price in
      // conversion
      std::string delimited((*iter).c_str());
      m_fileOut << delimited << std::endl;
      // m_fileOut << *iter << std::endl;
    }
    Close();
    if (bClearCache)
      ClearCache(); // Clear the cache so we don't write the file again
  }

  // Got an Error
  catch (Error &error) {
    Close();
    if (bClearCache)
      ClearCache();
    throw;
  }
}

// Add a complete line to the cache
void BaseFileSink::AddLine(const std::string &fileRec) {
  m_lineRecs.push_back(fileRec);
}

// Replace a complete line in the cache
void BaseFileSink::ReplaceLine(const std::string &fileRec, unsigned int nRec) {
  if (nRec < m_lineRecs.size())
    m_lineRecs[nRec] = fileRec;
}

////////////////////////////////////////
// Private methods
/////////////////
void BaseFileSink::Open(bool bAppend) {
  std::ios_base::openmode openMode = std::ios_base::out;
  if (bAppend)
    openMode = openMode | std::ios_base::app;
  m_fileOut.open(m_strFileName.c_str(), openMode);
  if (!m_fileOut)
    throw FileWriteError(_WHERE_, "Error opening " + m_strFileName);
}

void BaseFileSink::Close() { m_fileOut.close(); }

void BaseFileSink::ClearCache() { m_lineRecs.clear(); }

void rxdock::to_json(json &j, const BaseFileSink &baseFileSink) {
  j = json{{"line-recs", baseFileSink.m_lineRecs},
           {"file-name", baseFileSink.m_strFileName},
           // fileOut skipped
           {"append", baseFileSink.m_bAppend}};
}

void rxdock::from_json(const json &j, BaseFileSink &baseFileSink) {
  j.at("line-recs").get_to(baseFileSink.m_lineRecs);
  j.at("file-name").get_to(baseFileSink.m_strFileName);
  // fileOut skipped
  j.at("append").get_to(baseFileSink.m_bAppend);
}
