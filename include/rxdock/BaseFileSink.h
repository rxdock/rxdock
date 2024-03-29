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

// Abstract base class for text file-based output. Knows how to Open a file,
// Write a text cache to the file, and Close the file.
// It is up to derived classes to provide a Render method to populate the cache
// with text records of the appropriate format.

#ifndef _RBTBASEFILESINK_H_
#define _RBTBASEFILESINK_H_

#include <fstream>

#include "rxdock/Config.h"

#include <nlohmann/json.hpp>

using json = nlohmann::json;

namespace rxdock {

class BaseFileSink {
public:
  ////////////////////////////////////////
  // Constructors/destructors
  // BaseFileSink(const char* fileName);
  BaseFileSink(const std::string &fileName);

  virtual ~BaseFileSink(); // Default destructor

  friend void to_json(json &j, const BaseFileSink &baseFileSink);
  friend void from_json(const json &j, BaseFileSink &baseFileSink);

  ////////////////////////////////////////
  // Public methods
  ////////////////
  std::string GetFileName() const { return m_strFileName; }
  // void SetFileName(const char* fileName);
  void SetFileName(const std::string &fileName);
  bool StatusOK() { return Status().isOK(); }
  Error Status();

  // PURE VIRTUAL - MUST BE OVERRIDDEN IN DERIVED CLASSES
  virtual void Render() = 0;

protected:
  ////////////////////////////////////////
  // Protected methods
  ///////////////////
  // Write the cache to the file
  // DM 11 Feb 1999 - add flag to allow the cache to be written without clearing
  // it
  void Write(bool bClearCache = true);
  // Add a complete line to the cache
  void AddLine(const std::string &fileRec);
  // Replace a complete line in the cache
  void ReplaceLine(const std::string &fileRec, unsigned int nRec);
  // Is cache empty
  bool isCacheEmpty() const { return m_lineRecs.empty(); }

  // DM 06 Apr 1999 - append attribute is for derived class use only
  bool GetAppend() const {
    return m_bAppend;
  } // Get append status (true=append, false=overwrite)
  void SetAppend(bool bAppend) {
    m_bAppend = bAppend;
  } // Set append status (true=append, false=overwrite)

private:
  ////////////////////////////////////////
  // Private methods
  /////////////////

  BaseFileSink();                     // Disable default constructor
  BaseFileSink(const BaseFileSink &); // Copy constructor disabled by default
  BaseFileSink &
  operator=(const BaseFileSink &); // Copy assignment disabled by default

  void Open(bool bAppend = false);
  void Close();
  void ClearCache();

protected:
  ////////////////////////////////////////
  // Protected data
  ////////////////

private:
  ////////////////////////////////////////
  // Private data
  //////////////
  std::vector<std::string> m_lineRecs;
  std::string m_strFileName;
  std::ofstream m_fileOut;
  bool m_bAppend; // If true, Write() appends to file rather than overwriting
};

void to_json(json &j, const BaseFileSink &baseFileSink);
void from_json(const json &j, BaseFileSink &baseFileSink);

// Useful typedefs
typedef SmartPtr<BaseFileSink> BaseFileSinkPtr; // Smart pointer

} // namespace rxdock

#endif //_RBTBASEFILESINK_H_
